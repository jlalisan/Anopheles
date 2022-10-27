#!/bin/bash

# User: Nadja Brait
# Group: Virus evolution and ecology (Asst.Prof. Sebastian Lequime)
# System: Linux slequime-HP-Z8-G4-Workstation 5.8.0-38-generic #43~20.04.1-Ubuntu SMP Tue Jan 12 16:39:47 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux
# Date: 2021/02
# Project: "Screening of Orthomyxovirus sequences in Culicinae mosquitoes"

####################################################################################################

				###WORKFLOW###

####################################################################################################

# This script describes a screening process for a target virus against transcriptomice samples from the Sequence Read Archive (SRA).

### SOFTWARE ### please change software paths to your own

export PATH="$PATH:C:/Users/hp/Desktop/Tools/ray-build/Ray-2.3.1"
export PATH="$PATH:C:/Users/hp/Desktop/Tools/sra_toolkit/bin"
export PATH="$PATH:C:/Users/hp/Desktop/Tools/ncbi-blast/bin"
export BLASTDB="C:/Users/hp/Desktop/data/NCBI-Genbank/nr"
TRIMMOMATIC_DIR="C:/Users/hp/Desktop/Tools/Trimmomatic/Trimmomatic-0.39"
### PROJECT ###

PROJECT_DIR="$HOME" #please change the project path to your own
R_DIR="$HOME/R/dataframe_generation" # please change R directory to your own

#creating directories in PROJECT_DIR
mkdir -p $PROJECT_DIR/{DATA/SEQ,ANALYSIS,RESULTS/{POSITIVE_HITS,BLAST_RESULTS},REFERENCES/{HOST,VIRUS,ADAPTERS,KEYWORDS}}
### PATHS ###

REF_DIR="$PROJECT_DIR/REFERENCES"
ADAPTER_DIR="$REF_DIR/ADAPTERS" # put your adapter files into this folder
VIRUSREF_DIR="$REF_DIR/VIRUS" # input virus reference sequence as reference.fa
HOSTREF_DIR="$REF_DIR/HOST/CULEX" # input host reference sequence
KEYWORDS_DIR="$REF_DIR/KEYWORDS" # input keywords (species name)
FASTA_DIR="$PROJECT_DIR/DATA/SEQ"
ANALYSIS_DIR="$PROJECT_DIR/ANALYSIS"
READS_DIR="$ANALYSIS_DIR/HOST_UNMAPPED"
SRA_DIR="$HOME/PROJECTS/SRA_DEPOSIT" # SRA output directory has to be defined during SRA toolkit configuration
BLAST_DIR="$STORAGE_DIR/BLAST_RESULTS"
RESULTS_DIR="$STORAGE_DIR/POSITIVE_HITS"
STORAGE_DIR="$HOME" # long-term storage for bigger files
### INPUT ###
mkdir $STORAGE_DIR/unmapped

INPUT="C:/Users/hp/Documents/GitHub/Anopheles/Pipeline/scripts/mytable.txt" # Accession list with your SRA samples one per lane
#################################################################################################################################

### SAMPLE DOWNLOAD ###

while IFS="" read -u 9 SRR || [ -n "$SRR" ] #IFS=input file seperator #9:can be any random number except 1(stdout) and 2(stderr)
do
date
	echo "##### $SRR screening #####" >> $STORAGE_DIR/log_file
#dowloading the subset with prefetch.2.9.1 from SRA toolkit
	cd $SRA_DIR # SRA files downloaded to directory I am currently in (defined during SRA toolkit configuration)
	prefetch --max-size 25000000 ${SRR} #change max size if some files are too big
#extraction of FASTQ from SRA-accessions
        fasterq-dump $SRA_DIR/${SRR}/${SRR}.sra -O $FASTA_DIR
	gzip $FASTA_DIR/*.fastq
#only if fasta files present, continue with workflow, otherwise start with next sample

        if [ ! -e $FASTA_DIR/${SRR}_1.fastq.gz ] || [ ! -e $FASTA_DIR/${SRR}_2.fastq.gz ] && [ ! -e $FASTA_DIR/${SRR}.fastq.gz ];then #either paired or single file
                echo "$SRR: error in sample: $(date) " | tee -a $STORAGE_DIR/log_file
		n=$(awk -v x="$SRR" '$0~x {print NR}' $INPUT)
        	echo "--------------------------> $n analysed: $(date) " | tee -a $STORAGE_DIR/log_file
                continue
        fi
#################################################################################################################################

### TRIMMING - QUALITY, ADAPTER AND LENGTH ###

# Trimmmomatic for paired end sequences
        if [ -e $FASTA_DIR/${SRR}_1.fastq.gz ] && [ -e $FASTA_DIR/${SRR}_2.fastq.gz ];then
                java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE -threads 35 -summary $ANALYSIS_DIR/trim_metadata/${SRR}_summary.txt \
		-quiet $FASTA_DIR/${SRR}_1.fastq.gz $FASTA_DIR/${SRR}_2.fastq.gz \
                $ANALYSIS_DIR/${SRR}_1.trim.fastq.gz $ANALYSIS_DIR/${SRR}_1un.trim.fastq.gz \
                $ANALYSIS_DIR/${SRR}_2.trim.fastq.gz $ANALYSIS_DIR/${SRR}_2un.trim.fastq.gz \
                ILLUMINACLIP:$ADAPTER_DIR/all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45
                echo "${SRR} was PE trimmed: $(date) " | tee -a $STORAGE_DIR/log_file
# Trimmomatic for single end sequences
        else
                java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar SE -threads 35 -summary $ANALYSIS_DIR/trim_metadata/${SRR}_summary.txt \
		-quiet $FASTA_DIR/${SRR}.fastq.gz $ANALYSIS_DIR/${SRR}_single.fastq.gz \
                ILLUMINACLIP:$ADAPTER_DIR/all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45
                echo "${SRR} was SE trimmed: $(date) " | tee -a $STORAGE_DIR/log_file
        fi
#################################################################################################################################

### BOWTIE2 MAPPING ### Screening for putative Orthomyxovirus sequences

	cd $ANALYSIS_DIR
# Quaranja Reference Indexing
	if [ ! -e $VIRUSREF_DIR/reference-index.1.bt2 ];then
                bowtie2-build --threads 35 -f $VIRUSREF_DIR/reference.fa $VIRUSREF_DIR/reference-index
        fi
# concatenate unpaired reads from paired end sequences
	if [ -e ${SRR}_1un.trim.fastq.gz ];then
        	cat ${SRR}_1un.trim.fastq.gz ${SRR}_2un.trim.fastq.gz > ${SRR}_merged.trim.fastq.gz
	fi
# gunzip all files in directory
	gunzip *fastq.gz # bowtie cannot use zipped files
	echo "##### $SRR: mapping metadata #####" | tee -a mapped_log.file
#MAPPING WITH BOWTIE2 for paired end reads
	if [ -e ${SRR}_1.trim.fastq ] && [ -e ${SRR}_2.trim.fastq ];then
		bowtie2 -p 35 --end-to-end -x $VIRUSREF_DIR/reference-index --fr -1 ${SRR}_1.trim.fastq -2 ${SRR}_2.trim.fastq \
		-U ${SRR}_merged.trim.fastq --al-conc ./${SRR}_%.mapped.fastq --al ./${SRR}_unpaired.mapped.fastq 1> ${SRR}.sam 2>> mapped_log.file
		echo "${SRR}: PE Virus mapping finished $(date) " | tee -a $STORAGE_DIR/log_file
#MAPPING WITH BOWTIE2 for single end reads
	else
		bowtie2 -p 35 --end-to-end -x $VIRUSREF_DIR/reference-index ${SRR}_single.fastq --al ./${SRR}_single.mapped.fastq 1> ${SRR}.sam 2>> mapped_log.file
		echo "${SRR}: SE Virus mapping finished $(date) " | tee -a $STORAGE_DIR/log_file
	fi
## Mapping to host file to exclude host sequences - continuing with unmapped sequences
# for my project I had to do it for all samples, also for the ones with no Orthomyxo hit for future research, in general however smarter to only use the ones that had a positive hit

# Aedes reference Indexing
        if [ ! -e $HOSTREF_DIR/reference-index.1.bt2 ];then
                bowtie2-build --threads 35 -f $HOSTREF_DIR/VectorBase_AaegyptiLVP_AalbopictusFPA_Transcripts.fasta $HOSTREF_DIR/reference-index
        fi
#MAPPING WITH BOWTIE2 for paired end reads
	echo "##### $SRR: mapping metadata #####" | tee -a unmapped_log.file
        if  [ -e ${SRR}_1.trim.fastq ] && [ -e ${SRR}_2.trim.fastq ];then
                bowtie2 -p 35 --local -x $HOSTREF_DIR/reference-index --fr -1 ${SRR}_1.trim.fastq -2 ${SRR}_2.trim.fastq \
		-U ${SRR}_merged.trim.fastq --un-conc ./${SRR}_%.unmapped.fastq --un ./${SRR}_unpaired.unmapped.fastq 1> ./${SRR}.sam 2>> unmapped_log.file
		echo "${SRR}: PE Culex mapping finished $(date) " | tee -a $STORAGE_DIR/log_file
#MAPPING WITH BOWTIE2 for single end reads
        else
                bowtie2 -p 35 --local -x $HOSTREF_DIR/reference-index ${SRR}_single.fastq --un ./${SRR}_single.unmapped.fastq 1> ${SRR}.sam 2>> unmapped_log.file
		echo "${SRR}: SE Culex mapping finished $(date) " | tee -a $STORAGE_DIR/log_file
	fi
#DETECT MAPPED READS

#MOVING TARGET READS TO STORAGE
	mv $ANALYSIS_DIR/*unmapped.fastq  $STORAGE_DIR/unmapped
	if [ ! -s ${SRR}_1.mapped.fastq ] || [ ! -s ${SRR}_2.mapped.fastq ] && [ ! -s ./${SRR}_single.mapped.fastq ];then
                echo "$SRR: No mapped target reads: $(date) " | tee -a $STORAGE_DIR/log_file
		echo "$SRR" >> $STORAGE_DIR/zeromapping_log_file
		n=$(awk -v x="$SRR" '$0~x {print NR}' $INPUT)
	        echo "--------------------------> $n samples analysed: $(date) " | tee -a $STORAGE_DIR/log_file
		rm -r $ANALYSIS_DIR/*.sam $FASTA_DIR/*.fastq.gz $FASTA_DIR/*.fastq $SRA_DIR/${SRR}
		continue
	fi
	rm -r $ANALYSIS_DIR/*.sam $FASTA_DIR/*.fastq.gz $SRA_DIR/${SRR}
 	echo "$SRR: Mapped succesfully to target: $(date) " | tee -a $STORAGE_DIR/log_file
#################################################################################################################################

### DE-NOVO ASSEMBLY WITH RAY ###

#assembly with paired-end reads:
	cd $STORAGE_DIR/unmapped
	mkdir $STORAGE_DIR/${SRR}.forblast
	echo "$SRR error report:" >>denovo_log_file
        if [ -e ${SRR}_1.unmapped.fastq ] && [ -e ${SRR}_2.unmapped.fastq ];then
               mpiexec -n 20 Ray -p ${SRR}_1.unmapped.fastq ${SRR}_2.unmapped.fastq -o $STORAGE_DIR/${SRR}.forblast/paired 1>/dev/null 2>>denovo_log_file
               mpiexec -n 20 Ray -s ${SRR}_unpaired.unmapped.fastq -o $STORAGE_DIR/${SRR}.forblast/unpaired 1>/dev/null 2>>denovo_log_file
		if [ -e $STORAGE_DIR/${SRR}.forblast/unpaired/Contigs.fasta ];then
			sed 's/ /_U /' $STORAGE_DIR/${SRR}.forblast/unpaired/Contigs.fasta > $STORAGE_DIR/${SRR}.forblast/unpaired/ContigsU.fasta
			cat $STORAGE_DIR/${SRR}.forblast/unpaired/ContigsU.fasta $STORAGE_DIR/${SRR}.forblast/paired/Contigs.fasta > $STORAGE_DIR/${SRR}.forblast/Contigs.fasta
			echo "$SRR: PE denovo assembly: $(date)" | tee -a $STORAGE_DIR/log_file
		fi
##assembly with single-end reads:
        else
                mpiexec -n 20 Ray -s ${SRR}_single.unmapped.fastq -o $STORAGE_DIR/${SRR}.forblast/${SRR} 1>/dev/null 2>>denovo_log_file
		mv $STORAGE_DIR/${SRR}.forblast/${SRR}/Contigs.fasta $STORAGE_DIR/${SRR}.forblast
		echo "$SRR: SE denovo assembly: $(date) " | tee -a $STORAGE_DIR/log_file
        fi
#################################################################################################################################

### BLAST AND TARGET FILTERING ###

	if [ -e $STORAGE_DIR/${SRR}.forblast/Contigs.fasta ];then
                ~/PROJECTS/SOFTWARE/diamond blastx -d /media/nbrait/Data/diamond_blast_nr/nr_diamond \
                -q $STORAGE_DIR/${SRR}.forblast/Contigs.fasta --sensitive --quiet  \
                -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids skingdoms sscinames stitle \
                -k 1 -o $STORAGE_DIR/blast_diamond/${SRR}_matches.tsv
                echo "$SRR screened"
        else
                echo "$SRR not available"
        fi

	grep -f $KEYWORDS_DIR/keywords2 $STORAGE_DIR/blast_diamond/${SRR}_matches.tsv > $STORAGE_DIR/${SRR}_ortho_matches.tsv
        if [ -s $STORAGE_DIR/${SRR}_ortho_matches.tsv ];then
        	echo "$SRR: Positive blastx hit: $(date) " | tee -a $STORAGE_DIR/log_file $STORAGE_DIR/log_file
		awk '{print $2}' $STORAGE_DIR/${SRR}_ortho_matches.tsv | sort | uniq > $STORAGE_DIR/${SRR}_acc_hits.txt # for protein
	        cd $STORAGE_DIR
	        mkdir ${SRR}
	        while IFS="" read -u 8 Acc
                do
                        esearch -db protein -query $Acc | efetch -format fasta_cds_na > ${SRR}/${SRR}_${Acc}.fasta
                done 8< "${SRR}_acc_hits.txt"

        else
                echo "$SRR: no blastx hit: $(date) " | tee -a $STORAGE_DIR/log_file $STORAGE_DIR/false-positives_log_file
		mv ../${SRR}_contigs.fasta $STORAGE_DIR/
        fi
	rm  -r $ANALYSIS_DIR/*.fastq

#progress:
	n=$(awk -v x="$SRR" '$0~x {print NR}' $INPUT)
	echo "--------------------------> $n analysed: $(date) " | tee -a $STORAGE_DIR/log_file

done 9< "$INPUT"
################################################################################################################################
date
# THE END
