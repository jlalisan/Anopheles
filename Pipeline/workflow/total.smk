# Imports for python code within snakemake.
import os
import shutil
import subprocess
import re
import random

# Config file with all parameters.
configfile: "config.yaml"

# Prefetches all the files from the config.
rule prefetching:
    priority: 10
    output:
        # Output is put on temp so after they are no longer needed they are deleted.
        temp("sra_files/{sample}/{sample}.sra")
    log:
        "sra_files/logs/{sample}.log"
    shell:
        # Prefetches the files and outputs them to sra_files.
        """
        (prefetch {wildcards.sample} -O sra_files) > {log} 2>&1 && touch {output}
        echo {wildcards}
        """

# Downloads the fastq files from the config.
rule fasterqdump:
    priority: 9
    input:
        "sra_files/{sample}/{sample}.sra"
    output:
        "fastq_files/log/{sample}.log"
    shell:
        """
        (fasterq-dump -O fastq_files/ {input}) >{output} 2>&1 && touch {output}
        """

# Build the reference index for the Bowtie2 process.
rule bowtieindex:
    priority: 10
    input:
        # Reference genome from the config.
        refgen=config['refgen']
    output:
        "Bowtie2/btbuild.log"
    log:
        "Bowtie2/log/Bowtie.log"
    benchmark:
        "Bowtie2/benchmarks/Bowtiebench.csv"
    shell:
        """
        (bowtie2-build {input.refgen} Bowtie2/ref_genome_btindex > {output}) >{log} 2>&1
        """

# Trims the files according to the config settings.
rule trimmomatic:
    priority: 8
    input:
        "fastq_files/log/{sample}.log",
        "Bowtie2/btbuild.log"
    output:
        "trimmomatic/logs/{sample}_trimmed.log"
    params:
        # Connection to the trimmomatic jar.
        jar=config['trimmomatic']['jar'],
        # Connection to the trimmomatic adapter file.
        adapter=config['trimmomatic']['adapter'],
        # Connection to the trimmomatic settings for paired and single files.
        paired_config=config['trimmomatic']['paired'],
        single_config=config['trimmomatic']['single']
    message:
        "Started read trimming for {wildcards.sample}!"
    shell:
        # Checks for the files if they belong to single or paired and trims them accordingly.
        """
        if test -f "fastq_files/{wildcards.sample}_1.fastq"; then
            (java -jar {params.jar} PE -phred33 fastq_files/{wildcards.sample}_1.fastq fastq_files/{wildcards.sample}_2.fastq trimmomatic/{wildcards.sample}_1.trim.fastq trimmomatic/{wildcards.sample}_1.untrim.fastq trimmomatic/{wildcards.sample}_2.trim.fastq trimmomatic/{wildcards.sample}_2.untrim.fastq ILLUMINACLIP:{params.adapter}{params.paired_config}) 2>> {output}
        fi
        if test -f "fastq_files/{wildcards.sample}.fastq"; then
            (java -jar {params.jar} SE -phred33 fastq_files/{wildcards.sample}.fastq trimmomatic/{wildcards.sample}.trim.fastq ILLUMINACLIP:{params.adapter}{params.single_config}) 2>> {output}
        fi
        touch {output}
        """

# Merges the remainder of the paired files for the Bowtie2 process
rule bowtiemerger:
    input:
        "trimmomatic/logs/{sample}_trimmed.log"
    output:
        "trimmomatic/logs/{sample}_merged.log"
    shell:
        # Ignores the single files but will make an empty log file for them.
        """
        if test -f "trimmomatic/{wildcards.sample}_1.untrim.fastq"; then
            cat trimmomatic/{wildcards.sample}_1.untrim.fastq trimmomatic/{wildcards.sample}_2.untrim.fastq > trimmomatic/{wildcards.sample}_merged.fastq 2>> {output}
        else
            touch {output}
        fi
        """

# Maps the reads for paired and single files.
rule bowtie2:
    input:
        "trimmomatic/logs/{sample}_merged.log",
        "trimmomatic/logs/{sample}_trimmed.log"
    output:
        "Bowtie2/log/{sample}_unmapped.log"
    params:
        btindex=config['btindex']
    shell:
        # Tests if the files are paired or unpaired and maps them accordingly.
        """
        if test -f "trimmomatic/{wildcards.sample}_1.trim.fastq"; then
            bowtie2 -p 20 --end-to-end -x {params.btindex} --fr -1 trimmomatic/{wildcards.sample}_1.trim.fastq -2 trimmomatic/{wildcards.sample}_2.trim.fastq -U trimmomatic/{wildcards.sample}_merged.fastq --al-conc Bowtie2/{wildcards.sample}_unmapped.fastq --al Bowtie2/{wildcards.sample}_unpaired.unmapped.fastq 1> Bowtie2/{wildcards.sample}.sam 2>> {output} 2>&1
        fi
        if test -f "trimmomatic/{wildcards.sample}.trim.fastq"; then
            bowtie2 -p 10 --end-to-end -x {params.btindex} trimmomatic/{wildcards.sample}.trim.fastq --al Bowtie2/{wildcards.sample}_single.mapped.fastq 1> Bowtie2/{wildcards.sample}.sam 2>> {output}
        fi
        touch {output}
        """

# Get contigs for blasting. using RAY denovo assembly.
rule denovo:
    input:
        "Bowtie2/log/{sample}_unmapped.log"
    output:
        "denovo/log/{sample}_u.log.file",
        "denovo/log/{sample}_log.file"
    shell:
        # For single and paired files uses MPIrun for the denovo process.
        """
        if test -f "Bowtie2/{wildcards.sample}_unpaired.unmapped.fastq"; then
            mpirun -n 2 Ray -s Bowtie2/{wildcards.sample}_unpaired.unmapped.fastq -o denovo/{wildcards.sample}.forblast 1>/dev/null 2>> {output[0]}
        fi

        if test -f "Bowtie2/{wildcards.sample}_unmapped.1.fastq"; then
            mpirun -n 2 Ray -p Bowtie2/{wildcards.sample}_unmapped.1.fastq Bowtie2/{wildcards.sample}_unmapped.2.fastq -o denovo/{wildcards.sample}_u.forblast 1>/dev/null 2>> {output[1]}
        fi

        if test -f "Bowtie2/{wildcards.sample}_single.mapped.fastq"; then
            mpirun -n 2 Ray -s Bowtie2/{wildcards.sample}_single.mapped.fastq -o denovo/{wildcards.sample}.forblast 1>/dev/null 2>> {output[1]}
        fi

        touch {output}
        """

# Gathers all the contigs together and moves this to the correct location.
rule getcontigs:
    input:
        "denovo/log/{sample}_u.log.file",
        "denovo/log/{sample}_log.file"
    output:
        "contigs/{sample}.Contigs.fasta"
    shell:
        """
        if test -f "denovo/{wildcards.sample}_u.forblast/Contigs.fasta"; then
            cat denovo/{wildcards.sample}.forblast/Contigs.fasta denovo/{wildcards.sample}_u.forblast/Contigs.fasta > {output}
        else
            cat denovo/{wildcards.sample}.forblast/Contigs.fasta > {output} && touch {output}
        fi
        """

# Blasts the contigs to check for any hits.
rule blasting:
    input:
        "contigs/{sample}.Contigs.fasta"
    output:
        "blastoutput/{sample}_matches.tsv"
    params:
        diamond=config['diamond']
    log: "blastoutput/log/{sample}_log.file"
    run:
        if os.stat(f"{input}").st_size > 1:
            subprocess.call(f"./diamond blastx -d  {params.diamond} -q {input[0]} --sensitive --threads 8 --quiet -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids skingdoms sscinames stitle -k 1 -o {output} 2>> {log}",shell=True)
            subprocess.call(f"echo {wildcards.sample}",shell=True)
        else:
            subprocess.call(f"touch {output}",shell=True)


# Uses the Keyword.py script to create Keywords.txt according to the config.
rule keywords:
    output:
        "Keywords.txt"
    script:
        "keyword.py"

# Matches the blast results to the Keywords and puts this into a hits file.
rule blastmatcher:
    input:
        "blastoutput/{sample}_matches.tsv",
        "Keywords.txt"
    output:
        "blastoutput/{sample}_ortho_matches.tsv",
        "blastoutput/accessions/{sample}_acc_hits.txt"
    run:
        if os.stat(f"{input[0]}").st_size > 1:
            subprocess.call(f"grep -f {input[1]} {input[0]} > {output[0]}",shell=True)
            subprocess.call("awk '{{print $2}}' " + f"{output[0]} | sort | uniq > {output[1]} && touch {output}",shell=True)
        else:
            subprocess.call(f"touch {output}",shell=True)

# Fetches the hits and puts these into a seperate folder.
rule efetcher:
    priority: 2
    input:
        "blastoutput/accessions/{sample}_acc_hits.txt"
    benchmark:
        "blastoutput/fetched/bench/{sample}.bench.txt"
    log:
        "blastoutput/fetched/log/{sample}.log"
    params:
        # Api key. It will work without it, but it will give a long confusing output stream
        api=config['api-key']
    run:
        # To remember the SRR files and the Accession names
        my_srr = []
        my_acc = []

        for files in os.listdir("blastoutput//accessions"):
            if files.startswith("SRR"):
                my_srr.append(files)

        for srr in my_srr:
            with open(f"blastoutput//accessions//{srr}") as srr_file:
                for acc in srr_file:
                    # Incase any accession is named Uniprot-PQA4241 or anything. This will filter the PQA4241 name.
                    newacc = re.findall("[A-Z]+?[0-9 A-Z.]+[0-9.]+?",acc)
                    for myacc in newacc:
                        my_acc.append(myacc)
                        # Searches and fetches the accessions in nucleotide form
                        seqfetch = f"esearch -db protein -query {myacc} -api_key {params.api} | efetch -format fasta_cds_na > blastoutput//fetched//{srr.split('_')[0]}_{myacc}.fasta 2>> {log}"
                        subprocess.call(seqfetch,shell=True)

# Merges all the accessions per SRR file.
rule merge_acc:
    input:
        "blastoutput/fetched/bench/{sample}.bench.txt"
    output:
        "blastoutput/fetched/rework/{sample}_accession.fasta"
    shell:
        """
        if ! [[ -z $(grep '[^[:space:]]' blastoutput/{wildcards.sample}_matches.tsv) ]] ; then
            cat blastoutput/fetched/{wildcards.sample}* > {output}
        else
            touch {output}
        fi
        """

rule clean_acc:
    input:
        "blastoutput/fetched/rework/{sample}_accession.fasta"
    output:
        "blastoutput/fetched/total/{sample}_accession.fasta"
    shell:
        """
        sed -i '/^$/d' {input} | cat {input} > {output}
        """

# Build the new Bowtie index.
rule Bowtie_index:
    priority: 5
    input:
        "blastoutput/fetched/total/{sample}_accession.fasta"
    output:
        "Geneious/Bowtie2/{sample}/btbuild.log"
    log:
        "Geneious/Bowtie2/{sample}/log/Bowtie.log"
    benchmark:
        "Geneious/Bowtie2/{sample}/benchmark/Bowtiebench.csv"
    shell:
        """
        (if [ -s {input[0]} ]; then
            bowtie2-build {input} Geneious/Bowtie2/{wildcards.sample}/ref_genome_btindex > {output} 2>&1
        else
            touch {output}
        fi
        ) >{log} 2>&1"""

# Maps the contigs against the references (accessions).
rule Bowtie2:
    priority: 4
    input:
        "contigs/{sample}.Contigs.fasta",
        "Geneious/Bowtie2/{sample}/btbuild.log"
    params:
        btindex="Geneious/Bowtie2/{sample}/ref_genome_btindex"
    output:
        "Geneious/Bowtie2/mapped/{sample}.sam"
    log:
        "Geneious/Bowtie2/log/{sample}_unmapped.log"
    message:
        "Mapping file {wildcards.sample}"
    benchmark:
        "Geneious/Bowtie2/benchmark/{sample}_bench.txt"
    shell:
        """
        if [ -s {input[1]} ]; then
            bowtie2 -x {params.btindex} -fa {input[0]} --no-unal -S {output} 2>> {log}
        else
            touch {output}
        fi
        """

# Adjusts the sam in case of duplicate contig names.
rule adjust_sam:
    priority: 3
    input:
        "Geneious/Bowtie2/mapped/{sample}.sam"
    output:
        "Geneious/Bowtie2/done/{sample}.sam"
    run:
        subprocess.call(f"touch {output}",shell=True)
        # Opens the input files one by one.
        myfile = open(f'{input}')
        contigs = []
        # Opens the file in read and write mode.
        with open(f'{output}',"r+") as adjusted_sam:
            for line in myfile:
                mycontig = ""
                # Checks if the line starts with the correct prefix before starting.
                if line.startswith("@SQ"):
                    # Finds all contig names like contig-0
                    contig = (re.findall("[a-z]*-[0-9]*",line))
                    # If the contig is not inside the list then it gets added and writen to the file.
                    if contig not in contigs:
                        contigs.append(contig)
                        adjusted_sam.write(line)
                    # If the contig is a duplicate then the number at the end of the contig is changed.
                    else:
                        for item in contig:
                            # The number is changed to a random number between 0 and 10000.
                            adjusted_sam.write(re.sub("[a-z]*-[0-9]*","".join(
                                str(item).split("-")[0]) + "-" + str(random.randint(0,10000)),line,count=0,flags=0))
                else:
                    adjusted_sam.write(line)

# Creates the consensus sequence.
rule create_consensus:
    priority: 4
    input:
        "blastoutput/fetched/total/{sample}_accession.fasta",
        "Geneious/Bowtie2/done/{sample}.sam"
    output:
        "Geneious/chimeras/{sample}_consensus.fasta"
    log:
        "Geneious/chimeras/{sample}.log"
    shell:
        # Sorts the files, and calls them before making a consensus sequence.
        """
        if [ -s {input[1]} ]; then
            (
            samtools view -bS {input[1]}| samtools sort -o Geneious/chimeras/{wildcards.sample}_bowtie.bam
            bcftools mpileup -Ou -f {input[0]} Geneious/chimeras/{wildcards.sample}_bowtie.bam | bcftools call -mv -Oz -o Geneious/chimeras/{wildcards.sample}_calls.vcf.gz && bcftools index Geneious/chimeras/{wildcards.sample}_calls.vcf.gz
            bcftools norm -f {input[0]} Geneious/chimeras/{wildcards.sample}_calls.vcf.gz -Ob -o Geneious/chimeras/{wildcards.sample}_calls.norm.bcf
            bcftools filter --IndelGap 5 Geneious/chimeras/{wildcards.sample}_calls.norm.bcf -Ob -o Geneious/chimeras/{wildcards.sample}_calls.norm.flt-indels.bcf
            cat {input[0]} | bcftools consensus Geneious/chimeras/{wildcards.sample}_calls.vcf.gz > {output} && touch {output}
            ) >{log} 2>&1
        else
            touch {output}
        fi
        """

# Adjust the consensus, so it has the SRR name in it.
rule adjust_consensus:
    priority: 2
    input:
        "Geneious/chimeras/{sample}_consensus.fasta"
    output:
        "Postprocess/chimeras/{sample}_consensus.fasta"
    run:
        subprocess.call(f"touch {output} ",shell=True)
        myfile = open(f'{input}')
        with open(f'{output}',"r+") as adjusted_consensus:
            for line in myfile:
                # If the line is an header then append the SRR name to it.
                if line.startswith(">"):
                    adjusted_consensus.write("".join(f">{wildcards.sample}" + "_" + line.strip(">")))
                else:
                    adjusted_consensus.write(line)


rule Reference_builder:
    input:
        "Postprocess/chimeras/{sample}_consensus.fasta"
    output:
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    shell:
        """
        if [ -s {input[0]} ]; then
            grep -A 1 "{wildcards.sample}" {input} --no-group-separator > {output} && touch {output}
        else
            touch {output}
        fi
        """

rule Bowtie2_index:
    input:
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/references/{sample}/btbuild.log"
    log:
        "Postprocess/references/{sample}/{sample}_log.file"
    shell:
        """
        if [ -s {input[0]} ]; then
            (bowtie2-build --threads 35 -f {input} Postprocess/references/{wildcards.sample}/ref_genome_btindex > {output}) > {log} && touch {output} 2>&1
        else
            touch {output}
        fi
        """


rule Bowtie2_post:
    input:
        "Bowtie2/log/{sample}_unmapped.log",
        "Postprocess/references/{sample}/btbuild.log"
    output:
        "Postprocess/Bowtie2/{sample}/{sample}.sam"
    params:
        btindex="Postprocess/references/{sample}/ref_genome_btindex"
    log:
        "Postprocess/Bowtie2/{sample}/map_log_file.txt"
    shell:
        """
        if [ -s Postprocess/references/{wildcards.sample}/btbuild.log ]; then
            if [ -s Bowtie2/{wildcards.sample}_single.mapped.fastq ]; then
                bowtie2 -p 35 --very-sensitive -x {params.btindex} Bowtie2/{wildcards.sample}_single.mapped.fastq >> {output} 2>> {log} 
            else
                bowtie2 -p 35 --very-sensitive -x {params.btindex} --fr -1 Bowtie2/{wildcards.sample}_unmapped.1.fastq -2 Bowtie2/{wildcards.sample}_unmapped.2.fastq -U Bowtie2/{wildcards.sample}_unpaired.unmapped.fastq 1> {output} 2>> {log}
            fi
        else
            touch {output}
        fi  
        """

rule sam_to_bam:
    input:
        "Postprocess/Bowtie2/{sample}/{sample}.sam"
    output:
        "Postprocess/Bowtie2/{sample}/{sample}_sort.bam",
        "Postprocess/Bowtie2/{sample}/{sample}_sort.REF_unmapped.bam",
        "Postprocess/Bowtie2/{sample}/{sample}.bam"
    log:
        "Postprocess/Bowtie2/{sample}/{sample}_samtobam.log"
    shell:
        """
        (if [ -s Postprocess/Bowtie2/{wildcards.sample}/{wildcards.sample}.sam ]; then
            samtools view --threads 20 -bSh -o {output[2]} {input}
            samtools sort --threads 20 {output[2]} -o {output[0]}
            bamtools split -in {output[0]} -reference

        else
            touch {output}
        fi
        ) >{log} 2>&1
        """

rule bamrename:
    input:
        "Postprocess/Bowtie2/{sample}/{sample}_sort.REF_unmapped.bam",
    output:
        touch("logs/bams/{sample}_log.file"),
        directory("Postprocess/Bowtie2/{sample}/accessions")
    run:
        srrname = []
        accname = []
        os.makedirs(f"Postprocess/Bowtie2/{wildcards.sample}/accessions",exist_ok=True)
        for files in os.listdir(f"Postprocess/Bowtie2/{wildcards.sample}"):
            if files.startswith(f"{wildcards.sample}_sort.REF_{wildcards.sample}"):
                srrname.append(files.split('_')[0])
                accname.append(re.findall("[A-Z 0-9]*[.][0-9]_1",files))

        for index in range(len(srrname)):
            srr = srrname[index]
            acc = accname[index]

            for files in os.listdir(f"Postprocess/Bowtie2/{wildcards.sample}"):
                if ''.join(acc) in files and ''.join(srr) in files:
                    os.renames(f"Postprocess/Bowtie2/{wildcards.sample}/{files}",f"Postprocess/Bowtie2/{wildcards.sample}/accessions/{''.join(srr)}_{''.join(acc)}.bam")

rule ref:
    input:
        "Postprocess/references/{sample}/{sample}_reference.fasta",
        "logs/bams/{sample}_log.file"
    output:
        "Postprocess/references/{sample}_refs.txt"
    shell:
        """
        if [ -s {input[0]} ]; then
            grep ">" {input[0]} | sed 's/.*_//' > {output}
            awk '/^>/ {{OUT=substr($0,1) ".fa"}}; OUT {{print >OUT}}' {input[0]} && touch {output}
        else
            touch {output}
        fi
        """

rule filerefrename:
    input:
        "Postprocess/references/{sample}_refs.txt"
    output:
        touch("logs/reference/{sample}_log.file"),
        directory("Postprocess/references/{sample}/accessions")
    params:
        workdir=config['workdir']
    run:
        srrnames = []
        accnames = []
        os.makedirs(f"Postprocess/references/{wildcards.sample}/accessions",exist_ok=True)
        for items in os.listdir(f"{params.workdir}"):
            if items.startswith(f">{wildcards.sample}") or items.startswith(f"{wildcards.sample}"):
                srrnames.append(re.findall("[A-Z]...[0-9]*_",items))
                accnames.append(re.findall("[A-Z]...[0-9]*[.][0-9]_[0-9]",items))

        for index in range(len(srrnames)):
            srr = srrnames[index]
            acc = accnames[index]

            for items in os.listdir(f"{params.workdir}"):
                if ''.join(acc) in items and ''.join(srr) in items:
                    print(f"{''.join(srr)}{''.join(acc)}.fa")
                    os.renames(items,f"{''.join(srr)}{''.join(acc)}.fa")
                    shutil.move(f"{params.workdir}/{''.join(srr)}{''.join(acc)}.fa",f"Postprocess/references/{wildcards.sample}/accessions/{''.join(srr)}{''.join(acc)}.fa")

"-----------------------------------------------------------------------------------"

rule process:
    input:
        "Postprocess/references/{sample}_refs.txt",
        "Postprocess/Bowtie2/{sample}/accessions",
        "Postprocess/references/{sample}/accessions",
        "logs/reference/{sample}_log.file"
    output:
        "Postprocess/processed/{sample}_log.file",
        directory("Postprocess/processed/{sample}")
    run:
        os.makedirs(f"Postprocess/processed/{wildcards.sample}",exist_ok=True)
        for files in os.listdir(f"{input[1]}"):
            if files.endswith(".bam"):
                subprocess.call(f"samtools index {input[1]}/{files}",shell=True)
                subprocess.call(f"bedtools genomecov -d -ibam {input[1]}/{files} > {output[1]}/{files.split('.bam')[0]}_long.cov",shell=True)
                subprocess.call(f"grep {wildcards.sample} {output[1]}/{files.split('.bam')[0]}_long.cov > {output[1]}/{files.split('.bam')[0]}.cov",shell=True)
        for files in os.listdir(f"{input[2]}"):
            if files.endswith(".fa"):
                subprocess.call(f"samtools faidx {input[2]}/{files}",shell=True)
                subprocess.call(f"./lofreq call -f {input[2]}/{files} -o {output[1]}/{files.split('.fa')[0]}.vcf {input[1]}/{files.split('.fa')[0]}.bam",shell=True)
        subprocess.call(f"touch {output[0]}",shell=True)

rule deletefiles:
    input:
        "Postprocess/processed/{sample}_log.file"
    output:
        "Postprocess/finished/{sample}_log.file"
    shell:
        """
        find fastq_files/ \( -name '{wildcards.sample}*.fastq' \) -exec rm "{{}}" \;
        find trimmomatic/ \( -name '{wildcards.sample}*.fastq' \) -exec rm "{{}}" \;
        find Bowtie2/ \( -name '{wildcards.sample}*.fastq' -o -name '{wildcards.sample}*.sam' \) -exec rm "{{}}" \;
        echo files for {wildcards.sample} have been deleted > {output}
        """
