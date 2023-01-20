# Imports for python code within snakemake.
import os
import subprocess
import re

# Config file with all parameters.
configfile: "Pipeline/config/config.yaml"

rule Do_preprocess:
    input:
        expand("Postprocess/logs/{sample}.log.file", sample=config['samples'])

# Prefetches all the files from the config.
rule Prefetching:
    priority: 10
    output:
        # Output is put on temp, so after they are no longer needed they are deleted.
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
rule Fasterqdump:
    priority: 9
    input:
        "sra_files/{sample}/{sample}.sra"
    output:
        "fastq_files/log/{sample}.log"
    shell:
        """(
        fasterq-dump -O fastq_files/ {input} 
        ) >{output} 2>&1 && touch {output}"""

# Build the reference index for the Bowtie2 process.
rule Bowtieindex:
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
rule Trimmomatic:
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
rule Bowtie_merger:
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
rule Bowtie2:
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
            bowtie2 -p 20 --local --end-to-end -x {params.btindex} --fr -1 trimmomatic/{wildcards.sample}_1.trim.fastq -2 trimmomatic/{wildcards.sample}_2.trim.fastq -U trimmomatic/{wildcards.sample}_merged.fastq --al-conc Bowtie2/{wildcards.sample}_unmapped.fastq --al Bowtie2/{wildcards.sample}_unpaired.unmapped.fastq 1> Bowtie2/{wildcards.sample}.sam 2>> {output} 2>&1
        fi
        if test -f "trimmomatic/{wildcards.sample}.trim.fastq"; then
            bowtie2 -p 10 --local --end-to-end -x {params.btindex} trimmomatic/{wildcards.sample}.trim.fastq --al Bowtie2/{wildcards.sample}_single.mapped.fastq 1> Bowtie2/{wildcards.sample}.sam 2>> {output}
        fi
        touch {output}
        """

# Get contigs for blasting. using RAY denovo assembly.
rule Denovo:
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
rule Get_Contigs:
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
rule Blasting:
    input:
        "contigs/{sample}.Contigs.fasta"
    output:
        "blastoutput/{sample}_matches.tsv"
    params:
        diamond=config['diamond']
    log: "blastoutput/log/{sample}_log.file"
    run:
        if os.stat(f"{input}").st_size > 1:
            # Calls on diamond to do the blasting with the database.
            subprocess.call(f"./diamond blastx -d  {params.diamond} -q {input[0]} --sensitive --threads 8 --quiet -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids skingdoms sscinames stitle -k 1 -o {output} 2>> {log}",shell=True)
        else:
            subprocess.call(f"touch {output}",shell=True)


# Uses the Keyword.py script to create Keywords.txt according to the config.
rule Keywords:
    output:
        "Keywords.txt"
    script:
        "../scripts/keyword.py"

# Matches the blast results to the Keywords and puts this into a hits file.
rule Blast_matcher:
    input:
        "blastoutput/{sample}_matches.tsv",
        "Keywords.txt"
    output:
        "blastoutput/{sample}_ortho_matches.tsv",
        "blastoutput/accessions/{sample}_acc_hits.txt"
    run:
        if os.stat(f"{input[0]}").st_size > 1:
            # Grabs the correct keywords.
            subprocess.call(f"grep -f {input[1]} {input[0]} > {output[0]}",shell=True)
            # Matches these keys to the output and only matches get moved along.
            subprocess.call("awk '{{print $2}}' " + f"{output[0]} | sort | uniq > {output[1]} && touch {output}",shell=True)
        else:
            subprocess.call(f"touch {output}",shell=True)

# Fetches the hits and puts these into a separate folder.
rule Efetcher:
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
rule Merge_acc:
    input:
        "blastoutput/fetched/bench/{sample}.bench.txt"
    output:
        "blastoutput/fetched/rework/{sample}_accession.fasta"
    shell:
        """
        if test -f "blastoutput/{wildcards.sample}_ortho_matches.tsv"; then
            if ! [[ -z $(grep '[^[:space:]]' blastoutput/{wildcards.sample}_matches.tsv) ]] ; then
                cat blastoutput/fetched/{wildcards.sample}* > {output}
            else
                touch {output}
            fi
        else
            touch {output}
        fi
        """

# Removes any white space since SAMtools later on does not handle that well.
rule Clean_acc:
    input:
        "blastoutput/fetched/rework/{sample}_accession.fasta"
    output:
        "blastoutput/fetched/total/{sample}_accession.fasta"
    run:
        # Creates directory
        subprocess.call(f"touch {output} ",shell=True)
        myfile = open(f'{input}')
        with open(f'{output}',"r+") as accession:
            for line in myfile:
                # Sometimes the > is missing this is required.
                if line.startswith("lcl"):
                    accession.write("".join(f">" + line.strip()))
                else:
                    accession.write(line)
        # Cleans the file.
        subprocess.call(f"sed -i '/^$/d' {output}",shell=True)

# Creates directory for Geneious and cleans environment.
rule Clean_env:
    input:
        "blastoutput/fetched/total/{sample}_accession.fasta"
    output:
        "Postprocess/logs/{sample}.log.file"
    log:
        "Postprocess/logs/{sample}.log.file"
    shell:
        """
        if ! [ -s Postprocess/chimeras/ ]; then
            mkdir Postprocess/chimeras/
        else
            echo "Making correct directory"
        fi

        find fastq_files/ \( -name '{wildcards.sample}*.fastq' \) -exec rm "{{}}" \;
        find trimmomatic/ \( -name '{wildcards.sample}*.fastq' \) -exec rm "{{}}" \;
         >{log} 2>&1
        """

# Makes Dag diagram on success
onsuccess:
    with open("preprocess.txt","w") as f:
        f.writelines(str(workflow.persistence.dag))
    shell("cat preprocess_dag.txt | dot -Tpng > preprocess.png")