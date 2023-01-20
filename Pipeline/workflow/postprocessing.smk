# Imports for the Python code.
import os
import re
import shutil
import subprocess

# Imports the configuration file.
configfile: "Pipeline/config/config.yaml"

# Makes running the post-processing available from the terminal. Runs entire script
rule Do_postprocessing:
    input:
        expand("Postprocess/finished/{sample}_log.file", sample=config['samples'])

# Builds the new set of references.
rule Reference_builder:
    input:
        "Postprocess/logs/{sample}.log.file"
    output:
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    run:
        # Only looks for files in the build directory in case of the Geneious approach.
        for files in os.listdir(f"Postprocess/chimeras/"):
            if wildcards.sample in files:
                subprocess.call(f"grep -A 1 {wildcards.sample} Postprocess/chimeras/{files} --no-group-separator > {output} && touch {output}",shell=True)

# Build the Bowtie2 index for the post-processing.
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

# Runs the Bowtie2 mapping for the post-processing.
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

# Converts the Bowtie2 SAM files to BAM files.
rule Sam_to_bam:
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
        if [ -s {input[0]} ]; then
            (samtools view --threads 20 -bSh -o {output[2]} {input}
		    samtools sort --threads 20 {output[2]} -o {output[0]}
		    bamtools split -in {output[0]} -reference && touch {output}
            ) >{log} 2>&1
        else
            touch {output}
        fi
        """
# Renames the long file names the BAM files gained.
rule Bam_rename:
    input:
        "Postprocess/Bowtie2/{sample}/{sample}_sort.REF_unmapped.bam"
    output:
        touch("logs/bams/{sample}_log.file"),
        directory("Postprocess/Bowtie2/{sample}/accessions")
    run:
        # Storage for the SRR names and Accession names.
        srrname = []
        accname = []
        # Makes directory for files to go to.
        os.makedirs(f"Postprocess/Bowtie2/{wildcards.sample}/accessions",exist_ok=True)
        for files in os.listdir(f"Postprocess/Bowtie2/{wildcards.sample}"):
            # Only grab files with the correct name.
            if files.startswith(f"{wildcards.sample}_sort.REF_{wildcards.sample}"):
                srrname.append(files.split('_')[0])
                # Find the accession names and take these.
                accname.append(re.findall("[A-Z 0-9]*[.][0-9]_1",files))

        # Get ready to use the SRR and Accessions one by one.
        for index in range(len(srrname)):
            srr = srrname[index]
            acc = accname[index]
            # Renames any file to SRR_Accession.bam.
            for files in os.listdir(f"Postprocess/Bowtie2/{wildcards.sample}"):
                if ''.join(acc) in files and ''.join(srr) in files:
                    os.renames(f"Postprocess/Bowtie2/{wildcards.sample}/{files}",f"Postprocess/Bowtie2/{wildcards.sample}/accessions/{''.join(srr)}_{''.join(acc)}.bam")

# Creates separate references
rule Ref:
    input:
        "Postprocess/references/{sample}/{sample}_reference.fasta",
        "logs/bams/{sample}_log.file"
    output:
        "Postprocess/references/{sample}_refs.txt"
    shell:
        """
        grep ">" {input[0]} | sed 's/.*_//' > {output}
		awk '/^>/ {{OUT=substr($0,1) ".fa"}}; OUT {{print >OUT}}' {input[0]} && touch {output}
        """

# Renames the references after creation.
rule Ref_rename:
    input:
        "Postprocess/references/{sample}_refs.txt"
    output:
        touch("logs/reference/{sample}_log.file"),
        directory("Postprocess/references/{sample}/accessions")
    params:
        workdir=config['workdir']
    run:
        # Saves the SRR names and accessions.
        srrnames = []
        accnames = []
        # Similar to the Bamrename the files need proper names.
        os.makedirs(f"Postprocess/references/{wildcards.sample}/accessions",exist_ok=True)
        for items in os.listdir(f"{params.workdir}"):
            if items.startswith(f">{wildcards.sample}") or items.startswith(f"{wildcards.sample}"):
                # Finds the Srr file names and adds them to the list.
                srrnames.append(re.findall("[A-Z]...[0-9]*_",items))
                # Finds the accession names and adds them to the list.
                accnames.append(re.findall("[A-Z]...[0-9]*[.][0-9]_[0-9]",items))


        for index in range(len(srrnames)):
            srr = srrnames[index]
            acc = accnames[index]

            for items in os.listdir(f"{params.workdir}"):
                if ''.join(acc) in items and ''.join(srr) in items:
                    # Renames and moves the files for later usage.
                    os.renames(items,f"{''.join(srr)}{''.join(acc)}.fa")
                    shutil.move(f"{params.workdir}/{''.join(srr)}{''.join(acc)}.fa",f"Postprocess/references/{wildcards.sample}/accessions/{''.join(srr)}{''.join(acc)}.fa")

# Processes the files into VCF files
rule Process:
    input:
        "Postprocess/references/{sample}_refs.txt",
        "Postprocess/Bowtie2/{sample}/accessions",
        "Postprocess/references/{sample}/accessions",
        "logs/reference/{sample}_log.file"
    output:
        "Postprocess/processed/{sample}_log.file",
        directory("Postprocess/processed/{sample}")
    run:
        # Makes output directory.
        os.makedirs(f"Postprocess/processed/{wildcards.sample}",exist_ok=True)
        for files in os.listdir(f"{input[1]}"):
            if files.endswith(".bam"):
                # Indexes the files
                subprocess.call(f"samtools index {input[1]}/{files}",shell=True)
                # Creates the long coverage files to later be COV files.
                subprocess.call(f"bedtools genomecov -d -ibam {input[1]}/{files} > {output[1]}/{files.split('.bam')[0]}_long.cov",shell=True)
                # Grabs the coverage and makes COV files.
                subprocess.call(f"grep {wildcards.sample} {output[1]}/{files.split('.bam')[0]}_long.cov > {output[1]}/{files.split('.bam')[0]}.cov",shell=True)
        for files in os.listdir(f"{input[2]}"):
            if files.endswith(".fa"):
                # Indexes the references.
                subprocess.call(f"samtools faidx {input[2]}/{files}",shell=True)
                # Calls on Lofreq to build the VCF.
                subprocess.call(f"./lofreq call -f {input[2]}/{files} -o {output[1]}/{files.split('.fa')[0]}.vcf {input[1]}/{files.split('.fa')[0]}.bam",shell=True)
        subprocess.call(f"touch {output[0]}",shell=True)

# Deletes all files no longer used to avoid memory issues.
rule Delete_files:
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

# Makes Dag diagram on success
onsuccess:
    with open("postprocess.txt","w") as f:
        f.writelines(str(workflow.persistence.dag))
    shell("cat postprocess_dag.txt | dot -Tpng > postprocess.png")
