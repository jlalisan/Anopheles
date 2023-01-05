# Imports for python code within snakemake.
import os
import shutil
import subprocess
import re
import random

# Config file with all parameters.
configfile: "config.yaml"

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
        """(
        if [ -s Postprocess/Bowtie2/{wildcards.sample}/{wildcards.sample}.sam ]; then
            samtools view --threads 20 -bSh -o {output[2]} {input}
            samtools sort --threads 20 {output[2]} -o {output[0]}
            bamtools split -in {output[0]} -reference

        else
            touch {output}
        fi
        ) >{log} 2>&1"""

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
