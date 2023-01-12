# Imports for python code within snakemake.
import subprocess
import re
import random

# Config file with all parameters.
configfile: "Pipeline/config/config.yaml"

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
        ) >{log} 2>&1
        """

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
                # If the line is a header then append the SRR name to it.
                if line.startswith(">"):
                    adjusted_consensus.write("".join(f">{wildcards.sample}" + "_" + line.strip(">")))
                else:
                    adjusted_consensus.write(line)
