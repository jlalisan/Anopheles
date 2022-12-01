import subprocess
import random
import re

configfile: "config.yaml"

rule Bowtie_index_single:
    priority: 5
    input:
        "blastoutput/fetched/single/{sample}_accession.fasta",
    output:
        "Geneious/Bowtie2/single/{sample}/btbuild.log",
    log:
        "Geneious/Bowtie2/single/{sample}/log/Bowtie.log"
    benchmark:
        "Geneious/Bowtie2/single/{sample}/benchmark/Bowtiebench.csv"
    shell:
        """
        if [ -s {input[0]} ]; then
            (bowtie2-build {input} Geneious/Bowtie2/single/{wildcards.sample}/ref_genome_btindex > {output}) >{log} 2>&1
            echo {wildcards}
        else
            touch {output}
        fi
        """

rule Bowtie_index_paired:
    priority: 5
    input:
        "blastoutput/fetched/paired/{sample}_accession.fasta",
    output:
        "Geneious/Bowtie2/paired/{sample}/btbuild.log",
    log:
        "Geneious/Bowtie2/paired/{sample}/log/Bowtie.log"
    benchmark:
        "Geneious/Bowtie2/paired/{sample}/benchmark/Bowtiebench.csv"
    shell:
        """
        if [ -s {input[0]} ]; then
            (bowtie2-build {input} Geneious/Bowtie2/paired/{wildcards.sample}/ref_genome_btindex > {output}
            ) >{log} 2>&1
        else
            touch {output}
        fi
        """

rule Bowtie2_single:
    priority: 4
    input:
        "contigs/single/{sample}.Contigs.fasta",
        "Geneious/Bowtie2/single/{sample}/btbuild.log"
    params:
        btindex="Geneious/Bowtie2/single/{sample}/ref_genome_btindex"
    output:
        "Geneious/Bowtie2/single/mapped/{sample}.sam"
    log:
        "Geneious/Bowtie2/single/log/{sample}_unmapped.log"
    message:
        "Mapping single file {wildcards.sample}"
    benchmark:
        "Geneious/Bowtie2/single/benchmark/{sample}_bench.txt"
    shell:
        """
        if [ -s {input[0]} ]; then
            bowtie2 -x {params.btindex} -fa {input[0]} -S {output} 2>> {log} && touch {output}
        else
            touch {output}
        fi
        """

rule Bowtie2_paired:
    priority: 4
    input:
        "contigs/paired/{sample}.Contigs.fasta",
        "Geneious/Bowtie2/paired/{sample}/btbuild.log"
    params:
        btindex="Geneious/Bowtie2/paired/{sample}/ref_genome_btindex"
    output:
        "Geneious/Bowtie2/paired/mapped/{sample}.sam",
    log:
        "Geneious/Bowtie2/paired/log/{sample}_unmapped.log"
    message:
        "Mapping paired file {wildcards.sample}"
    benchmark:
        "Geneious/Bowtie2/paired/benchmark/{sample}_bench.txt"
    shell:
        """
        if [ -s {input[0]} ]; then
            bowtie2 -x {params.btindex} -fa {input[0]} -S {output} 2>> {log} && touch {output}
        else
            touch {output}
        fi
        """

rule adjust_sam_paired:
    priority: 3
    input:
        "Geneious/Bowtie2/paired/mapped/{sample}.sam"
    output:
        "Geneious/Bowtie2/paired/done/{sample}.sam"
    run:
        subprocess.call(f"touch {output}", shell=True)
        myfile = open(f'{input}')
        contigs = []
        with open(f'{output}',"r+") as newtest:
            for line in myfile:
                mycontig = ""
                if line.startswith("@SQ"):
                    contig = (re.findall("[a-z]*-[0-9]*",line))
                    if contig not in contigs:
                        contigs.append(contig)
                        newtest.write(line)
                    else:
                        for item in contig:
                            newtest.write(re.sub("[a-z]*-[0-9]*","".join(
                                str(item).split("-")[0]) + "-" + str(random.randint(0,10000)),line,count=0,flags=0))
                else:
                    newtest.write(line)

rule adjust_sam_single:
    priority: 3
    input:
        "Geneious/Bowtie2/single/mapped/{sample}.sam"
    output:
        "Geneious/Bowtie2/single/done/{sample}.sam"
    run:
        subprocess.call(f"touch {output}", shell=True)
        myfile = open(f'{input}')
        contigs = []
        with open(f'{output}',"r+") as newtest:
            for line in myfile:
                mycontig = ""
                if line.startswith("@SQ"):
                    contig = (re.findall("[a-z]*-[0-9]*",line))
                    if contig not in contigs:
                        contigs.append(contig)
                        newtest.write(line)
                    else:
                        for item in contig:
                            newtest.write(re.sub("[a-z]*-[0-9]*","".join(
                                str(item).split("-")[0]) + "-" + str(random.randint(0,10000)),line,count=0,flags=0))
                else:
                    newtest.write(line)

rule create_consencus_paired:
    priority: 4
    input:
        "blastoutput/fetched/paired/{sample}_accession.fasta",
        "Geneious/Bowtie2/paired/done/{sample}.sam"
    output:
        "Geneious/chimeras/paired/{sample}_concensus.fasta"
    log:
        "Geneious/chimeras/paired/{sample}.log"
    shell:
        """
        if [ -s {input[0]} ]; then
            (
            samtools view -bS {input[1]}| samtools sort -o Geneious/chimeras/paired/{wildcards.sample}_bowtie.bam
            bcftools mpileup -Ou -f {input[0]} Geneious/chimeras/paired/{wildcards.sample}_bowtie.bam | bcftools call -mv -Oz -o Geneious/chimeras/paired/{wildcards.sample}_calls.vcf.gz && bcftools index Geneious/chimeras/paired/{wildcards.sample}_calls.vcf.gz
            bcftools norm -f {input[0]} Geneious/chimeras/paired/{wildcards.sample}_calls.vcf.gz -Ob -o Geneious/chimeras/paired/{wildcards.sample}_calls.norm.bcf
            bcftools filter --IndelGap 5 Geneious/chimeras/paired/{wildcards.sample}_calls.norm.bcf -Ob -o Geneious/chimeras/paired/{wildcards.sample}_calls.norm.flt-indels.bcf
            cat {input[0]} | bcftools consensus Geneious/chimeras/paired/{wildcards.sample}_calls.vcf.gz > {output} && touch {output}
            ) >{log} 2>&1
        else
            touch {output}
        fi
        """


rule create_concensus_single:
    priority: 4
    input:
        "blastoutput/fetched/single/{sample}_accession.fasta",
        "Geneious/Bowtie2/single/done/{sample}.sam"
    output:
        "Geneious/chimeras/single/{sample}_concensus.fasta"
    shell:
        """
        if [ -s {input[0]} ]; then
            samtools view -bS {input[1]}| samtools sort -o Geneious/chimeras/single/{wildcards.sample}_bowtie.bam
            bcftools mpileup -Ou -f {input[0]} Geneious/chimeras/single/{wildcards.sample}_bowtie.bam | bcftools call -mv -Oz -o Geneious/chimeras/single/{wildcards.sample}_calls.vcf.gz && bcftools index Geneious/chimeras/single/{wildcards.sample}_calls.vcf.gz
            bcftools norm -f {input[0]} Geneious/chimeras/single/{wildcards.sample}_calls.vcf.gz -Ob -o Geneious/chimeras/single/{wildcards.sample}_calls.norm.bcf
            bcftools filter --IndelGap 5 Geneious/chimeras/single/{wildcards.sample}_calls.norm.bcf -Ob -o Geneious/chimeras/single/{wildcards.sample}_calls.norm.flt-indels.bcf
            cat {input[0]} | bcftools consensus Geneious/chimeras/single/{wildcards.sample}_calls.vcf.gz > {output} && touch {output}
        else
            touch {output}
        fi
        """

rule adjust_concensus_paired:
    priority: 2
    input:
        "Geneious/chimeras/paired/{sample}_concensus.fasta"
    output:
        "Postprocess/chimeras/paired/{sample}_concensus.fasta"
    run:
        subprocess.call(f"touch {output} ", shell=True)
        myfile = open(f'{input}')
        with open(f'{output}',"r+") as newtest:
            for line in myfile:
                if line.startswith(">"):
                    newtest.write("".join(f">{wildcards.sample}" + "_" + line.strip(">")))
                else:
                    newtest.write(line)

rule adjust_concensus_single:
    priority: 2
    input:
        "Geneious/chimeras/single/{sample}_concensus.fasta"
    output:
        "Postprocess/chimeras/single/{sample}_concensus.fasta"
    run:
        subprocess.call(f"touch {output} ", shell=True)
        myfile = open(f'{input}')
        with open(f'{output}',"r+") as newtest:
            for line in myfile:
                if line.startswith(">"):
                    newtest.write("".join(f">{wildcards.sample}" + "_" + line.strip(">")))
                else:
                    newtest.write(line)


onsuccess:
    with open("process_dag.txt","w") as f:
        f.writelines(str(workflow.persistence.dag))
    shell("cat process_dag.txt | dot -Tpng > Geneious.png | rm process_dag.txt")
