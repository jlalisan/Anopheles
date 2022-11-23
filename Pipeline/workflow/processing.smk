rule all:
    input:
        expand()

rule Bowtie_index_single:
    input:
        "blastoutput/fetched/single/{sample}_accession.fasta",
        "contigs/single/{sample}.Contigs.fasta",
    output:
        "Geneious/Bowtie2/single/btbuild.log",
    log:
        "Geneious/Bowtie2/single/log/Bowtie.log"
    benchmark:
        "Geneious/Bowtie2/single/benchmark/Bowtiebench.csv"
    shell:
        """(
        bowtie2-build {input[1]} Bowtie2/single/{wildcards.sample}/ref_genome_btindex > {output}
        ) >{log} 2>&1"""

rule Bowtie_index_paired:
    input:
        "blastoutput/fetched/paired/{sample}_accession.fasta",
        "contigs/paired/{sample}.Contigs.fasta",
    output:
        "Geneious/Bowtie2/paired/btbuild.log",
    log:
        "Geneious/Bowtie2/paired/log/Bowtie.log"
    benchmark:
        "Geneious/Bowtie2/paired/benchmark/Bowtiebench.csv"
    shell:
        """(
        bowtie2-build {input[1]} Bowtie2/paired/{wildcards.sample}/ref_genome_btindex > {output}
        ) >{log} 2>&1"""

rule Bowtie2_single:
    input:
        "contigs/single/{sample}.Contigs.fasta"
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
        bowtie2 -x {params.btindex} {input} {output}  2>> {log} && touch {output}
        """

rule Bowtie2_paired:
    input:
        "contigs/paired/{sample}.Contigs.fasta"
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
        bowtie2 -x {params.btindex} {input} {output}  2>> {log} && touch {output}
        """

rule create_concensus_paired:
    input:
        "blastoutput/fetched/paired/{sample}_accession.fasta",
        "Geneious/Bowtie2/paired/mapped/{sample}.sam"
    output:
        "Geneious/chimeras/paired/{sample}_concensus.fasta"
    shell:
        """
        samtools view -bS {input[1]}| samtools sort - -o {wildcards.sample}_bowtie.bam
        bcftools mpileup -Ou -f {input[0]} {wildcards.sample}_bowtie.bam | bcftools call -mv -Oz -o {wildcards.sample}_calls.vcf.gz | bcftools index {wildcards.sample}_calls.vcf.gz
        bcftools norm -f {input[0]} {wildcards.sample}_calls.vcf.gz -Ob -o {wildcards.sample}_calls.norm.bcf
        bcftools filter --IndelGap 5 {wildcards.sample}_calls.norm.bcf -Ob -o {wildcards.sample}_calls.norm.flt-indels.bcf
        cat {input[0]} | bcftools consensus {wildcards.sample}_calls.vcf.gz > {output} && touch {output}
        """

rule create_concensus_single:
    input:
        "blastoutput/fetched/single/{sample}_accession.fasta",
        "Geneious/Bowtie2/single/mapped/{sample}.sam"
    output:
        "Geneious/chimeras/single/{sample}_concensus.fasta"
    shell:
        """
        samtools view -bS {input[1]}| samtools sort - -o {wildcards.sample}_bowtie.bam
        bcftools mpileup -Ou -f {input[0]} {wildcards.sample}_bowtie.bam | bcftools call -mv -Oz -o {wildcards.sample}_calls.vcf.gz | bcftools index {wildcards.sample}_calls.vcf.gz
        bcftools norm -f {input[0]} {wildcards.sample}_calls.vcf.gz -Ob -o {wildcards.sample}_calls.norm.bcf
        bcftools filter --IndelGap 5 {wildcards.sample}_calls.norm.bcf -Ob -o {wildcards.sample}_calls.norm.flt-indels.bcf
        cat {input[0]} | bcftools consensus {wildcards.sample}_calls.vcf.gz > {output} && touch {output}
        """