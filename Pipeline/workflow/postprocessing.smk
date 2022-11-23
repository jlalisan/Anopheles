configfile: "config.yaml"

rule all:
    input:
        expand()

rule Reference_builder_single:
    input:
        "contigs/single/{sample}.Contigs.fasta",
        "Geneious/chimeras/single/{sample}_chimera.fasta"
    output:
        "Postprocess/References/single/log/{sample}.log.file"
    shell:
        """
        f"grep -A 1 {input[0]} {input[1]} --no-group-separator > {wildcards.sample}/reference.fa
        """

rule Reference_builder_paired:
    input:
        "contigs/paired/{sample}.Contigs.fasta",
        "Geneious/chimeras/paired/{sample}_chimera.fasta"
    output:
        "Postprocess/References/paired/log/{sample}.log.file"
    shell:
        """
        f"grep -A 1 {input[0]} {input[1]} --no-group-separator > {wildcards.sample}/reference.fa
        """

rule Bowtie2_single:
    input:
        ""
    output:
        ""
    shell:
        """

        """
rule Bowtie2_paired:
    input:
        ""
    output:
        ""
    shell:
        """

        """
rule SamtoBam:
    input:
        ""
    output:
        ""
    shell:
        """

        """
rule Reference_grabber:
    input:
        ""
    output:
        ""
    shell:
        """

        """
rule Process_segments:
    input:
        ""
    output:
        ""
    shell:
        """

        """
rule Rcaller:
    input:
        ""
    output:
        ""
    shell:
        """

        """