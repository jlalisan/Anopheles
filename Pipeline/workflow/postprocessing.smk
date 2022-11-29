configfile: "config.yaml"

rule alldo:
    input:
        expand("Postprocess/processed/{sample}.vcf",sample=config['paired'])

rule Reference_builder_paired:
    input:
        "Postprocess/chimeras/paired/{sample}_concensus.fasta"
    output:
        "Postprocess/references/paired/{sample}/{sample}_reference.fasta"
    shell:
        """
        grep -A 1 "{wildcards.sample}" {input} --no-group-separator > {output} && touch {output}
        """

rule Reference_builder_single:
    input:
        "Postprocess/chimeras/single/{sample}_concensus.fasta"
    output:
        "Postprocess/references/single/{sample}/{sample}_reference.fasta"
    shell:
        """
        grep -A 1 "{wildcards.sample}" {input} --no-group-separator > {output} && touch {output}
        """

rule Bowtie2_index_paired:
    input:
        "Postprocess/references/paired/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/references/paired/{sample}/btbuild.log"
    log:
        "Postprocess/references/paired/{sample}/{sample}.log"
    shell:
        """(
        bowtie2-build --threads 35 -f --quiet {input} Postprocess/references/paired/{wildcards.sample}/ref_genome_btindex > {output}
        ) >{log} 2>&1"""

rule Bowtie2_paired_post:
    input:
        "Bowtie2/paired/{sample}_unmapped.1.fastq",
        "Bowtie2/paired/{sample}_unmapped.2.fastq",
        "Bowtie2/paired/{sample}_unpaired.unmapped.fastq",
        "Postprocess/references/paired/{sample}/btbuild.log"
    output:
        "Postprocess/Bowtie2/paired/{sample}/{sample}.sam"
    params:
        btindex="Postprocess/references/paired/{sample}/ref_genome_btindex"
    log:
        "Postprocess/Bowtie2/paired/{sample}/map_log_file.txt"
    shell:
        """
        bowtie2 -p 35 --very-sensitive -x {params.btindex} --fr -1 {input[0]} -2 {input[1]} -U {input[2]} 1> {output} 2>> {log}
        """

rule sam_to_bam_paired:
    input:
        "Postprocess/Bowtie2/paired/{sample}/{sample}.sam"
    output:
        temp("Postprocess/Bowtie2/paired/{sample}/{sample}_sort.bam"),
        "Postprocess/Bowtie2/paired/{sample}/{sample}_sort.REF_unmapped.bam"
    log:
        "Postprocess/Bowtie2/paired/{sample}/{sample}_samtobam.log"
    shell:
        """(
		samtools view --threads 20 -bSh -o {wildcards.sample}.bam {input}
		samtools sort --threads 20 {wildcards.sample}.bam -o {output[0]}
		bamtools split -in {output[0]} -reference
        ) >{log} 2>&1"""

rule ref_paired:
    input:
        "Postprocess/references/paired/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/references/paired/{sample}/{sample}_refs.txt"
    shell:
        """
        grep ">" {input} | sed 's/.*_//' > {output}
		awk '/^>/ {{OUT=substr($0,2) ".fa"}}; OUT {{print >OUT}}' {input}
        """

rule process_paired:
    input:
        "Postprocess/references/paired/{sample}/{sample}_refs.txt",
        "Postprocess/Bowtie2/paired/{sample}/{sample}_sort.REF_unmapped.bam",
        "Postprocess/references/paired/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/processed/{sample}_long.cov",
        "Postprocess/processed/{sample}.cov",
        "Postprocess/processed/{sample}.vcf"

    shell:
        """
        samtools index {input[1]}
        bedtools genomecov -d -ibam {input[1]} > {output[0]}
        grep {wildcards.sample} {output[0]} > {output[1]}
        samtools faidx {input[2]}
        ./lofreq call -f {input[2]} -o {output[2]} {input[1]}
        """

rule Rcaller_cov:
    input:
        ""
    output:
        ""
    shell:
        """

        """

rule Rcaller_cor:
    input:
        ""
    output:
        ""
    shell:
        """

        """