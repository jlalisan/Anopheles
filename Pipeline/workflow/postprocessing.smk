configfile: "config.yaml"

rule alldo:
    input:
        expand("Postprocess/processed/cov_{sample}.png", sample=config['samples'])

rule Reference_builder:
    input:
        "Postprocess/chimeras/single/{sample}_concensus.fasta",
        "Postprocess/chimeras/paired/{sample}_concensus.fasta"
    output:
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    shell:
        """
        grep -A 1 "{wildcards.sample}" {input} --no-group-separator > {output} && touch {output}
        """

rule Bowtie2_index:
    input:
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/references/{sample}/btbuild.log"
    log:
        "Postprocess/references/{sample}/{sample}.log"
    shell:
        """(
        bowtie2-build --threads 35 -f --quiet {input} Postprocess/references/{wildcards.sample}/ref_genome_btindex > {output} && touch {output}
        ) >{log} 2>&1"""


rule Bowtie2_paired_post:
    input:
        "Bowtie2/paired/{sample}_unmapped.1.fastq",
        "Bowtie2/paired/{sample}_unmapped.2.fastq",
        "Bowtie2/paired/{sample}_unpaired.unmapped.fastq",
        "Postprocess/references/{sample}/btbuild.log"
    output:
        "Postprocess/Bowtie2/{sample}/{sample}.sam"
    params:
        btindex="Postprocess/references/{sample}/ref_genome_btindex"
    log:
        "Postprocess/Bowtie2/{sample}/map_log_file.txt"
    shell:
        """
        bowtie2 -p 35 --very-sensitive -x {params.btindex} --fr -1 {input[0]} -2 {input[1]} -U {input[2]} 1> {output} 2>> {log} && touch {output}
        """

rule Bowtie2_single_post:
    input:
        "Bowtie2/single/{sample}_single.mapped.fastq",
        "Postprocess/references/{sample}/btbuild.log"
    output:
        "Postprocess/Bowtie2/{sample}/{sample}.sam"
    params:
        btindex="Postprocess/references/{sample}/ref_genome_btindex"
    log:
        "Postprocess/Bowtie2/{sample}/map_log_file.txt"
    shell:
        """
		bowtie2 -p 35 --very-sensitive -x {params.btindex} {input[0]} 1> {output} >> {log} && touch {output}  
		"""

rule sam_to_bam:
    input:
        "Postprocess/Bowtie2/{sample}/{sample}.sam"
    output:
        "Postprocess/Bowtie2/{sample}/{sample}_sort.bam",
        "Postprocess/Bowtie2/{sample}/{sample}_sort.REF_unmapped.bam"
    log:
        "Postprocess/Bowtie2/{sample}/{sample}_samtobam.log"
    shell:
        """(
		samtools view --threads 20 -bSh -o {wildcards.sample}.bam {input}
		samtools sort --threads 20 {wildcards.sample}.bam -o {output[0]}
		bamtools split -in {output[0]} -reference && touch {output}
        ) >{log} 2>&1"""

rule ref:
    input:
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/references/{sample}/{sample}_refs.txt"
    shell:
        """
        grep ">" {input} | sed 's/.*_//' > {output}
		awk '/^>/ {{OUT=substr($0,2) ".fa"}}; OUT {{print >OUT}}' {input} && touch {output}
        """

rule process:
    input:
        "Postprocess/references/{sample}/{sample}_refs.txt",
        "Postprocess/Bowtie2/{sample}/{sample}_sort.REF_unmapped.bam",
        "Postprocess/references/{sample}/{sample}_reference.fasta"
    output:
        "Postprocess/processed/{sample}/{sample}_long.cov",
        "Postprocess/processed/{sample}/{sample}.cov",
        "Postprocess/processed/{sample}/{sample}.vcf"
    log:
        "Postprocess/processed/log/{sample}.log"
    shell:
        """
        if [ -s {input[1]} ]; then
            samtools index {input[1]}
            bedtools genomecov -d -ibam {input[1]} > {output[0]}
            grep {wildcards.sample} {output[0]} > {output[1]}
            samtools faidx {input[2]}
            ./lofreq call -f {input[2]} -o {output[2]} {input[1]} && touch {output} 2> {log}
            2>&1
        else
            touch {output}
        fi

        """


rule Rscriptcall:
    input:
        file1="Postprocess/processed/{sample}/{sample}.cov",
        file2="Postprocess/processed/{sample}/{sample}.vcf"
    output:
        outfile="Postprocess/processed/cov_{sample}.png",
        outfile2="Postprocess/processed/vcf_{sample}.png"
    script:
        "#./newtest.R" # Needs to be build
