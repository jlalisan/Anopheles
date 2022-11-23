import os
import subprocess
import re

configfile: "config.yaml"

rule all:
    input:
        expand("blastoutput/fetched/single/{sample}_accession.fasta", sample=config['single']),
        expand("blastoutput/fetched/paired/{sample}_accession.fasta", sample=config['paired']),


rule prefetching:
    priority: 10
    output:
        "sra_files/{sample}/{sample}.sra"
    log:
        "sra_files/logs/{sample}.log"
    shell:
        """
        (prefetch {wildcards.sample} -O sra_files) > {log} 2>&1 && touch {output}
        """

rule fastqdump:
    priority: 9
    input:
        "sra_files/{sample}/{sample}.sra"
    output:
        "fastq_files/log/{sample}.log",
        "fastq_files/{sample}.fastq",
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq"
    shell:
        """
        fasterq-dump -S -O fastq_files/ -t fastq_files/ --split-3 {input} > {output} 2>&1 && touch {output}
        """

rule bowtieindex:
    priority: 10
    input:
        "refgen_t.fasta"
    output:
        "Bowtie2/btbuild.log",
    log:
        "Bowtie2/log/Bowtie.log"
    benchmark:
        "Bowtie2/benchmarks/Bowtiebench.csv"
    shell:
        """
        (bowtie2-build {input} Bowtie2/ref_genome_btindex > {output}) >{log} 2>&1
        """

rule trimmomatic_paired:
    priority: 8
    input:
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq",
        "Bowtie2/btbuild.log",
    output:
        "trimmomatic/paired/{sample}_1.trim.fastq",
        "trimmomatic/paired/{sample}_1.untrim.fastq",
        "trimmomatic/paired/{sample}_2.trim.fastq",
        "trimmomatic/paired/{sample}_2.untrim.fastq"
    message:
        "Started read trimming for {wildcards.sample}!"
    log:
        "trimmomatic/paired/logs/{sample}_trimmed.log"
    shell:
        """
        (java -jar trimmomatic-0.39.jar PE -phred33 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45) >{log} 2>&1 && touch {output}
        echo Btindex build see {input[2]}
        """
rule trimmomatic_single:
    priority: 8
    input:
        "fastq_files/{sample}.fastq",
        "Bowtie2/btbuild.log"
    output:
        "trimmomatic/single/{sample}.trim.fastq"
    message:
        "Started read trimming for {wildcards.sample}!"
    log:
        "trimmomatic/single/logs/{sample}_trimmed.log"
    shell:
        """
        (java -jar trimmomatic-0.39.jar SE -phred33 {input[0]} {output} ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15) >{log} 2>&1 && touch {output}
        echo {input[1]}
        """
rule bowtiemerger:
    priority: 7
    input:
        "trimmomatic/paired/{sample}_1.untrim.fastq",
        "trimmomatic/paired/{sample}_2.untrim.fastq"
    output:
        "trimmomatic/merged/{sample}_merged.fastq"
    shell:
        """
        cat {input[0]} {input[1]} > {output}
        """
rule bowtiesingle:
    priority: 6
    input:
        "trimmomatic/single/{sample}.trim.fastq"
    output:
        "Bowtie2/single/{sample}_single.mapped.fastq",
        "Bowtie2/single/{sample}.sam"
    params:
        btindex="Bowtie2/ref_genome_btindex"
    log:
        "Bowtie2/log/{sample}_unmapped.log"
    shell:
        """
        bowtie2 -p 10 --end-to-end -x {params.btindex} {input} --al {output[0]} 1> {output[1]} 2>> {log} && touch {output}
        """

rule bowtiepaired:
    priority: 6
    input:
        "trimmomatic/paired/{sample}_1.trim.fastq",
        "trimmomatic/paired/{sample}_2.trim.fastq",
        "trimmomatic/merged/{sample}_merged.fastq",
    output:
        "Bowtie2/paired/{sample}_unmapped.fastq",
        "Bowtie2/paired/{sample}_unpaired.unmapped.fastq",
        "Bowtie2/paired/{sample}.sam",
        "Bowtie2/paired/{sample}_unmapped.1.fastq",
        "Bowtie2/paired/{sample}_unmapped.2.fastq",
    params:
        btindex="Bowtie2/ref_genome_btindex"
    log:
        "Bowtie2/log/{sample}_unmapped.log"
    shell:
        """
        bowtie2 -p 20 --end-to-end -x {params.btindex} --fr -1 {input[0]} -2 {input[1]} -U {input[2]} --al-conc {output[0]} --al {output[1]} 1> {output[2]} 2>> {log} 2>&1 && touch {output}
        """

rule denovo_single:
    priority: 5
    input:
        "Bowtie2/single/{sample}_single.mapped.fastq"
    output:
        directory("denovo/single/{sample}.forblast"),
        "contigs/single/{sample}.Contigs.fasta"
    log:
        "denovo/single/log/{sample}_log.file"
    shell:
        """
        mpirun -n 2 Ray -s {input} -o {output} 1>/dev/null 2>> {log} && touch {output[0]}
        cat {output[0]}/Contigs.fasta > {output[1]} && touch {output[1]}
        """

# Get contigs for blast
rule denovo_paired:
    priority: 5
    input:
        "Bowtie2/paired/{sample}_unpaired.unmapped.fastq",
        "Bowtie2/paired/{sample}_unmapped.1.fastq",
        "Bowtie2/paired/{sample}_unmapped.2.fastq"
    output:
        directory("denovo/paired/{sample}.forblast"),
        directory("denovo/paired/{sample}_u.forblast"),
        "contigs/paired/{sample}.Contigs.fasta"
    log:
        "denovo/paired/log/{sample}_u.log.file",
        "denovo/paired/log/{sample}.log.file"
    shell:
        """
        mpiexec  -n 2 Ray -s {input[0]} -o {output[0]} 1>/dev/null 2>> {log[0]} && touch {output[0]}
        mpiexec  -n 2 Ray -p {input[1]} {input[2]} -o {output[1]} 1>/dev/null 2>> {log[1]} && touch {output[1]}
        cat {output[0]}/Contigs.fasta {output[1]}/Contigs.fasta > {output[2]} && touch {output[2]}
        """

rule blasting_single:
    priority: 4
    input:
        "contigs/single/{sample}.Contigs.fasta"
    output:
        "blastoutput/single/{sample}_matches.tsv"
    params:
        diamond="spikeddb.dmnd"
    run:
        if os.stat(f"{input}").st_size > 1:
            subprocess.call(f"./diamond blastx -d  {params.diamond} -q {input} --sensitive --quiet -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids skingdoms sscinames stitle -k 1 -o {output} && touch {output}", shell=True)
        else:
            subprocess.call(f"touch {output}", shell=True)

rule blasting_paired:
    priority: 4
    input:
        "contigs/paired/{sample}.Contigs.fasta"
    output:
        "blastoutput/paired/{sample}_matches.tsv"
    params:
        diamond="spikeddb.dmnd"
    run:
        if os.stat(f"{input}").st_size > 1:
            subprocess.call(f"./diamond blastx -d  {params.diamond} -q {input} --sensitive --quiet -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids skingdoms sscinames stitle -k 1 -o {output} && touch {output}", shell=True)
        else:
            subprocess.call(f"touch {output}", shell=True)

rule keywords:
    priority: 10
    output:
        "Keywords.txt"
    script:
        "keyword.py"

rule blastmatcher_single:
    priority: 3
    input:
        "blastoutput/single/{sample}_matches.tsv",
        "Keywords.txt"
    output:
        "blastoutput/single/{sample}_ortho_matches.tsv",
        "blastoutput/accessions/single/{sample}_acc_hits.txt"
    run:
        if os.stat(f"{input[0]}").st_size > 1:
            subprocess.call(f"grep -f {input[1]} {input[0]} > {output[0]}", shell=True)
            subprocess.call("awk '{{print $2}}' " + f"{output[0]} | sort | uniq > {output[1]} && touch {output}", shell=True)
        else:
            subprocess.call(f"touch {output}", shell=True)

rule blastmatcher_paired:
    priority: 3
    input:
        "blastoutput/paired/{sample}_matches.tsv",
        "Keywords.txt"
    output:
        "blastoutput/paired/{sample}_ortho_matches.tsv",
        "blastoutput/accessions/paired/{sample}_acc_hits.txt"
    run:
        if os.stat(f"{input[0]}").st_size > 1:
            subprocess.call(f"grep -f {input[1]} {input[0]} > {output[0]}", shell=True)
            subprocess.call("awk '{{print $2}}' " + f"{output[0]} | sort | uniq > {output[1]} && touch {output}", shell=True)
        else:
            subprocess.call(f"touch {output}", shell=True)


rule efetcher_paired:
    priority: 2
    input:
        "blastoutput/accessions/paired/{sample}_acc_hits.txt"
    benchmark:
        "blastoutput/fetched/paired/bench/{sample}.bench.txt"
    run:
        my_srr = []
        my_acc = []

        for files in os.listdir("blastoutput//accessions//paired"):
            if files.startswith("SRR"):
                my_srr.append(files)

        for srr in my_srr:
            with open(f"blastoutput//accessions//paired//{srr}") as srr_file:
                for acc in srr_file:
                    newacc = re.findall("[A-Z]+?[0-9 A-Z.]+[0-9.]+?",acc)
                    for myacc in newacc:
                        my_acc.append(myacc)
                        seqfetch = f"esearch -db nucleotide -query {myacc} | efetch -format fasta > blastoutput//fetched//paired//{srr.split('_')[0]}_{myacc}.fasta"
                        subprocess.call(seqfetch,shell=True)


rule efetcher_single:
    priority: 2
    input:
        "blastoutput/accessions/single/{sample}_acc_hits.txt"
    benchmark:
        "blastoutput/fetched/single/bench/{sample}.bench.txt"
    run:
        my_srr = []
        my_acc = []

        for files in os.listdir("blastoutput//accessions//single"):
            if files.startswith("SRR"):
                my_srr.append(files)

        for srr in my_srr:
            with open(f"blastoutput//accessions//single//{srr}") as srr_file:
                for acc in srr_file:
                    newacc = re.findall("[A-Z]+?[0-9 A-Z.]+[0-9.]+?",acc)
                    for myacc in newacc:
                        my_acc.append(myacc)
                        seqfetch = f"esearch -db nucleotde -query {myacc} | efetch -format fasta > blastoutput//fetched//single//{srr.split('_')[0]}_{myacc}.fasta"
                        subprocess.call(seqfetch,shell=True)

rule merge_acc_single:
    input:
        "blastoutput/fetched/single/bench/{sample}.bench.txt"
    output:
        "blastoutput/fetched/single/{sample}_accession.fasta"
    run:
        if os.path.exists(f"blastoutput/fetched/single/{wildcards.sample}*"):
            print(os.stat("blastoutput/fetched/single/").st_size)
            subprocess.call(f"cat blastoutput/fetched/single/{wildcards.sample}* >> {output} && touch {output}")
        else:
            subprocess.call(f"touch {output}",shell=True)

rule merge_acc_paired:
    input:
        "blastoutput/fetched/paired/bench/{sample}.bench.txt"
    output:
        "blastoutput/fetched/paired/{sample}_accession.fasta"
    run:
        if os.path.exists(f"blastoutput/fetched/paired/{wildcards.sample}*"):
            subprocess.call(f"cat blastoutput/fetched/paired/{wildcards.sample}* >> {output} && touch {output}")
        else:
            subprocess.call(f"touch {output}",shell=True)

onsuccess:
    with open("dag.txt","w") as f:
        f.writelines(str(workflow.persistence.dag))
    shell("cat dag.txt | dot -Tpng > dag.png | rm dag.txt")
