#! usr/bin/env python3

"""
Config script
"""
path_to_tools = dict(
    sra_tools="/media/studentVEE/Data/Lisan/Tools/sratoolkit.3.0.0-ubuntu64/bin",
    ray="/media/studentVEE/Data/Lisan/Tools/Ray-2.3.1",
    diamondblast="/media/studentVEE/Data/Lisan/Tools/Diamond_blast/diamond"
)

path_to_dirs = dict(
    workdir="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/Pipeline/scripts/Python/",
    basepipe="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/Pipeline/",
    sra_dir="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/Pipeline/scripts/Python/sra_files/",
    fastqdir="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/Pipeline/scripts/Python/fastq_files",
    trimjar="/media/studentVEE/Data/Lisan/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar",
    adapter="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/all_adapters.fa",
    bowindex="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/Pipeline/scripts/Python/Bowtie2/ref_genome_btindex",
    denovodir="/media/studentVEE/Data/Lisan/PycharmProjects/Pipeline/Pipeline/scripts/Python/denovo/",
    blastdir="/media/studentVEE/Data/Lisan/Tools/Diamond_blast/",
    diamnond="/media/studentVEE/Data/diamond_blast_nr"
)

path_to_files = dict(
    samplefile=path_to_dirs['basepipe'] + "resources/testtable",
    path_paired=path_to_dirs['workdir'] + "trimmed_reads/paired",
    path_single=path_to_dirs['workdir'] + "trimmed_reads/single",
    refgen=path_to_dirs['basepipe'] + "resources/refgen_t.fasta",
    mapped_p=path_to_dirs['workdir'] + "trimmed_reads/paired/mapped",
    mapped_s=path_to_dirs['workdir'] + "trimmed_reads/single/mapped",
    denovo_s=path_to_dirs['denovodir'] + "single",
    denovo_p=path_to_dirs['denovodir'] + "paired",
    contigs=path_to_dirs['workdir'] + "contigs"

)

params = dict(
    single_trim=":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15",
    paired_trim=":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45"
)