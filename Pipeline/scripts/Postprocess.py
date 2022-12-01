#! usr/bin/env python3

"""
Post processing script >> Currently skeleton script
"""

import os
import subprocess


def refmaker():
    """ Makes the reference files for used for Bowtie2 """
    workdir = os.listdir("path/to/files")  # reads from denovo and chimera needed here
    os.chdir("path/to/files")

    srr = []
    chimera = []
    for item in workdir:
        if item.startswith("SRR"):
            srr.append(item)
        else:
            chimera.append(item)

    for index in range(len(srr)):
        mysrr = srr[index]
        # Grabs the references.
        makeref = f"grep -A 1 {mysrr} chimeras/chimeras.fasta --no-group-separator > {mysrr}/reference.fa"
        subprocess.call(makeref, shell=True)


def Bowtie2():
    """ Maps the single and paired reads against the index. """
    workdir = os.listdir("path/to/files")
    os.chdir("path/to/files")

    srrlist = []

    for item in workdir:
        if item.startswith("SRR"):
            srrlist.append(item)

    for index in range(len(srrlist)):
        srr = srrlist[index]

        # Paired first
        pairedbowtie2 = f"bowtie2 -p 35 --very-sensitive -x $PROC_DIR/REF/${srr}/reference-index " \
                        f"--fr -1 {srr}_1.unmapped.fastq -2 {srr}_2.unmapped.fastq -U " \
                        f"{srr}_unpaired.unmapped.fastq 1> $PROC_DIR/{srr}.sam 2>> $PROC_DIR/map_log_file.txt"
        # Single second
        singlebowtie = f"bowtie2 -p 35 --very-sensitive -x $PROC_DIR/REF/${srr}/reference-index ${srr}_single.unmapped.fastq " \
                       f"1> $PROC_DIR/${srr}.sam 2>> $PROC_DIR/map_log_file.txt"

        print(f"Mapping {srr}")
        subprocess.call(pairedbowtie2, shell=True)
        print(f"Mapping {srr}")
        subprocess.call(singlebowtie, shell=True)


def samtobam():
    """ Converts sam files to bam files """
    workdir = os.listdir("path/to/files") # Sam files
    os.chdir("path/to/files")

    samfiles = []

    for items in workdir:
        if items.endswith(".sam"):
            samfiles.append(items)

    for index in range(len(samfiles)):
        samfile = samfiles[index]
        # SEE IF THIS CAN BE PIPED INTO ONE COMMAND
        samunique = f"cp {samfile.split('.')[0]}.sam ./{samfile.split('.')[0]}_unique.sam"
        subprocess.call(samunique, shell=True)

        samview = f"samtools view --threads 20 -bSh -o {samfile.split('.')[0]}.bam {samfile.split('.')[0]}_unique.sam"
        subprocess.call(samview, shell=True)

        samsort = f"samtools sort --threads 20 ${samfile.split('.')[0]}.bam -o ${samfile.split('.')[0]}_sort.bam"
        subprocess.call(samsort, shell=True)

        bamsplit = f"bamtools split -in ${samfile.split('.')[0]}_sort.bam -reference | rm {samfile.split('.')[0]}_unique.sam {samfile.split('.')[0]}_sort.bam"
        subprocess.call(bamsplit, shell=True)


def referencegrabber():
    """ Grabs the correct references and splits these."""
    workdir = os.listdir("path/to/files") # Reference files
    os.chdir("path/to/files")

    srr_names = []

    for items in workdir:
        if items.startswith("SRR"):
            srr_names.append(items)

    for index in range(len(srr_names)):
        srr = srr_names[index]

        grab = f"grep > reference.fa | sed 's/.*_//' > ../{srr}_refs.txt"
        split = "'/^>/ {OUT=substr($0,2) .fa}; OUT {print >OUT}' reference.fa"

        subprocess.call(grab, shell=True)
        subprocess.call(split, shell=True)

def process_segments():
    """ Processes the chimeras and the srr files. """
    workdir = os.listdir("path/to/files") # chimera / sorted ref files
    os.chdir("path/to/files")

    chimeras = []
    srrfiles = []

    for items in workdir:
        if items.startswith("SRR"):
            srrfiles.append(items)
        if "chimera" in items:
            chimeras.append(items)

    chimeras.sort()
    srrfiles.sort()

    for index in range(len(srrfiles)):
        srr = srrfiles[index]
        chim = chimeras[index]

        filemover = f"mv {srr}_sort.REF_{srr}_{chim}.bam {srr}_{chim}_sort.bam"
        subprocess.call(filemover, shell=True)

        indexer = f"samtools index {srr}_{chim}_sort.bam | bedtools genomecov -d -ibam {srr}_{chim}_sort.bam > " \
                  f"{srr}_{chim}_long.cov | grep {srr}_{chim} {srr}_{chim}_long.cov > {srr}_{chim}.cov | " \
                  f"rm {srr}_{chim}_long.cov | samtools faidx {srr}_{chim}.fa"
        subprocess.call(indexer, shell=True)

        lofreq = f"lofreq call -f {srr}_{chim}.fa -o {srr}_{chim}.vcf {srr}_{chim}_sort.bam"
        subprocess.call(lofreq, shell=True)

        mover = f"mv {srr}_{chim}.vcf $PROC_DIR / PROCESSED / DEPTH_VARIANTS /{srr}_{chim}.vcf | cp $PROC_DIR / DIR_${srr} / *.cov, *.vcf $PROC_DIR / PROCESSED / CORRECTION | cp $PROC_DIR / REF /${srr} / *.fa $PROC_DIR / PROCESSED / CORRECTION"
        subprocess.call(mover, shell=True)


def Rcaller():
    """ Calls upon R scripts to make visualisations. """
    coverage = "Rscript coverage.R 2>> R_log_file.txt"
    subprocess.call(coverage, shell=True)
    listsample = "ls *.cov | sed 's/.cov//' >> $PROC_DIR/PROCESSED/CORRECTION/listsample.txt"
    subprocess.call(listsample, shell=True)
    variants = "Rscript $PROC_DIR/correction.R 2>> $PROC_DIR/PROCESSED/CORRECTION/R_log_file.txt"
    subprocess.call(variants, shell=True)


def main():
    """ Main function to call all the other functions. """
    refmaker()
    Bowtie2()
    samtobam()
    referencegrabber()
    process_segments()
    Rcaller()


if __name__ == "__main__":
    main()
