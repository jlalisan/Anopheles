#! usr/bin/env python3

"""
Script to recreate the Geneious steps originally done manually
"""

import os
import subprocess
import re

def esearcher():
    """ Esearches and efetches the accession numbers. """
    my_srr = []

    for line in os.listdir("blastoutput//accessions"):
        if line.startswith("SRR"):
            my_srr.append(line)

    for srr in my_srr:
        with open(f"blastoutput//accessions/{srr}") as srr_file:
            for acc in srr_file:
                newacc = re.findall("[A-Z]+?[0-9 A-Z.]+[0-9.]+?", acc)
                for myacc in newacc:
                    seqfetch = f"esearch -db nucleotide -query {myacc} | efetch -format fasta > blastoutput/fetched/{srr.split('_')[0]}_{myacc}.fasta"
                    subprocess.call(seqfetch, shell=True)


def Bowtie2():
    refs = "Input_references.fasta"
    accession = "SRR6155880_exo_virus_contigs.fasta"
    subprocess.call(f"bowtie2-build {refs} Bowtie2/ref_genome_btindex > Bowtie2/btbuild.log")
    subprocess.call(f"bowtie2 -p 35 --end-to-end -x Bowtie2/ref_genome_btindex {accession} -S {accession.split('_')[0]}_single_unmapped.fastq 2>> mapped_log.file")


def main():
    Bowtie2()


if __name__ == "__main__":
    main()
