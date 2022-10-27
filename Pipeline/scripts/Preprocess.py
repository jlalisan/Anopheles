#! usr/bin/env python3

"""
Script to preprocess Virus data from NCBI accession numbers.

"""

import config
import os
import subprocess
import glob


def get_accession():
    """ Function to retrieve all the accession numbers and return these. """
    accession = {f.strip() for f in open(config.path_to_files["samplefile"])}
    return accession


def create_output_dirs():
    """ Creates all the output directories. """

    print("creating directories")
    # Config is a seperate file for the paths that can be adjusted
    work_dir = config.path_to_dirs["workdir"]

    # Gives them a place in memory using f-strings
    sra_files = f"{work_dir}sra_files/logs"
    fastq_files = f"{work_dir}fastq_files/logs"
    trimmed_s = f"{work_dir}trimmed_reads/single/mapped"
    trimmed_p = f"{work_dir}trimmed_reads/paired/mapped"
    trimmed_l = f"{work_dir}trimmed_reads/logs"
    bowtie_i = f"{work_dir}Bowtie2/results"
    denovo_s = f"{work_dir}/denovo/single"
    denovo_p = f"{work_dir}/denovo/paired"
    contigs = f"{work_dir}/contigs"

    # Makes all the directories
    mylist = (sra_files, fastq_files, trimmed_s, trimmed_p, trimmed_l, bowtie_i,
              denovo_s, denovo_p, contigs)
    for i in mylist:
        os.makedirs(i, exist_ok=True)

    print("Directories successfully build")


def bowtie2_btbuild():
    """ Builds the references for the genome. """
    bowtieindex = "bowtie2-build " + config.path_to_files[
        "refgen"] + " Bowtie2/ref_genome_btindex > Bowtie2/results/btbuild.log"
    subprocess.call(bowtieindex, shell=True)


def prefetch(sample):
    """ Prefetches the samples given from the get_accession function. """
    for sra_id in sample:
        print("Currently downloading: " + sra_id)
        # Makes directory is not there.
        prefetcher = "prefetch " + sra_id + " --max-size 250GB -O sra_files 2>> sra_files/logs/" + sra_id + "_log.file"
        # Command line function to call the prefetch
        subprocess.call(prefetcher, shell=True)


def fastqdump(sample):
    """ Downloads the files with the SRA accessions from prefetch and get_accession. """
    for sra_id in sample:
        print("Generating fastq for: " + sra_id)
        fastq_dump = "fastq-dump --outdir fastq_files --split-3 " + config.path_to_dirs[
            "sra_dir"] + sra_id + "/" + sra_id + ".sra 2>> fastq_files/logs/" + sra_id + "_log.file"
        subprocess.call(fastq_dump, shell=True)


def trimmomatic():
    """ Single and paired trimmomatic"""
    # Single
    rawdir = os.listdir(config.path_to_dirs["fastqdir"])
    rl1 = []
    rl2 = []
    for item in rawdir:
        # Filters out the single files by checking for _
        if "_" not in item and item.startswith("SRR"):
            trimmomatic_single = "java -jar " + config.path_to_dirs[
                "trimjar"] + " SE " + "fastq_files/" + item + " trimmed_reads/single/" + item.split(".")[
                                     0] + "_single.fastq ILLUMINACLIP:" + \
                                 config.path_to_dirs["adapter"] + config.params[
                                     "single_trim"] + " 2>> trimmed_reads/logs/" + item.split(".")[0] + "_log.file"
            subprocess.call(trimmomatic_single, shell=True)
            print("Trimming:", item.split(".")[0])
        # If there is an _ it determines if its part 1 or 2
        else:
            if "_1" in item:
                rl1.append(item)
            if "_2" in item:
                rl2.append(item)
    # After adding it to the list it sorts so all the files are aligned
    rl1.sort()
    rl2.sort()
    # Paired
    for index in range(0, len(rl1)):
        # Picks the correct index for alignment
        r1 = rl1[index]
        r2 = rl2[index]

        print("Trimming:", r1.split(".")[0], "and", r2.split(".")[0])

        trimmomatic_paired = "java -jar " + config.path_to_dirs["trimjar"] + " PE fastq_files/"\
                             + r1 + " fastq_files/" + r2 + " trimmed_reads/paired/" + r1.split(".")[0]\
                             + ".trim.fastq trimmed_reads/paired/" + r1.split(".")[0]\
                             + "un.trim.fastq trimmed_reads/paired/" + r2.split(".")[0]\
                             + ".trim.fastq trimmed_reads/paired/" + r2.split(".")[0]\
                             + "un.trim.fastq ILLUMINACLIP:" + config.path_to_dirs["adapter"] \
                             + config.params["paired_trim"] + " 2>> trimmed_reads/logs/paired_log.file"
        subprocess.call(trimmomatic_paired, shell=True)


def unzipping():
    """ Unzips GZ files in case those are present. """
    os.chdir(config.path_to_files["path_single"])
    single = os.listdir(config.path_to_files["path_single"])
    if ".gz" in single:
        subprocess.call("gunzip *fastq.gz", shell=True)
    os.chdir(config.path_to_files["path_paired"])
    paired = os.listdir(config.path_to_files["path_paired"])
    if ".gz" in paired:
        subprocess.call("gunzip *fastq.gz", shell=True)
    else:
        print("No files to unzip, moving on.")


def bowtiemerge():
    """ Merges unpaired files. """
    os.chdir(config.path_to_files["path_paired"])
    mergedir = os.listdir(config.path_to_files["path_paired"])
    file1 = []
    file2 = []
    # Makes sure the paired dir is not empty
    if len(mergedir) != 0:
        for file in mergedir:
            if "un" in file and "_1" in file:
                file1.append(file)
            if "un" in file and "_2" in file:
                file2.append(file)
    file1.sort()
    file2.sort()

    for index in range(0, len(file1)):
        r1 = file1[index]
        r2 = file2[index]
        print("Merging:", r1.split("un")[0], "with", r2.split("un")[0])
        subprocess.call("cat " + r1 + " " + r2 + " 1> " + r1.split("_")[0] + "_merged_trim.fastq", shell=True)


def bowtie2():
    """ Maps reads with Bowtie2. """
    rawsingledir = os.listdir(config.path_to_files["path_single"])
    rawpaireddir = os.listdir(config.path_to_files["path_paired"])
    rl1 = []
    rl2 = []
    mergedl = []
    # Makes sure there are paired files otherwise skips.
    if len(rawpaireddir) != 0:
        for item in rawpaireddir:
            if "_1" in item and "un" not in item:
                rl1.append(item)
            if "_2" in item and "un" not in item:
                rl2.append(item)
            if "merged" in item:
                mergedl.append(item)

    rl1.sort()
    rl2.sort()
    mergedl.sort()

    for index in range(0, len(rl1)):
        r1 = rl1[index]
        r2 = rl2[index]
        merged = mergedl[index]

        os.chdir(config.path_to_files["path_paired"])

        print("Mapping:", r1.split(".")[0], "with", r2.split(".")[0], "and", merged.split("_trim")[0])

        bowtiepaired = "bowtie2 -p 35 --end-to-end -x " + config.path_to_dirs["bowindex"] + " --fr -1 " + r1 + " -2 " +\
                       r2 + " -U " + merged + " --al-conc " + config.path_to_files["mapped_p"] + "/" + \
                       r1.split("_")[0] + "_%.unmapped.fastq --al " + config.path_to_files["mapped_p"] + "/" + \
                       r1.split("_")[0] + "_unpaired.unmapped.fastq 1> " + config.path_to_files[
                           "mapped_p"] + "/" + r1.split("_")[0] + ".sam 2>> " \
                       + config.path_to_files["mapped_p"] + "/unmapped_log.file"

        subprocess.call(bowtiepaired, shell=True)

    # Maps single files
    for item in rawsingledir:
        os.chdir(config.path_to_files["path_single"])
        if "SRR" in item:
            print("Mapping:", item.split("_")[0])
            bowtie2single = "bowtie2 -p 35 --end-to-end -x " + config.path_to_dirs[
                "bowindex"] + " " + item + " -S " + config.path_to_files["mapped_s"] + "/" + item.split("_")[
                                0] + "_single_unmapped.fastq 2>> " + config.path_to_files[
                                "mapped_s"] + "/mapped_log.file"
            subprocess.call(bowtie2single, shell=True)


def denovo():
    """ Creates the contigs for the blast. """
    # single
    os.chdir(config.path_to_files["mapped_s"])
    rawsingledir = os.listdir(config.path_to_files["mapped_s"])

    correct = []

    for file in rawsingledir:
        if "_single" in file:
            correct.append(file)

    for index in range(0, len(correct)):
        myfile = correct[index]
        print("Assembling:", myfile.split("_")[0])
        novosingle = "mpiexec -n 20 Ray -s " + myfile + " -o " + config.path_to_files["denovo_s"] + "/" + \
                     myfile.split("_")[0] + ".forblast 1>/dev/null 2>> " + config.path_to_files[
                         "denovo_s"] + "/denovo_log_file"
        subprocess.call(novosingle, shell=True)

    # paired
    file1 = []
    file2 = []
    others = []

    os.chdir(config.path_to_files["mapped_p"])
    rawpaireddir = os.listdir(config.path_to_files["mapped_p"])
    for file in rawpaireddir:
        if "_1.unmapped" in file:
            file1.append(file)
        if "_2.unmapped" in file:
            file2.append(file)
        if "unpaired" in file:
            others.append(file)

    file1.sort()
    file2.sort()
    others.sort()

    for index in range(0, len(file1)):
        r1 = file1[index]
        r2 = file2[index]
        other = others[index]

        print("Assembling:", r1.split("."[0]), "with", r2.split("."[0]), "and", other.split("."[0]))

        myray = "mpiexec -n 20 Ray -p " + r1 + " " + r2 + " -o " + config.path_to_files["denovo_p"] + "/" + \
                r1.split("_")[0] + ".forblast 1>/dev/null 2>> " + config.path_to_files["denovo_p"] + "/denovo_log_file"
        myrayu = "mpiexec -n 20 Ray -s " + other + " -o " + config.path_to_files["denovo_p"] + "/" + other.split("_")[
            0] + ".forblast_u 1>/dev/null 2>> " + config.path_to_files["denovo_p"] + "/denovo_log_file"

        subprocess.call(myray, shell=True)
        subprocess.call(myrayu, shell=True)


def contigmerger():
    os.chdir(config.path_to_files["denovo_p"])
    mergedir = os.listdir(config.path_to_files["denovo_p"])

    SRR_u = []
    SRR = []

    for items in mergedir:
        if items.endswith("u"):
            SRR_u.append(items)
        if items.endswith("forblast"):
            SRR.append(items)

    SRR.sort()
    SRR_u.sort()

    for index in range(0, len(SRR_u)):
        myfile_u = SRR_u[index]
        myfile_p = SRR[index]

        statement2 = "cat " + myfile_u + "/Contigs.fasta " + myfile_p + "/Contigs.fasta > " + config.path_to_dirs[
            'workdir'] + "contigs/" + myfile_u.split(".")[0] + ".Contigs.fasta"
        subprocess.call(statement2, shell=True)


def blaster():
    myblast = ""
    subprocess.call(myblast)
    pass


def main():
    """ Main function to call all the other functions. """
    create_output_dirs()
    bowtie2_btbuild()
    # The get_accession is a required argument here.
    prefetch(get_accession())
    fastqdump(get_accession())
    trimmomatic()
    unzipping()
    bowtiemerge()
    bowtie2()
    denovo()
    contigmerger()


if __name__ == "__main__":
    main()
