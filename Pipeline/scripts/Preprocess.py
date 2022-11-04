#! usr/bin/env python3

"""
Script to preprocess Virus data from NCBI accession numbers.

"""
import glob
import shutil
import os
import subprocess
import re
import requests

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import config


def get_accession():
    """ Function to retrieve all the accession numbers and return these. """
    # Uses every accession number logged in the file
    accession = {f.strip() for f in open(config.path_to_files["samplefile"])}
    return accession


def downloader():
    """ Downloads the HTML page for the virus family of choice"""
    # Gets the correct URL and Virus family from the config
    url = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={config.path_to_words['virusfamily']}"
    data = requests.get(url)

    # Writes the entire HTML page into a txt format
    with open(f"{config.path_to_dirs['basepipe']}resources/{config.path_to_words['virusfamily']}.txt", "w") as out_f:
        out_f.write(data.text.strip())


def getter():
    """ Reworks the HTML txt file into all the viruses belonging to the chosen family"""
    viruses = []
    with open(f"{config.path_to_dirs['basepipe']}resources/{config.path_to_words['virusfamily']}.txt") as htmlfile:
        for line in htmlfile:
            # All viruses in NCBI have links so here we filter all the links.
            if line.startswith("<LI"):
                # Almost all links have the same start and before line 115 ~ 120 no names are mentioned
                first = line.strip()[115:].split(">",)[1]
                # Uses a regular expression to remove anything behind the virus name such as (Arizona/2019)
                second = re.sub("[\(\[].*?[\)\]]+", "", first.split("<")[0])
                # Append all the viruses to a list for later writing.
                viruses.append(" ".join(second.split()[:8]))
                # Clears the list from all duplicates.
                viruses = list(dict.fromkeys(viruses))

    # Opens a new file for the writing of all the virusnames
    with open(f"{config.path_to_dirs['basepipe']}resources/Keywords.txt", "w") as out_file:
        for index, line in enumerate(viruses):
            # Makes sure to skip empty lines that .strip does not find.
            if len(line) >= 2:
                # Adds an enter after every virus so they dont clutter together.
                out_file.write(line.strip() + "\n")
        # Closes the file as per standard regulations.
        out_file.close()
        htmlfile.close()
    # Removes the HTML txt file since it is now obsolete
    subprocess.call(f"rm {config.path_to_dirs['basepipe']}resources/{config.path_to_words['virusfamily']}.txt", shell=True)


def databasedownload():
    """ Downloads the databases required for this project and converts the nucleotides to aminoacids if required. """
    # All Aminoacids except for any that could be misread as an nucleotide to avoid the use of a big dict.
    protein = "RDQEHILKMFPSWYV"
    os.chdir(config.path_to_dirs['databases'])
    databasedir = os.listdir(config.path_to_dirs['databases'])
    # If the databasedir is empty then download the files, otherwise assume they present from previous run.
    if len(databasedir) <= 1:
        database = f"wget -O {config.path_to_dirs['databases']}uniref90.fasta.gz ftp://ftp.uniprot.org/pub/databases" \
                   f"/uniprot/uniref/uniref90/uniref90.fasta.gz & " \
                   f"wget -O {config.path_to_dirs['databases']}U-RVDBv24.1.fasta.gz " \
                   f"https://rvdb.dbi.udel.edu/download/U-RVDBv24.1.fasta.gz"

        subprocess.call(database, shell=True)
    # Unzip any databases for zipped onces are of no use.
    if ".gz" in databasedir:
        unzipping = f"gunzip *.fasta.gz"
        subprocess.call(unzipping, shell=True)

    for files in databasedir:
        file = open(files)
        # Read the first lines to determine if the file is nucleotide or protein
        head = [next(file) for x in range(2)]
        # Always assume nucleotide first.
        amino = False
        # Skip first line in file which is the header
        for line in head[1:]:
            # Checks if a protein is found within the line
            for aa in protein:
                # If a protein is found the amino is put to True and the assumption of a protein file is made.
                if aa in line:
                    amino = True
        # To be done if the file is not a protein file.
        if not amino:
            # To alert the user that there is a nucleotide file that is going to be converted.
            print(f"converting {files} to aminoacids")
            # Makes new file for the creation of the new protein sequence file
            with open(f"{config.path_to_dirs['databases']}virusprotdatabase.fasta", 'w') as aa_fa:
                # Parses through the files that are nucleotides
                for dna_record in SeqIO.parse(f"{config.path_to_dirs['databases']}{files}",
                                              'fasta'):
                    # Use both forward and reverce sequences
                    dna_seqs = [dna_record.seq, dna_record.seq.reverse_complement()]
                    # Generate all translation frames
                    aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in dna_seqs)
                    # Select the longest one
                    max_aa = max(aa_seqs, key=len)
                    # Write new record
                    aa_record = SeqRecord(max_aa, id=dna_record.id, description="translated sequence")
                    SeqIO.write(aa_record, aa_fa, 'fasta')
                    # Deletes the old file since this is obsolete now
                    deleteold = f"rm {config.path_to_dirs['databases']}{files}"
                subprocess.call(deleteold, shell=True)


def merger():
    """ Merges both databases creating a virus spiked one"""
    os.chdir(config.path_to_dirs['databases'])
    spikeddir = os.listdir(f"{config.path_to_dirs['blastdir']}spiked_database/")
    # If the spiked database has not been build then do this
    if "spikeddatabase.fasta" not in spikeddir:
        with open(f"{config.path_to_dirs['blastdir']}spiked_database/spikeddatabase.fasta", 'wb') as outfile:
            # Selects all the files ending with .fasta to merge
            for filename in glob.glob('*.fasta'):
                if filename == outfile:
                    # Do not want to copy the output into the wrong file
                    continue
                with open(filename, 'rb') as readfile:
                    # Merges the file
                    shutil.copyfileobj(readfile, outfile)

        print("Database has been build")
        # Creates the diamond database
        os.chdir(f"{config.path_to_dirs['blastdir']}spiked_database/")
        dmndbase = "./diamond makedb --in spikeddatabase.fasta -d spikeddb --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp"
        subprocess.call(dmndbase, shell=True)


def create_output_dirs():
    """ Creates all the output directories. """

    print("Creating directories.")
    # Config is a separate file for the paths that can be adjusted
    work_dir = config.path_to_dirs["workdir"]

    # Gives them a place in memory using f-strings
    sra_files = f"{work_dir}sra_files/logs"
    fastq_files = f"{work_dir}fastq_files/logs"
    trimmed_s = f"{work_dir}trimmed_reads/single/mapped"
    trimmed_p = f"{work_dir}trimmed_reads/paired/mapped"
    trimmed_l = f"{work_dir}trimmed_reads/logs"
    bowtie_i = f"{work_dir}Bowtie2"
    denovo_s = f"{work_dir}/denovo/single"
    denovo_p = f"{work_dir}/denovo/paired"
    contigs = f"{work_dir}/contigs"
    blastoutput = f"{work_dir}/blastoutput"
    orthomatch = f"{blastoutput}/orthomatch"
    accessions = f"{blastoutput}/accessions/fetched"

    # Makes all the directories
    mylist = (sra_files, fastq_files, trimmed_s, trimmed_p, trimmed_l, bowtie_i,
              denovo_s, denovo_p, contigs, blastoutput, orthomatch, accessions)
    for i in mylist:
        # Exist_ok is required incase the user reruns the program and just deletes or moves data
        os.makedirs(i, exist_ok=True)

    print("Directories successfully build")


def bowtie2_btbuild():
    """ Builds the references for the genome. """
    bowtieindex = f"bowtie2-build {config.path_to_files['refgen']} Bowtie2/ref_genome_btindex > Bowtie2/btbuild.log"
    subprocess.call(bowtieindex, shell=True)


def prefetch(sample):
    """ Prefetches the samples given from the get_accession function. """
    os.chdir(config.path_to_dirs['workdir'])
    for sra_id in sample:
        print(f"Currently downloading {sra_id}")
        # Makes directory if not there.
        prefetcher = f"prefetch {sra_id} --max-size 250GB -O sra_files 2>> sra_files/logs/{sra_id}_log.file"
        # Command line function to call the prefetch
        subprocess.call(prefetcher, shell=True)


def fastqdump(sample):
    """ Downloads the files with the SRA accessions from prefetch and get_accession. """
    for sra_id in sample:
        print(f"Generating fastq for: {sra_id}")
        # The --gzip can be added here and then all the files will be zipped incase required.
        fastq_dump = f"fastq-dump --outdir fastq_files --split-3 {config.path_to_dirs['sra_dir']}/{sra_id}/{sra_id}.sra 2>> fastq_files/logs/{sra_id}_log.file"
        subprocess.call(fastq_dump, shell=True)


def trimmomatic():
    """ Single and paired trimmomatic"""
    # Opens the directory as a list.
    rawdir = os.listdir(config.path_to_dirs["fastqdir"])
    rl1 = []
    rl2 = []
    # Single files trimmomatic
    for item in rawdir:
        # Filters out the single files by checking for _
        if "_" not in item and item.startswith("SRR"):
            trimmomatic_single = f"java -jar {config.path_to_dirs['trimjar']} SE fastq_files/{item} trimmed_reads/single/{item.split('.')[0]}_single.fastq ILLUMINACLIP:{config.path_to_dirs['adapter']}{config.params['single_trim']} 2>> trimmed_reads/logs/{item.split('.')[0]}_log.file"
            subprocess.call(trimmomatic_single, shell=True)
            print(f"Trimming: {item.split('.')[0]}")
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

        print(f"Trimming: {r1.split('.')[0]} and {r2.split('.')[0]}")
        # Paired trimmomatic
        trimmomatic_paired = f"java -jar {config.path_to_dirs['trimjar']} PE fastq_files/{r1} fastq_files/{r2} trimmed_reads/paired/{r1.split('.')[0]}.trim.fastq trimmed_reads/paired/{r1.split('.')[0]}un.trim.fastq trimmed_reads/paired/{r2.split('.')[0]}.trim.fastq trimmed_reads/paired/{r2.split('.')[0]}un.trim.fastq ILLUMINACLIP:{config.path_to_dirs['adapter']}{config.params['paired_trim']} 2>> trimmed_reads/logs/paired_log.file"
        subprocess.call(trimmomatic_paired, shell=True)


def unzipping():
    """ Unzips GZ files in case those are present. """
    # Makes the working directory to the path it requires.
    os.chdir(config.path_to_files["path_single"])
    single = os.listdir(config.path_to_files["path_single"])
    # Makes sure if there are any files that need unzipping
    if ".gz" in single:
        subprocess.call("gunzip *fastq.gz", shell=True)
    os.chdir(config.path_to_files["path_paired"])
    paired = os.listdir(config.path_to_files["path_paired"])
    if ".gz" in paired:
        subprocess.call("gunzip *fastq.gz", shell=True)
    else:
        # Lets user know that there are no files to unzip
        print("No files to unzip, moving on.")


def bowtiemerge():
    """ Merges paired - unpaired files. """
    os.chdir(config.path_to_files["path_paired"])
    mergedir = os.listdir(config.path_to_files["path_paired"])
    file1 = []
    file2 = []
    # Makes sure the paired dir is not empty
    if len(mergedir) != 0:
        for file in mergedir:
            # Works with the exact output from Trimmomatic
            if "un" in file and "_1" in file:
                file1.append(file)
            if "un" in file and "_2" in file:
                file2.append(file)
    # Sorts both lists so all the _1 and _2 parts match with their counterparts.
    file1.sort()
    file2.sort()

    # Loops through everything one by one so all get merged
    for index in range(0, len(file1)):
        r1 = file1[index]
        r2 = file2[index]
        # Splits on certain parts so the user only sees the SRR number (not required but looks better)
        print(f"Merging: {r1.split('un')[0]} with {r2.split('un')[0]}")
        subprocess.call(f"cat {r1} {r2} 1> {r1.split('_')[0]}_merged_trim.fastq", shell=True)


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

        print(f"Mapping: {r1.split('.')[0]} with {r2.split('.')[0]} and {merged.split('_trim')[0]}")
        # Bowindex path can be changed in config. Bowtie2 is required, Bowtie wont work
        bowtiepaired = f"bowtie2 -p 35 --end-to-end -x {config.path_to_dirs['bowindex']} --fr -1 {r1} -2 {r2} -U {merged} --al-conc {config.path_to_files['mapped_p']}/{r1.split('_')[0]}_%.unmapped.fastq --al {config.path_to_files['mapped_p']}/{r1.split('_')[0]}_unpaired.unmapped.fastq 1> {config.path_to_files['mapped_p']}/{r1.split('_')[0]}.sam 2>> {config.path_to_files['mapped_p']}/unmapped_log.file"

        subprocess.call(bowtiepaired, shell=True)

    # Maps single files
    for item in rawsingledir:
        # Change the workdirectory or Bowtie2 wont find the files.
        os.chdir(config.path_to_files["path_single"])
        if "SRR" in item:
            print(f"Mapping: {item.split('_')[0]}")
            bowtie2single = f"bowtie2 -p 35 --end-to-end -x {config.path_to_dirs['bowindex']} {item} -S {config.path_to_files['mapped_s']}/{item.split('_')[0]}_single_unmapped.fastq 2>> {config.path_to_files['mapped_s']}/mapped_log.file"
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
        # Alerts the user to the start of the assembly process.
        print(f"Assembling: {myfile.split('_')[0]}")
        novosingle = f"mpiexec -n 20 Ray -s {myfile} -o {config.path_to_files['denovo_s']}/{myfile.split('_')[0]}.forblast 1>/dev/null 2>> {config.path_to_files['denovo_s']}/denovo_log_file"
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

        print(f"Assembling: {r1.split('.')[0]} with {r2.split('.')[0]} and {other.split('.')[0]}")
        # Assembles all the files that are paired.
        myray = f"mpiexec -n 20 Ray -p {r1} {r2} -o {config.path_to_files['denovo_p']}/{r1.split('_')[0]}.forblast 1>/dev/null 2>> {config.path_to_files['denovo_p']}/denovo_log_file"
        myrayu = f"mpiexec -n 20 Ray -s {other} -o {config.path_to_files['denovo_p']}/{other.split('_')[0]}.forblast_u 1>/dev/null 2>> {config.path_to_files['denovo_p']}/denovo_log_file"

        subprocess.call(myray, shell=True)
        subprocess.call(myrayu, shell=True)


def contigmerger():
    """ Merges the contigs from the Denovo assemly paired and moves the single ones for blast use."""
    os.chdir(config.path_to_files["denovo_p"])
    mergedir = os.listdir(config.path_to_files["denovo_p"])
    singledir = os.listdir(config.path_to_files["denovo_s"])

    SRR_u = []
    SRR = []
    SRR_S = []

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
        # Merges the contigs from the ray for the paired files.
        contigmerge = f"cat {myfile_u}/Contigs.fasta {myfile_p}/Contigs.fasta > {config.path_to_dirs['workdir']}contigs/{myfile_u.split('.')[0]}.Contigs.fasta"
        subprocess.call(contigmerge, shell=True)

    os.chdir(config.path_to_files["denovo_s"])
    for items in singledir:
        if items.endswith("forblast"):
            SRR_S.append(items)

    SRR_S.sort()

    for index in range(0, len(SRR_S)):
        myfile_s = SRR_S[index]
        # Moves the contigs from the single files
        contig_single = f"cat {myfile_s}/Contigs.fasta > {config.path_to_dirs['workdir']}contigs/{myfile_s.split('.')[0]}.Contigs.fasta"
        subprocess.call(contig_single, shell=True)


def blaster():
    """ Incorporates Diamond blast for all the contigs with the spiked database. """
    os.chdir(config.path_to_files["contigs"])
    contig_dir = os.listdir(config.path_to_files["contigs"])

    all_contigs = []
    for items in contig_dir:
        # Makes sure the contig file is not empty
        if not os.stat(items).st_size == 0:
            all_contigs.append(items)

    for index in range(0, len(all_contigs)):
        mycontig = all_contigs[index]
        print(f"Blasting: {mycontig.split('.')[0]}")
        os.chdir(config.path_to_tools["diamondblast"])
        myblast = f"./diamond blastx -d  {config.path_to_dirs['blastdir']}spiked_database/spikeddb.dmnd -q {config.path_to_files['contigs']}/{mycontig} --sensitive --quiet -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids skingdoms sscinames stitle -k 1 -o {config.path_to_files['blastout']}/{mycontig.split('.')[0]}_matches.tsv"
        subprocess.call(myblast, shell=True)


def blastmatcher():
    """ Checks the blast matches against the known keywords """
    workdir = os.listdir(config.path_to_files['blastout'])
    blastmatches = []
    for match in workdir:
        if match.endswith(".tsv"):
            blastmatches.append(match)

    for index in range(0, len(blastmatches)):
        mymatch = blastmatches[index].split("_")[0]
        # Move to correct location after checking against keywords.
        keyword = f"grep -f {config.path_to_dirs['basepipe']}/resources/Keywords.txt {config.path_to_files['blastout']}/{mymatch}_matches.tsv > {config.path_to_dirs['workdir']}blastoutput/orthomatch/{mymatch}_ortho_matches.tsv"
        subprocess.call(keyword, shell=True)

        os.chdir(f"{config.path_to_dirs['workdir']}blastoutput/orthomatch")
        # Picks the second item out of the file and sees if this is unique and moves it.
        uniquecheck = "awk '{print $2}' " + f"{mymatch}_ortho_matches.tsv | sort | uniq > {config.path_to_files['blastout']}/accessions/{mymatch}_acc_hits.txt"
        subprocess.call(uniquecheck, shell=True)


def esearcher():
    """ Esearches and efetches the accession numbers. """
    os.chdir(f"{config.path_to_files['blastout']}/accessions")
    my_srr = []

    for line in os.listdir(f"{config.path_to_files['blastout']}/accessions"):
        if line.startswith("SRR"):
            my_srr.append(line)

    for srr in my_srr:
        with open(srr) as srr_file:
            for acc in srr_file:
                print(f"Finding sequence for {acc}".strip().split('_')[1])
                seqfetch = f"esearch -db protein -query {acc.strip().split('_')[1]} | efetch -format fasta > fetched/{srr.split('_')[0]}_{acc.strip().split('_')[1]}.fasta"
                subprocess.call(seqfetch, shell=True)



def main():
    """ Main function to call all the other functions. """
    #create_output_dirs()
    #bowtie2_btbuild()
    #downloader()
    #getter()
    #databasedownload()
    #merger()
    # The get_accession is a required argument here.

    #prefetch(get_accession())
    #fastqdump(get_accession())
    #trimmomatic()
    #unzipping()
    #bowtiemerge()
    #bowtie2()
    #denovo()
    #contigmerger()
    #blaster()
    blastmatcher()
    esearcher()


if __name__ == "__main__":
    main()
