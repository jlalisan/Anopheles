## Pipeline ReadMe
For the pipeline there are a lot of rules that are included, in order to make it more clear this separate readme has been written to include all possible useful parameters and some extra information on the rules.

## Table of content

- [Rules explained](#rules-explained)
    * [Prefetching and Fasterq dumping](#prefetching-and-fasterq-dumping)
    * [Trimmomatic](#trimmomatic)
    * [Bowtie2](#bowtie2)
    * [Denovo](#denovo)
    * [Blasting](#blasting)
    * [Efetching](#efetching)
    * [Sam and Bamtools](#sam-and-bamtools)
    * [Lofreq](#lofreq)
- [Side notes](#side-notes)
- [Contact](#contact)

## Rules explained
For all the rules there are parameters or 'need to know' facts that can help the user run the pipeline more efficiently, please follow the table of content if you are looking for a specific tool that is used in a rule. Some of this will overlap with information inside the repository ReadMe.

### Prefetching and Fasterq dumping
For the fasterq-dump and prefetching the file names need to be inside the config.yaml file, all these files will be downloaded. Note that even if the names of the files are incorrect the pipeline will still download them or attempt to download them, with files that do not exist this will throw an error. There are no parameters that need to be filled out for the fasterq-dumping or prefetching process. The table below however shows any parameters that could be added to either in case it is required. For the prefetching stage only one parameters is of significance which is the `--max-size` which can be set in order to force the pipeline to download only files up to the max size.

|Fasterq-dump parameters |Explanation                                   |
|---                     |---                                           |
|--gzip                  |Gzips all the files that are downloaded       |
|--split-3               |Splits the files in case of paired files      |
|--concatenated          |Concatenates the files.                       |
|--include-technical     |Includes all reads. --exclude is the opposite.|

For the prefetching and fasterq-dumping the SRA toolkit is required as it states in the repositories ReadMe, be sure to always run the ```export PATH=$PATH:/path/to/sratoolkit/bin``` command in the command line. If this command is forgotten the pipeline might not be able to find the sratools and throw an error. If the SRA toolkit is not yet installed on the device that you are using then please [click here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to download the latest version, the Github repository here also has an guide on how to [install this](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration). After this is done the ```prefetch``` command should be tested out to see if the command can be found, if this is not the case please repeat the steps within the Github tutorial, note that on windows it has been noticed that it does not always work well.

### Trimmomatic
The Trimmomatic process is done via the Trimmomatic-0.39 jar, if the jar is not downloaded on the device or installed in the incorrect place or in case of an incomplete directory it can be found [here](http://www.usadellab.org/cms/?page=trimmomatic). The Trimmomatic command has several parameters that could be added to the statements if it is necessary. The table below shows some of the Trimmomatic parameters.

|Trimmomatic parameters |Explanation                                   |
|---                    |---                                           |
|ADAPTER                | Trimming of adapter sequences from the reads |
|MINLEN                 | Minimum length of reads after clipping       |
|SLIDINGWINDOW          | Trimming of low quality reads                |
|TRAILING               | Trimming of low quality bases from the 3' ends of the reads                                                           |
|LEADING                | Trimming of low quality bases from the 5' ends of the reads                                                           |

The input of the Trimmomatic is not to be changed since these are the files that are downloaded by the Fasterq-dump. The PHRED encoding can be changed to 64 if needed or removed entirely if all the files have a defined PHRED score. If the PHRED score is not known and not defined in the parameters the file will not be trimmed. The parameters for Trimmomatic are explained more [here](https://wiki.bits.vib.be/index.php/Parameters_of_Trimmomatic)

### Bowtie2
The Bowtie2 is a mapping tool that can be installed via the command line on linux or on mac, this is not a tool that can be used for Windows and will throw an error upon installation. Installation can be achieved with the ```conda install -c bioconda Bowtie2``` command if Bioconda is not on the device please install that first.

For the Bowtie2 tool there are a lot of parameters that could be used when opening up the code, the most important ones that are used are listed below.

|Bowtie2 parameters|Explanation                                         |
|---               |---                                                 |
|--sensitive       | Slows down the mapping to be more sensitive.       |
|-x                | The base name of the index for the reference genome|
|-S                | File to write SAM alignments                       |

There are of course a lot more parameters available, however not all are required, the ones above in table three are used for their usefulness, others can be found in the [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

The Bowtie2 mapping makes use of the Bowtie2 indexing, this indexes the reference genome that is noted within the config.yaml file. With the indexing there are six different index files that are created, these are all used with the mapping.

### Denovo
The Denovo assembly is facilitated by RAY which uses MPI to run the assembly. Note that a C++ compiler is needed in order to run Ray for this project g++ was used which can be obtained with the following command ```sudo apt install g++```.
The MPI (Message passing interface) required can be downloaded anywhere, the one used for this project was OpenMPI which was downloaded [here](https://www.open-mpi.org/). The Ray.tar was downloaded from the [Github page](https://github.com/sebhtml/Ray-Releases).

Ray has the -n option which amounts to the amount of cores the program uses, for this project the cores have been put on a low number of two since there were multiple devices testing the code. However if there are more cores available then more cores are recommended. Do note that Ray is not perfect and if it runs too fast then some files do not get contigs. for a complete overview of ray their website has a [manual](https://denovoassembler.sourceforge.net/manual.html) which contains all the parameters.

### Blasting
The Blast was done with a Diamond blast database, since this is quicker then a normal blast. The NR database was used for this project because it was available. The pipeline has also been tested with a spiked (24.1) Uniprot database (Uniprot 90), the Uniprot database can be found [here](https://www.uniprot.org/help/downloads), and the virus database can be downloaded [here](https://rvdb.dbi.udel.edu/). A reference genome can be picked after installing diamond, which can be done by [downloading](https://github.com/bbuchfink/diamond) the software, after which the database can be build, since the taxonomy data is required an argument similar to ```	diamond makedb --in nr.fa -d database --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp``` needs to be used. With this the database file will be build resulting in a 'database.dmnd' file that can be inserted in the config file.

|Blast parameters   | Explanation                                      | 
|---               |---                                               |
|sensitive         | Makes sure the blast is sensitive                |
|threads           | The threads chosen so the blast wont get killed  |
|quiet             | Disable all terminal output                      |
|qseqid            | Query Sequence - id                              |
|sseqid            | Subject Sequence - id                            |
|pident            | Percentage of identical matches                  |
|length            | Alignment length                                 |
|mismatch          | Number of mismatches                             |
|gapopen           | Number of gap openings                           |
|qstart            | Start of alignment in query                      |
|qend              | End of alignment in query                        |
|sstart            | Start of alignment in subject                    |
|send              | End of alignment in subject                      |
|evalue            | Expect value                                     |
|bitscore          | Bit score                                        |
|staxids           | Unique Subject Taxonomy ID                       |
|skingdoms         | Unique Subject Kingdom                           |
|sscinames         | Unique Subject Scientific Name                   |

For a complete overview and an extra explanation of the diamond Blast parameters go to the [Github page](https://github.com/bbuchfink/diamond), here a full explanation is given including parameters.

### Efetching
The efetching is done by making use of the bio.entrez package that can be viewed [here](https://biopython.org/docs/1.76/api/Bio.Entrez.html) this can be installed with the command line by using the ```conda install -c bioconda entrez-direct``` command. With this the Efetch and Esearch options open up. With the Efetch there can be a choice made what data needs to be fetched, this can be protein if that is required for the scope of this pipeline it is put on 'fasta_cds_na' for the nucleotide sequence. Also the databases can be changed from protein, this is however not recommended.

### Sam and Bamtools
SAM, BAM, and Bedtools are used to make the final VCF file and earlier in the pipeline as well, in order for these arguments to function properly the correct packages need to be installed. For the SAM and BAM tools HTSlib needs to be [downloaded](http://www.htslib.org/download/) this contains both packages and has an user guide on the site. For the bedtools, ```apt-get install bedtools``` is sufficient. The parameters or arguments used should not be changed.

### Lofreq
Lofreq is a fast and sensitive variant-caller for inferring SNVs and indels from next-generation sequencing data. It makes full use of base-call qualities and other sources of errors inherent in sequencing, which are usually ignored by other methods or only used for filtering and can be obtained on their [website](https://csb5.github.io/lofreq/). No arguments need to be changed here.

## Side notes
For all the software and programs the latest versions should be used unless it is specified in the doc of the program that a needed requirement has been removed. The software all runs independently of each other but it is recommended that it is put in the same directory for example: 'Tools'. 

## Contact
* Lisan Eisinga
  * j.l.a.eisinga@st.hanze.nl 

