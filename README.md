# Anopheles mosquitos #
This projects main focus lies with the research into the anopheles mosquito, more specifically which viruses they could be carrying with them. This research into the mosquito and its possible viruses has been done with data collected from the NCBI website and in coordination with the NCBI databases in order to screen for endogenous viral elements (EVEs)

## Table of content

- [Project description](#project-description)
- [Installation](#installation)
    * [Prerequisites](#prerequisites)
    * [Packages](#packages)
    * [Scripts and pipeline](#scripts-and-pipeline)
- [Contact](#contact)

## Project description
One of the main goals of this project will focus on building/updating an already existing pipeline and make it so that it is able to screen anopheles data and map which viruses could be within the mosquito. This will also look for EVEs, but this is done in an already pre-existing pipeline.

Another goal for this project is that the mosquito data will be run through the pipeline and analysed, this will be done in order to see if any viruses are present within the mosquito and which viruses they may be

The final product that will reside within this repository will be the pipeline that can analyse the anopheles mosquito data and gives a representation in the form of a dag file in order to show the user what happened.


## Installation
For the complete repository please use the following piece of code inside the terminal:
```git clone https://github.com/jlalisan/Anopheles```

In order to run the project, please install the packages that are required and check their version, it may occur that a different version of a package does not have the same functions.

Note that the databases are not present within this repository, these can be manually chosen. A NCBI database can be chosen to use or a spiked database can be made. An installation guide for the NCBI database will be shown below.

### Prerequisites
* Python (3.7 or higher)
* Snakemake (7.14 or higher)
* SRA-toolkit
* Bowtie2
* Ray-tools
* Trimmomatic (0.39)
* Diamond blast

Snakemake can be installed via a commandline input, this can be done for any engine, this does not have to be conda. ```conda install -c bioconda snakemake``` Furthermore an complete installation guide van be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

The SRA-toolkit can be cloned from the github via [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) here a full installation guide is also shown for the first time user of this program. Do note if there are multiple users and for the first time, the following argument needs to be used. 
```export PATH=$PATH:/path/to/sratoolkit/bin```
This needs to be done if the program has trouble finding the prefetch or the fasterq-dump application needed.

Bowtie2 can be installed via the commandline if this is not yet installed. Do note that Bowtie2 is *not* available for windows. The command to install this is: ```conda install -c bioconda bowtie2``` Mamba or another engine can also be used for the installation process.

Ray-tools are used for the assembly process within the pipeline and can be found [here](https://github.com/sebhtml/Ray-Releases) do note that the latest release was used, and another version could potentially not work.

For the trimming process the trimmomatic 0.39 jar is required along with the adapter sequences unless the user has their own adapters. The jar can be downloaded [here](http://www.usadellab.org/cms/?page=trimmomatic) Another version or an later version can be used but are not garuanteed to give the same or an correct answer.

Diamond blast is used to see if there are viruses within the data and for a full guide and download information the [github link](https://github.com/bbuchfink/diamond) or directly via [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) Any required database can also be downloaded via this same link.

### Packages
|Name                                   |Version              |   
|---                                    |---                  |
|os                                     |Varies per device    |
|re                                     |2022.6.2             |
|subprocess                             |3.7                  |

The packeges used within this pipeline are used with a python version 3.7, any version above 3.7 should also work. Do note that these packages are for the snakemake scripts and not the python scripts within the repository. 

### Scripts and pipeline
The pipeline consists out of three main scripts. Pre-processing, processing and post-processing. These are all included within Main.smk. To call the scripts enter the following ```snakemake --snakefile main.smk --jobs --keep-going``` with this command it is ensured the entire process finishes even if a file cannot be found due to an human error with the input.

## Contact

* Lisan Eisinga
  * j.l.a.eisinga@st.hanze.nl 
