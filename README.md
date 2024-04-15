**GenepowerX Bioinformatics pipeline for Whole Exome Sequencing (WES) data analysis.** 

The pipeline follows the pre-processing of the raw data, variants and InDels identification followed by VEP annotations, Tabular separation of VCF file. The pipeline is tested under Ubuntu 21.04 LTE, CentOS Linux 8.0. Part of the Python script for visualizing Insert size frequency and shell scripts for data selection are also included.

**Prerequisites**
Linux/CentOS/ operating system
This pipeline is built with Docker
For optimal computational performance, it is recommended to use multi threading (Multiple cores) with sufficient Memory (Minimum 32 RAM is recommended), and 500Gb storage is required. 
This repository has the next generation of the Whole Exome Sequencing (WES) analysis from raw data fastq-files to vcf generation , annotation, and filtration for reporting the variants based on condition specific. The contents of this repository are 100% open source. 

**Installation**:

Manual installations using source codes or installing binaries for FastQC, Cutadapt, Trimmmomatic, Bowtie2, Samtools, VarScan, snpsift, VEP. The installation instructions were provided in the installations.txt file.  

**Dependencies**

The automated pipeline consists of multiple open-source tools with their dependencies so docker installation is recommended.
Python 3 with libraries like numpy and pandas needs to be installed. 

**Dockerfiles creation and Installations:**

We have created the Dockerfiles for each tool with the latest versions and these dockerfiles can be found in the “Dockerfiles.zip”. We recommend using them once you are done with downloading and extraction. Instructions were provided in the “Docker_installations.sh” file.   

**Pipeline Execution:**

bash Genepowerx_Bioinformatics_WES_Pipeline.sh <samples name/bar code> <Number of Threads>
Example: “bash Genepowerx_Bioinformatics_WES_Pipeline.sh Genepowerx_tutorial 24”
If file names have Genepowerx_tutorial.R1.fastq.gz & Genepowerx_tutorial.R2.fastq.gz script requires the first part of the file name for creating directories and output file names for each sample.

**Data Availability:**

The data requirements and test data along with detailed installations are available in the githiub repository. 

**Directories creation and prerequisites:** 

After installation of Dockers and successful validation, the pipeline is allowed to run. But we recommend creating a mother/base directory named as "pipeline" in the HOME directory where the system has enough storage > 300Gb for a single WES sample to complete the run.
Once you are done with creation of pipeline directory, download the "Genepowerx_Bioinformatics_WES_Pipeline.sh" script into it along with “get_trim_value.py” (This python script generates the head and tail cropping parameters for trimmomatic from fastqc result files)
python script. Next we recommend creating a sub directory in the pipeline folder named as "data". And within the data directory create a sub-directory named "raw_samples". Once we are done with successful creation of these directories, we recommend copying the paired-end raw data fastq files into this directory. Now, setting up pre-requisites for initiating the pipeline like alignment index, reference genome fasta, clinvar, dbsnp files.

These are optional, because we have already downloaded these files and uploaded them to github/dropbox/google_drive. Download them and use them for now.
 
**Download hg38 fasta:** 

wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -k hg38.fa.gz

**Download ClinVar:**

https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

**Download dbSNP:**

https://ftp.ncbi.nih.gov/snp/latest_release/VCF/

Build index files for the hg38 genome assembly using bowtie2 build
https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer

**Raw data for tutorial**

Forward file
https://www.dropbox.com/scl/fi/27azwd52f4n5tf4d3tt73/Genepowerx_tutorial.R1.fastq.gz?rlkey=o3ymovaklcsx0xb08fp50fdjw&dl=0

Click on the above link >>> File >>> Download >>> Or Continue with download only

Reverse file
https://www.dropbox.com/scl/fi/xj8jvb5lnublgchjc0d6j/Genepowerx_tutorial.R2.fastq.gz?rlkey=o2u0sxh6f0g837fvlsoiqxvw1&dl=0

Click on the above link >>> File >>> Download >>> Or Continue with download only

After downloading these prerequisites, open the terminal from the pipeline directory and initiate the pipeline as mentioned in the section “Pipeline Execution”. Then, the pipeline creates all the sub directories for each step along with result files. This pipeline has been tested on more than 3000 samples with xeon processors (Server configurations were mentioned in the main paper). The run time was ~90 minutes for 50 million paired-end reads. In addition to that we also included the log file creation where the user can find the run time in seconds for each tool. After completion of the entire pipeline we will get the total run time for all the tools. Once you are done with successful completion of running the pipeline, the end result will be generated in the “final_vcf” directory with sample/bar code name along with annotations.
 
Here. We have added additional information like splitting the VCF file and converting the vcf file format into tabular separated file. The “INFO_splitting_VEP.py” considers the annotated final vcf file and creates an excel file with column wise splitted Quality Check (QC) parameters i.e., read depth, variant allele fraction etc. In addition to that the annotated information also splitted from vcf INFO column i.e., variant associated gene name, variant consequences, clinical condition, Intron/Exon number, population based allele frequencies, etc. The script can be modified for specific field requirements accordingly. 
 
**Troubleshooting:**

1.  The file permissions need to be changed once the end results are done.
Ex: sudo chown user:user <directory/file name>. This command will change the ownership of the file which allows accessing the files. Note: We have included this in the main script. 

2.  Mounting the volumes to the docker containers in the pipeline is a validation check.
Usually, most of the errors thrown out are “No such file or directory”. 
Solution: Mount the directory using -v option followed by the path to the directory on the host machine and the path where you want to mount it in the container. 
Here's the basic syntax: 
docker run -v /path/on/host:/path/in/container <image_name> <command> 
Replace /path/on/host with the directory path on the host machine that you want to mount, /path/in/container with the directory path where you want to mount it inside the container, image_name with the name of the Docker image you want to run, and command with the command you want to execute inside the container.

3.  Sometimes the log files can be created inside the output files which leads to errors in generating the final output. 
Solution: In such cases the script needs to be modified accordingly. Ex: VarScan mpileup output file. We have modified the file accordingly using shell commands. Please refer to the pipeline execution script for more details. 

4.  Running VEP - MSG: ERROR: Multiple assemblies found for cache version 
Solution: If there are multiple cache versions of an assembly, Specify the assembly to be used or move the unused assembly from the directory.  (“--assembly [assembly]”)

