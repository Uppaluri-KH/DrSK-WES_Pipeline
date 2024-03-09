Download and install software and tools

This section includes the steps of installing the K&H_WES pipeline with docker 
Installing Docker: https://docs.docker.com/engine/install/
Installing Docker-compose: https://docs.docker.com/compose/install/
These manuals has the installation and execution commands with different operating systems.
After the installation is successful with docker and docker-compose, the next step is to create dockerfiles for the tools used in this pipeline with required packages for successful executions i.e., FastQC, Cutadapt, Trimmomatic, Bowtie2, Samtools, VarScan2, SnpSift, and VEP. 
Then, create a docker-compose.yml script to build all these dockers into our local machine. 
Now, once the docker images are created, check whether the dockers are running or not. 


Here's the docker execution commands used in the pipeline. 

1. FastQC: docker run --rm -v <path to mount the directory> -w <working_directory> <fastqc_docker_image_name> -t <no.of threads to use> <input_fastq_file> -o <output_directory>

2. Cutadapt: docker run --rm -v <path to mount the directory> -w <working_directory> <cutadapt_docker_image_name> --cores <no.of threads to use> -a <Adapter_Sequences> -A <Adapter_Sequences> -o <output_file.R1.fastq.gz> -p <output_file.R2.fastq.gz> <input_file.R1.fastq.gz> <input_file.R2.fastq.gz>

2a. InHouse scripts to get the trimming parameters from the results of fastqc for trimmomatic execution. 
3. Trimmomatic: docker run --rm -v <path to mount the directory> -w <working_directory> <trimmomatic_docker_image_name> trimmomatic PE -threads <no.of threads to use> <forward_reads.fastq> <reverse_reads.fastq> <forward_trim.fq> <forward_untrim.fq> <reverse_trim.fq> <reverse_untrim.fq> <trimming parameters for HEADCROP, CROP>

4.Bowtie2: docker run --rm -v <path to mount the directory> -w <working_directory> <bowtie2_docker_image_name> bowtie2-build <reference_sequence.fasta> <index_name>
docker run --rm -v <path to mount the directory> -w <working_directory> <bowtie2_docker_image_name> bowtie2 -p (No. of threads) -x <index_name> -1 <trimmed_forward_file.fastq> -2 <trimmed_reverse_file.fastq> -S <aligned_output.sam>

5. Samtools: SAM2BAM, SORT, MPILEUP
docker run --rm -v <path to mount the directory> -w <working_directory> <samtools_docker_image_name> samtools view -bS <aligned_file.sam> > <output_file.bam>
docker run --rm -v <path to mount the directory> -w <working_directory> <samtools_docker_image_name> samtools sort <bam_file> > <sorted_bam> 
docker run --rm -v <path to mount the directory> -w <working_directory> <samtools_docker_image_name> samtools mpileup -E -f <human_ref.fa> <sorted_bam> > <mpileup_file>

6. Varscan2: SNPs, Indels calling & Filtering SNPs, Indels
docker run --rm -v <path to mount the directory> -w <working_directory> <varscan_docker_image_name> mpileup2snp <mpileup_file> --output-vcf 1 > <snp_output.vcf> 
docker run --rm -v <path to mount the directory> -w <working_directory> <varscan_docker_image_name mpileup2indel <mpileup_file> --output-vcf 1 > <indel_output.vcf> 
docker run --rm -v <path to mount the directory> -w <working_directory> <varscan_docker_image_name> filter <snp_file.vcf> --indel-file <indel_file.vcf> > <filtered_snp_vcf> 
docker run --rm -v <path to mount the directory> -w <working_directory> <varscan_docker_image_name> filter <indel_vcf_varscan> > <filtered_indel_vcf> 

7. SnpSift: Annotating SNPs & Indels
docker run --rm -v <path to mount the directory> -w <working_directory> <snpsift_docker_image_name> annotate <database_file> <filtered_snp_vcf> > <annotated_snp_vcf>
docker run --rm -i -v <path to mount the directory> -w <working_directory> <snpsift_docker_image_name> annotate <dbsnp_file> <filtered_indel_vcf> > <annotated_indel_vcf>

8. Ensembl-VEP: 
docker run -i -v -v <path to mount the directory> <vep_docker_image_name> vep --format vcf --vcf --cache --offline --force_overwrite --dir_cache <cache_directory> --dir_plugins <Plugins_directory> --input_file <input_file.vcf> --output_file <output_vep.vcf> --no_stats --fork <no.of threads> --everything --custom <customised database annotations with required fields> 

Resources download
1. Reference genome: wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
unzip hg38.fa

2. Databases:
CliVar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
dbSNP: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/
1KG: https://www.internationalgenome.org/data
gnomAD: https://gnomad.broadinstitute.org/downloads
CADD: https://cadd.gs.washington.edu/download
dbNSFP: http://database.liulab.science/dbNSFP#database
OMIM: https://www.omim.org/


Troubleshooting
1. Most of the times problems arises while installing the tools such as dependency errors. So, While creating dockerfiles make sure that to keep the commands of required dependencies then the docker image installations and exectution of the tool is successful without any possible errors. 

2. Trimmomatic execution errors: Trimmomatic execution mainly focusses on eliminating the noise/bad quality reads. There are more possible chances of removing good quality reads. So consideration of trimming parameters i.e., soft, hard trimming, headcrop, crop, minimum length, sliding window, average quality etc are to be prepared wisely. Always note that the trimmomatic works in a sequential manner based on the provided trimming parameters in the command.

3. Adapter removal is necessary before proceeding to alignment otherwise it will mislead with false alignments and it will misguide the variants detected by the variant caller.

4. Variants called from the selective variant calling algorithm are to be annotated with latest databases information, as the variant information keeps updating with databases like ClinVar, dbSNP and Ensembl-VEP. 

5. For low quality read depth samples, some of the parameters needs to be tweaked/modified accordingly in variant calling steps such as vaf, maf, read depth etc., Otherwise setting a threshold of read coverage, total number of reads, alignemnt percentage for each sample is a very good move before initiating the WES pipeline.

6. Always keep on a track that removing the cacahe/temporary files after exectuting the pipeline with multiple samples. 
Because it stores the cache information and there are possibility of throwing an error like broken pipe/unsuccessful pipeline execution. 


