#!/bin/bash

# Instructions before executing the pipeline: 
# Make sure the docker containers are installed and the verify the installations
# This pipeline is automated with the dockerized commands and inhouse scripts
# Input: Raw_FastQ files(R1, R2) should be placed in the raw_samples directory. 
# Final_Output: Annotated VCF file. 
# Execution of the pipeline: bash Genepowerx_Bioinformatics_WES_Pipeline.sh <File_Name> <Threads> 

if [[ $# -eq 0 ]]; then
	echo "No arguments Given"
	exit 1
fi

if [ -z "$2" ]; then
        echo "Number of threads not provided, the program will exit as erroneous"
        exit 1
fi

# Directories which contain the files
base_dir="./data"
raw_samples="$base_dir/raw_samples"
ca_raw_samples="$base_dir/ca_raw_samples"
fastqc_dir="$base_dir/fastqc_reports/$1"
trim_params_dir="$ca_raw_samples/trim_params"
trimmed_samples="$base_dir/trimmed_samples/$1"
bowtie_dir="$base_dir/bowtie2/$1"
samtools_dir="$base_dir/samtools/$1"
varscan_dir="$base_dir/varscan/$1"
snpsift_dir="$base_dir/snpsift/$1"
vep_dir="./vep_data"
vep_input="$vep_dir/input"
vep_output="$vep_dir/output"
vcf_files="$base_dir/final_vcf/$1"

echo "Creating Directories, if they already do not exist"
mkdir -p $base_dir
mkdir -p $raw_samples
mkdir -p $trimmed_samples
mkdir -p $bowtie_dir
mkdir -p $samtools_dir
mkdir -p $varscan_dir
mkdir -p $snpsift_dir
mkdir -p $vcf_files
mkdir -p $ca_raw_samples
mkdir -p $trim_params_dir
mkdir -p $fastqc_dir
ls $raw_samples/$1*
mv $raw_samples/$1*R1*.fastq.gz $raw_samples/$1.R1.fastq.gz
mv $raw_samples/$1*R2*.fastq.gz $raw_samples/$1.R2.fastq.gz

#Filenames for sample
echo "Initializing filenames for $1"
raw_file_forward="$raw_samples/$1.R1.fastq.gz"
raw_file_reverse="$raw_samples/$1.R2.fastq.gz"
ca_raw_file_forward="$ca_raw_samples/$1.R1.fastq.gz"
ca_raw_file_reverse="$ca_raw_samples/$1.R2.fastq.gz"
trim_param_file="$trim_params_dir/$1_trim_params.txt"
trimmed_forward="$trimmed_samples/$1_trim_forward.fastq.gz"
trimmed_reverse="$trimmed_samples/$1_trim_reverse.fastq.gz"
untrimmed_forward="$trimmed_samples/$1_untrim_forward.fastq.gz"
untrimmed_reverse="$trimmed_samples/$1_untrim_reverse.fastq.gz"
human_hg38="$base_dir/human_ref/Human_hg38/Human_hg38"
sam_file="$bowtie_dir/$1_aligned.sam"
sorted_bam="$samtools_dir/$1_sorted.bam"
human_hg_fasta="$base_dir/human_ref/Human_hg38/hg38.fa"
mpileup_file="$samtools_dir/$1.mpileup"
snp_vcf_varscan="$varscan_dir/$1_snp.vcf"
indel_vcf_varscan="$varscan_dir/$1_indel.vcf"
snp_vcf_varscan_m="$varscan_dir/$1_snp_modified.vcf"
indel_vcf_varscan_m="$varscan_dir/$1_indel_modified.vcf"
filtered_snp_vcf="$varscan_dir/$1_filtered_snp.vcf"
filtered_indel_vcf="$varscan_dir/$1_filtered_indel.vcf"
filtered_snp_vcf_m="$varscan_dir/$1_filtered_snp_modified.vcf"
filtered_indel_vcf_m="$varscan_dir/$1_filtered_indel_modified.vcf"
annotated_snp_vcf="$snpsift_dir/$1_annotated_snp.vcf"
annotated_indel_vcf="$snpsift_dir/$1_annotated_indel.vcf"
final_vcf="$vcf_files/$1_final.vcf"
vep_input_file="$vep_input/$1.kh_sample.vcf"
vep_output_file="$vep_output/$1.kh_sample.vcf"

# Creating empty files for outputs
touch $ca_raw_file_forward
touch $ca_raw_file_reverse
touch $trim_param_file
touch $trimmed_forward
touch $trimmed_reverse
touch $untrimmed_forward
touch $untrimmed_reverse
touch $sam_file
touch $sorted_bam
touch $mpileup_file
touch $snp_vcf_varscan
touch $indel_vcf_varscan
touch $filtered_snp_vcf
touch $filtered_indel_vcf
touch $annotated_snp_vcf
touch $annotated_indel_vcf
touch $final_vcf
touch $vep_input_file

# Database file paths
database_dir="$vep_dir/databases"
dbsnp_file="$database_dir/dbSNP.vcf.gz"

log_dir="$base_dir/log/$1"
cutadapt_log="$log_dir/cutadapt.log"
fastqc_log="$log_dir/fastqc.log"
trimmomatic_log="$log_dir/trimmomatic.log"
bowtie2_log="$log_dir/bowtie2.log"
sam_bam_conversion_log="$log_dir/sam_bam_conversion.log"
bam_sorting_log="$log_dir/bam_sorting.log"
mpileup_log="$log_dir/mpileup.log"
mpileup_to_snp_log="$log_dir/mpileup_to_snp.log"
mpileup_to_indel_log="$log_dir/mpileup_to_indel.log"
snp_filter_log="$log_dir/snp_filter.log"
indel_filter_log="$log_dir/indel_filter.log"
annotate_snp_log="$log_dir/annotate_snp.log"
annotate_indel_log="$log_dir/annotate_indel.log"
vep_log="$log_dir/vep.log"
timings_log="$log_dir/time_taken.log"

echo "Creating log directories"

mkdir -p $log_dir

echo -e "$(tput setaf 13)Starting GenepoweRx WES pipeline.... $(tput sgr0)";

################ Step1: CutAdapt #################################################################
	
	read1="$folders/*R1_001.fastq.gz"
	read2="$folders/*R2_001.fastq.gz"
	adaptor1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	adaptor2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
	echo -e "$(tput setaf 6)Cutadapt tool is running to adapter removal.... $(tput sgr0)";
	cutadapt_st=$(date +%s)
	docker run -it --rm -v $(pwd):/data/ -w /data/ cutadapt:4.6 --cores $2 -a $adaptor1 -A $adaptor2 -o $ca_raw_samples/$1.R1.fastq.gz -p $ca_raw_samples/$1.R2.fastq.gz $raw_file_forward $raw_file_reverse
	cutadapt_et=$(date +%s)
################ Step2: FastQC ################################################################## 
	
	echo -e "$(tput setaf 6)fastqc processing is starting now.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data/ -w /data/ fastqc:0.12.1 fastqc -t 12  $ca_raw_file_forward $ca_raw_file_reverse -o $fastqc_dir
	fastqc_et=$(date +%s)

	cd $fastqc_dir
	unzip $1.R1_fastqc.zip 
	unzip $1.R2_fastqc.zip 
	cd ../../..
         ########### Extracting QC Information from the FastQC Reports #############
	R1_Reads=$(grep -w 'Total Sequences' $fastqc_dir/$1.R1_fastqc/fastqc_data.txt | awk -F " " '{print $3}') 
	R2_Reads=$(grep -w 'Total Sequences' $fastqc_dir/$1.R2_fastqc/fastqc_data.txt | awk -F " " '{print $3}')
	R1_R2=`echo "x=${R1_Reads}; y=${R2_Reads};x+y" |bc`
	Depth=`echo "x=${R1_R2}*151; y=36035818; x/y" | bc`
	echo "R1_Reads are: $R1_Reads"
	echo "R2_reads are: $R2_Reads"
	echo "R1+R2_Reads are: $R1_R2"
	echo "Depth value is: $Depth"

#########################################################################################################
	
	#Getting trim parameters - Python program
	python3 get_trim_value.py $fastqc_dir/$1.R1_fastqc/fastqc_data.txt $fastqc_dir/$1.R2_fastqc/fastqc_data.txt > $trim_param_file
	trim_param=$(<$trim_param_file)
	echo "trimming parameters are: $trim_param"
	HEADCROP=$(echo $trim_param | grep -E "HEADCROP:" | awk -F " " '{print $1}')
	### If conditions checks whether trim_parameter has only headcrop / both values
	### If crop value has the trim_parameter then it will subtract the trim value from 151
	if echo "$trim_param" | grep -q " CROP:" ; then
			CROP_Value=$(echo $trim_param | awk -F " " '{print $2}' | awk -F ":" '{print $2}')
			CROP=$(echo "151 - $CROP_Value" | bc)
			CROP1=" CROP:$CROP"
			trim_param1="$HEADCROP$CROP1"
	else
			trim_param1="$HEADCROP"
	fi
	echo $trim_param1
	
################ Step3: Trimmomatic ##################################################################
	
	if [ -n "$2" ]; then
		echo -e "$(tput setaf 6)Calling Trimmomatic to trim the raw files for better quality.... $(tput sgr0)";
		trimmo_st=$(date +%s)
		docker run -it --rm -v $(pwd):/data/ -w /data/ trimmomatic:0.39 PE -threads $2 $ca_raw_file_forward $ca_raw_file_reverse $trimmed_forward $untrimmed_forward $trimmed_reverse $untrimmed_reverse $trim_param1
	else
		echo "Trimming not being performed"
		trimmed_forward="$ca_raw_file_forward"
		trimmed_reverse="$ca_raw_file_reverse"
		echo "trimmed forward  is $trimmed_forward and \n trimmed reverse is $trimmed_reverse \n these are set for bowtie2"
	fi
	trimmo_et=$(date +%s)
	
################ Step4: Bowtie2 ##################################################################
	
	if [ -n "$2" ]; then
		echo -e "$(tput setaf 6)Bowtie2 Alignment starting in a moment.... $(tput sgr0)";
		docker run -it --rm -v $(pwd):/data/ -w /data/ bowtie2:2.5.3 -p $2 -x $human_hg38 -1 $trimmed_forward -2 $trimmed_reverse -S $sam_file 2> $bowtie2_log
	fi
	bow_et=$(date +%s)
	
################ Step5: Samtools ##################################################################
	
	bow_et_et=$(date +%s)
	echo -e "$(tput setaf 6)sorting the sam file.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data/ -w /data/ samtools:1.19 samtools sort --threads $2 $sam_file --output-fmt BAM -o $sorted_bam 
	sort_sam_et=$(date +%s)
	echo -e "$(tput setaf 6)Generating mpileup file.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data/ -w /data/ samtools:1.19 samtools mpileup -E -f $human_hg_fasta $sorted_bam -o $mpileup_file 
	mpileup_et=$(date +%s)
	
############### Step6: Varscan ##################################################################
	
	echo -e "$(tput setaf 6)mpileup to snp.... $(tput sgr0)";
	# sed -i 1d $mpileup_file
	docker run -it --rm -v $(pwd):/data -w /data varscan:2.4.6 mpileup2snp $mpileup_file --output-vcf 1 > $snp_vcf_varscan
	mpileup_snp_et=$(date +%s)
	echo -e "$(tput setaf 6)mpileup to indel.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data -w /data varscan:2.4.6 mpileup2indel $mpileup_file --output-vcf 1 > $indel_vcf_varscan
	mpileup_indel_et=$(date +%s)
	head -n 8 $snp_vcf_varscan >> $mpileup_to_snp_log; tail -n 4 $snp_vcf_varscan >> $mpileup_to_snp_log
	head -n 8 $indel_vcf_varscan >> $mpileup_to_indel_log; tail -n 4 $indel_vcf_varscan >> $mpileup_to_indel_log
	grep -E '^#|^chr' $snp_vcf_varscan > $snp_vcf_varscan_m 
	grep -E '^#|^chr' $indel_vcf_varscan > $indel_vcf_varscan_m
	# cat $snp_vcf_varscan | tail -n +9 | head -n -4 > $snp_vcf_varscan_m 
	# cat $indel_vcf_varscan | tail -n +9 | head -n -4 > $indel_vcf_varscan_m  
	echo -e "$(tput setaf 6)Filtering snps.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data -w /data varscan:2.4.6 filter $snp_vcf_varscan_m --indel-file $indel_vcf_varscan > $filtered_snp_vcf
	filter_snp_et=$(date +%s)
	echo -e "$(tput setaf 6)Filtering indels.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data -w /data varscan:2.4.6 filter $indel_vcf_varscan_m > $filtered_indel_vcf
	filter_indel_et=$(date +%s)
	grep -E '^#|^chr' $filtered_snp_vcf > $filtered_snp_vcf_m
	grep -E '^#|^chr' $filtered_indel_vcf > $filtered_indel_vcf_m
 
################ Step7: SnpSift ##################################################################
	
	echo -e "$(tput setaf 6)Annotating with dbsnp for snp and indel files.... $(tput sgr0)";
	docker run -it --rm -v $(pwd):/data/ -w /data/ snpsift:5.1 annotate $dbsnp_file $filtered_snp_vcf_m > $annotated_snp_vcf
	annotate_snp_et=$(date +%s)
	docker run -it --rm -v $(pwd):/data/ -w /data/ snpsift:5.1 annotate $dbsnp_file $filtered_indel_vcf_m > $annotated_indel_vcf
	annotate_indel_et=$(date +%s)

################ Step8: Ensembl-VEP ##################################################################

	annotation_snpsift_et=$(date +%s)
	export SAMPLE="$1"
	cp $annotated_snp_vcf $vep_input_file
	echo -e "$(tput setaf 6)vep for annotating with clinvar.... $(tput sgr0)";
	docker run -it -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep --format vcf --vcf --cache --offline --force_overwrite --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/ --input_file /opt/vep/.vep/input/$SAMPLE.kh_sample.vcf --output_file /opt/vep/.vep/output/$SAMPLE.kh_sample_final.vcf --no_stats --fork $2 --everything --custom /opt/vep/.vep/databases/clinvar_20240317.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN
	mv $vep_output/$SAMPLE.kh_sample_final.vcf $final_vcf
	vep_et=$(date +%s)

######################################################################################################

	# TODO Check the time elapsed for each step and the total pipeline
	echo "Cutadapt took :  $(( $cutadapt_et - $cutadapt_st )) " > $timings_log
	echo "FastQC took :  $(( $fastqc_et - $cutadapt_et )) " >> $timings_log
	echo "Trimmomatic took : $(( $trimmo_et - $trimmo_st )) " >> $timings_log
	echo "Bowtie2 took : $(( $bow_et - $trimmo_et )) " >> $timings_log
	echo "sorting the sam file took : $(( $sort_sam_et - $bow_et )) " >> $timings_log
	echo "Generating mpileup took : $(( $mpileup_et - $sort_bam_et )) " >> $timings_log
	echo "mpileup to snp  took : $(( $mpileup_snp_et - $mpileup_et )) " >> $timings_log
	echo "mpileup to indel  took : $(( $mpileup_indel_et - $mpileup_snp_et )) " >> $timings_log
	echo "Filtering snps took : $(( $filter_snp_et - $mpileup_indel_et )) " >> $timings_log
	echo "Filtering indel took : $(( $filter_indel_et - $filter_snp_et )) " >> $timings_log
	echo "annotating snps and indels took : $(( $annotation_snpsift_et - $filter_indel_et )) " >> $timings_log
	echo "vep to generate vcf file took : $(( $vep_et - $annotation_snpsift_et )) " >> $timings_log
	echo "Total time till now for partially completed pipeline is $(( $vep_et - $cutadapt_st )) " >> $timings_log
	echo "The pipeline is complete for sample $1"
	echo 
	echo -e "$(tput setaf 13)Thanks for using GenepoweRx WES pipeline. Have a good Day!!!.... $(tput sgr0)";
	echo

## After generation of the final VCF the info and depth columns splitting.  
python3 INFO_splitting_VEP.py ${1}
echo "Info column splitting has been done and saved the output file for the sample - ${1}"

###########################################################################################
This block of commands will change the ownership of the directories...
username=$(whoami)
hostname=$(hostname)
sudo chown ${username}:${hostname} ./data/
sudo chown ${username}:${hostname} ./data/*
sudo chown ${username}:${hostname} ./data/*/*
sudo chown ${username}:${hostname} ./data/*/*/*
sudo chown ${username}:${hostname} ./data/*/*/*/*
