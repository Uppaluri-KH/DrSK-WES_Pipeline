The readme file includes the list of tools and their installation steps with instructions were given here. 




This section includes the steps of installing the K&H_WES pipeline with docker 
Installing Docker: https://docs.docker.com/engine/install/
Installing Docker-compose: https://docs.docker.com/compose/install/
These manuals has the installation and execution commands with different operating systems.
After the installation is successful with docker and docker-compose, the next step is to create dockerfiles for the tools used in this pipeline with required packages for successful executions i.e., FastQC, Cutadapt, Trimmomatic, Bowtie2, Samtools, VarScan2, SnpSift, and VEP. 
Then, create a docker-compose.yml script to build all these dockers into our local machine. 
Now, once the docker images are created, check whether the dockers are running or not. 

##### Splitting the depth columns from the FORMAT and SAMPLE columns for Down stream analysis ###########################################

1. The raw vcf file from the wes pipe line had FORMAT and SAMPLE columns for depth ("**GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR**")information.
2. Splitted these key value pairs and created new depth columns(**GT,GQ,SDP,DP,RD,AD,FREQ,PVAL,RDF,RDR,ADF,ADR**).
3. Exported the vcf file along with the depth column as tsv file for Downstrem Analysis.
