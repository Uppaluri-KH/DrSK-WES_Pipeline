The readme file includes the list of tools and their installation steps with instructions were given here. 




This section includes the steps of installing the K&H_WES pipeline with docker 
Installing Docker: https://docs.docker.com/engine/install/
Installing Docker-compose: https://docs.docker.com/compose/install/
These manuals has the installation and execution commands with different operating systems.
After the installation is successful with docker and docker-compose, the next step is to create dockerfiles for the tools used in this pipeline with required packages for successful executions i.e., FastQC, Cutadapt, Trimmomatic, Bowtie2, Samtools, VarScan2, SnpSift, and VEP. 
Then, create a docker-compose.yml script to build all these dockers into our local machine. 
Now, once the docker images are created, check whether the dockers are running or not. 


**######################## Splitting the depth columns from the FORMAT and SAMPLE columns for Down stream analysis ###########################################
**
1. 
