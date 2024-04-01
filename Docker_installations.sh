## Installing Docker engine to local machine
# sudo apt install docker.io

## Checking docker is running actively or not. If not enable it
sudo systemctl status docker 
# sudo systemctl enable docker
sudo -S chmod 666 /var/run/docker.sock


### Creating Dockerfiles for each tool
We have created dockerfiles for each tool. These docker files can be found the Docker directory. 

### Building the docker images from the respective dockerfiles.
docker build -t cutadapt:4.6 ./cutadapt/
docker build -t fastqc:0.12.1 ./Docker_files/fastqc/
docker build -t trimmomatic:0.39 ./Docker_files/trimmomatic/
docker build -t bowtie2:2.5.3 ./Docker_files/bowtie2/
docker build -t samtools:1.19 ./Docker_files/samtools/
docker build -t varscan:2.4.6 ./Docker_files/varscan/
docker build -t snpsift:5.1 ./Docker_files/snpsift/

### ENSEMBL-VEP Dockerization. 
mkdir $HOME/vep_data
sudo chmod 777 -R $HOME/vep_data
docker pull ensemblorg/ensembl-vep
docker run -t -i -v $HOME/vep_data:/data ensemblorg/ensembl-vep INSTALL.pl -a c
sudo docker run -t -i -v $HOME/vep_data:/data ensemblorg/ensembl-vep INSTALL.pl -a cf -s homo_sapiens -y GRCh38
