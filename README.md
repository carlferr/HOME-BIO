# SHOME-BIO

SHOtgun MEtagenomic analysis of BIOlogical entities (SHOME-BIO) is a dockerized solution for metagenomics. Inside a docker container, we installed UBUNTU with python 3.7 and Anaconda 3 (V. 2020/02).
Inside an anaconda base environment, we installed two common metagenomics software, KAIJU (http://kaiju.binf.ku.dk/) and KRAKEN 2 (https://ccb.jhu.edu/software/kraken2/), and other mandatory softwares for the pipeline. For the complete list, please read our paper here (link)


## 1 - DOCKER INSTALLATION

To use our pipeline, it is mandatory to install and run DOCKER (https://hub.docker.com/). Here (https://hub.docker.com/search?q=&type=edition&offering=community) you can find the correct version of DOCKER for your OS and all the infos about how to install and run docker. 


## 2 - PIPELINE INSTALLATION

Run DOCKER on your PC or server and  install our Docker image pulling it from here (https://hub.docker.com/r/biohaz/metagenomic) or just type in your console:
```
docker pull biohaz/shome_bio:latest
```

Now, just download this repository or, in your console, type: 

```
git clone https://github.com/carlferr/SHOME-BIO.git
```


## 3 - BEGIN THE ANALYSIS

### 3.1 Prepare the Config file

Before run SHOME-BIO, manually change the "config_file.txt". In this way, it is possible to choose the correct options for your analysis.

### 3.2 - Download the databases

In order to run a metagenomic analysis, you should download different genome reference for bacteria, protozoa or virus.
We provide them and you can download them directly from here (https://drive.google.com/drive/folders/17PrBIJAjBP6XApBEvfBOsxfYliMsyVgf?usp=sharing). Inside the "Input" folder, it is also possible to download two .fastq test files in order to test the pipeline (this is optional).

### 3.3 - Run the pipeline

It is possible, now, to run SHOME-BIO just typing in the console:

```
docker run -it --rm -v /your_path/SHOME-BIO/Script.py:/home/Script.py:ro -v /your_path/SHOME-BIO/config_file.txt:/home/config.txt:ro -v /your_fastq_path/Fastq_folder:/home/Input:ro -v /your_output_folder:/home/Output:rw -v /your_hg19_path/Bowtie2Index:/home/Genome:ro -v /your_path/KRAKENdb_bacteria:/home/Db_Kraken2_Kaiju_bacteria:ro -v /your_path/KRAKENdb_protozoa:/home/Db_Kraken2_Kaiju_protozoa:ro -v /your_path/KRAKENdb_viruses:/home/Db_Kraken2_Kaiju_viruses:ro -v /your_path/KAIJUdb:/home/Db_Kaiju:ro -v /your_path/KAIJUdb_virus:/home/Db_Kaiju_virus:ro   biohaz/shome_bio
```
Please change "your_path" with your exact path. Each path after the -v option will be imported in the Docker container.
Running this command in your console, automatically it will call the Docker container and it will launch the analysis.

## 4 - LICENSE
This is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
