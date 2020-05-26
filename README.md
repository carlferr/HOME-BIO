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

Now, just download this repository or 

```
docker pull biohaz/shome_bio:latest
```


## 3 - BEGIN THE ANALYSIS

### 3.1 Prepare the Config file

To run SHOME-BIO, just download this repository and manually change the "config.txt" file in it. In this way, it is possible to choose the correct options for your analysis.

### 3.2 - Download the databases

In order to run a metagenomic analysis, you should download different genome reference for bacteria, protozoa or virus.
We provide them and you can download them directly from here ().
We also provide a test data folder with two .fastq files and you can download them here () (this is optional).

### 3.3 - Run the Script.py

It is possible now to run SHOME-BIO just typing in the console

```
python ./Script.py
```

The script automatically will call the docker container and it will launch the analysis.

## 4 - LICENSE
This is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
