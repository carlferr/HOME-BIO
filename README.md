# HOME-BIO

Despite the large amount of metagenomic shotgun data produced, there is a lack of a comprehensive and easy-use pipeline for data analysis that avoid annoying and complicated bioinformatics steps. Here we present **HOME-BIO** (sHOtgun MEtagenomic analysis of BIOlogical entities), a modular and exhaustive pipeline for analysis of biological entity estimation, specific designed for shotgun sequenced clinical samples. SHOME-BIO analysis provides comprehensive taxonomy classification by querying different source database and carry out main steps in metagenomic investigation.

HOME-BIO is a dockerized solution for metagenomics. Inside a DOCKER image, we installed UBUNTU with python 3.7 and Anaconda 3 (V. 2020/02). Inside an anaconda base environment, we installed two common metagenomics software, [KAIJU](http://kaiju.binf.ku.dk/) and [KRAKEN 2](https://ccb.jhu.edu/software/kraken2/), and other mandatory softwares for the pipeline. For the complete list, please read our paper here (link)


## 1 - DOCKER INSTALLATION

To use our pipeline, it is mandatory to install and run DOCKER (https://hub.docker.com/). [Here](https://hub.docker.com/search?q=&type=edition&offering=community) you can find the correct version of DOCKER for your OS and all the infos about how to install and run DOCKER. 


## 2 - PIPELINE INSTALLATION

Run DOCKER on your PC or server and  install our DOCKER image pulling it from [here](https://hub.docker.com/r/biohaz/shome_bio) or just type in your console:
```
docker pull biohaz/shome_bio:latest
```

Now, just download this repository or, in your console, type: 

```
git clone https://github.com/carlferr/SHOME-BIO.git
```


## 3 - BEGIN THE ANALYSIS

### 3.1 Prepare the Config file

Before run HOME-BIO, manually change the "config_file.txt". In this way, it is possible to choose the correct options for your analysis. The user can modify this file following the written examples. For some options is mandatory a complete path, while others just require a "yes" or "no". Only the "Adapter" option requires a sequence to use as adapter during the trimming step.
If the contaminat filtering is not requested the "Path contaminant genome" is not used.
Please change "your_path" with your absolute path.

### 3.2 - Download the databases

In order to run a metagenomic analysis, you should download different genome reference for bacteria, protozoa or virus.
We provide a previously indexed version of all of them [here](https://drive.google.com/drive/folders/17PrBIJAjBP6XApBEvfBOsxfYliMsyVgf?usp=sharing). It is possible to download the entire zipped database folder, but, before run the analysis, please unzip it on your machine.
Moreover, inside the "Input" folder, it is also possible to download two paried-end .fastq test files in order to test the pipeline on your system.

### 3.3 - Run the pipeline

HOME-BIO uses fastq (or fastq.gz) files as input. Remenber to put all your data in one unique folder ("Fastq_folder" in the code below)
It is possible, now, to run SHOME-BIO just typing in the console:

```
docker run -it --rm -v /your_path/HOME-BIO/Script.py:/home/Script.py:ro -v /your_path/HOME-BIO/config_file.txt:/home/config_file.txt:ro -v /your_fastq_path/Fastq_folder:/home/Input:ro -v /your_output_folder:/home/Output:rw -v /your_hg19_path/Bowtie2Index:/home/Genome:ro -v /your_path/KRAKENdb_bacteria:/home/Db_Kraken2_Kaiju_bacteria:ro -v /your_path/KRAKENdb_protozoa:/home/Db_Kraken2_Kaiju_protozoa:ro -v /your_path/KRAKENdb_viruses:/home/Db_Kraken2_Kaiju_viruses:ro -v /your_path/KAIJUdb:/home/Db_Kaiju:ro -v /your_path/KAIJUdb_virus:/home/Db_Kaiju_virus:ro   biohaz/shome_bio
```
**Please change "your_path" with your exact path.** Each path after the -v option will be imported in the Docker container.
Running this command in your console, it automatically will call the Docker container and it will launch the analysis.

## 4 - OUTPUT

HOME-BIO modules (Shotgun metagenomic and Assembly de novo) will produce output files in tab and graphic format (.png). In the Output folder, is it possible to find, for each sample, two folders. 

The first one with the output numbered folders for each tool used (e.g. "3_Bowtie2" with all the files generated from the alignment). 

The second one with output files generated from the analysis in tab and .png format.
In particular:
- .png files show abundance estimation of top 15 species classified in the queried databases (e.g. bacteria, virus, protozoi)
- report.txt file show classification data. All the entities classified are reported with information about protein validation.
- output.txt file show count of genera and species in each sample.
