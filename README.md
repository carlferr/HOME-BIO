# HOME-BIO

Despite the large amount of metagenomic shotgun data produced, there is a lack of a comprehensive and easy-use pipeline for data analysis that avoid annoying and complicated bioinformatics steps. Here we present **HOME-BIO** (sHOtgun MEtagenomic analysis of BIOlogical entities), a modular and exhaustive pipeline for analysis of biological entity estimation, specific designed for shotgun sequenced clinical samples. HOME-BIO analysis provides comprehensive taxonomy classification by querying different source database and carry out main steps in metagenomic investigation.

HOME-BIO is a dockerized solution for metagenomics. Inside a DOCKER image, we installed UBUNTU with python 3.7, python 2.7, and Anaconda 3 (V. 2020/02). Inside an anaconda base environment, we installed two common metagenomics software, [KAIJU](http://kaiju.binf.ku.dk/) and [KRAKEN 2](https://ccb.jhu.edu/software/kraken2/), and other mandatory softwares for the pipeline. For the complete list, please read our paper here (link)


## 1 - DOCKER & PYTHON INSTALLATION

To use our pipeline, it is mandatory to install and run DOCKER (https://hub.docker.com/), and to install and use python 2.7.

- [Here](https://hub.docker.com/search?q=&type=edition&offering=community) you can find the correct version of DOCKER for your OS and all the infos about how to install and run DOCKER.

- To install python 2.7 there are several options depending on your OS:

#### Linux\Debian
```
sudo apt install python2.7
```

#### Centos\RHEL\Fedora
```
sudo yum install python27
```

#### Windows\MacOS

Just download and install it from [here](https://www.python.org/downloads/release/python-2718/)


## 2 - PIPELINE INSTALLATION

Run DOCKER on your PC or server and  install our DOCKER image pulling it from [here](https://hub.docker.com/r/biohaz/home_bio) or just type in your console:
```
docker pull biohaz/home_bio:latest
```

Now, just download this repository or, in your console, type: 

```
git clone https://github.com/carlferr/HOME-BIO.git
```


## 3 - BEGIN THE ANALYSIS

### 3.1 Prepare the Config file

Before run HOME-BIO, manually change the "config_file.txt". In this way, it is possible to choose the correct options for your analysis. The user can modify this file following the written examples. For some options is mandatory a complete path, while others just require a "yes" or "no". Only the "Adapter" option requires a sequence to use as adapter during the trimming step.
If the contaminat filtering is not requested the "Path contaminant genome" is not used.
Please change "your_path" with your absolute path.

### 3.2 - Download the databases

In order to run a metagenomic analysis, you should download different genome reference for bacteria, protozoa or virus.
We provide a previously indexed version of all of them [here](https://drive.google.com/drive/folders/17PrBIJAjBP6XApBEvfBOsxfYliMsyVgf?usp=sharing) (Db_Kraken2_Kaiju.zip). It is possible to download the entire zipped database folder, but, before run the analysis, please unzip it on your machine.
Moreover, inside the "Input" folder, it is also possible to download two paried-end .fastq.gz test files in order to test the pipeline on your system.

### 3.3 - Run the pipeline

HOME-BIO uses fastq (or fastq.gz) files as input. Remember to put all your data in one unique folder (e.g. "Fastq_folder", "Data", "Input").
It is possible, now, to run HOME-BIO just typing in the console:

```
python2.7 HOME_Bio.py -c config_file.txt
```
Running this command in your console, it automatically will call the Docker container, read the path from config_file.txt and it will launch the analysis.

## 4 - OUTPUT

HOME-BIO modules (Shotgun metagenomic and Assembly de novo) will produce output files in tab and graphic format (.png). In your output folder, it is possible to find, for each sample, some folders. 

Most of them are output numbered folders for each tool used (e.g. "3_Bowtie2" with all the files generated from the alignment). 

The "Results" forlder has all the output files generated from the analysis in tabular and .png format.
In particular  the ‘Metagenomic shotgun’ module generates a table containing the Kranken2 taxonomy profile and related Kaiju protein-validation information. A given taxon is considered protein-validated when both tools classify and assign reads to it. In addition, HOME-BIO generates output pie-charts in .png format with top 15 represented species with estimation of the relative abundance.
