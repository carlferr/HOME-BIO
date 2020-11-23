# HOME-BIO

Despite the large amount of metagenomic shotgun data produced, there is a lack of a comprehensive and easy-use pipeline for data analysis that avoids annoying and complicated bioinformatics steps. Here, we present **HOME-BIO** (sHOtgun MEtagenomic analysis of BIOlogical entities), a modular and exhaustive pipeline for analysis of biological entity estimation, specifically designed for shotgun sequenced clinical samples. HOME-BIO analysis provides a comprehensive taxonomy classification by querying different source database and carry out the main steps in the metagenomic investigation.

HOME-BIO is a dockerized solution for metagenomics. Inside a DOCKER image, we installed UBUNTU with Python 3.8 and Anaconda 3 (V. 2020/02). Inside an anaconda base environment, we installed two common metagenomics software, [KAIJU](http://kaiju.binf.ku.dk/) and [KRAKEN 2](https://ccb.jhu.edu/software/kraken2/), and other mandatory software for the pipeline. For more details, please read our paper here (link)


## 1 - DOCKER & PYTHON INSTALLATION

To use our pipeline, it is mandatory to install and run DOCKER (https://hub.docker.com/) and to install and use Python 3.

- [Here](https://hub.docker.com/search?q=&type=edition&offering=community) you can find the correct version of DOCKER for your OS and all the info about how to install and run DOCKER.

- To install Python 3 there are several options depending on your OS:

#### Linux\Debian
```
sudo apt update
sudo apt-get install python3.8
```

#### Centos\RHEL\Fedora
```
sudo yum update -y
sudo yum install -y python3
```

#### Windows\MacOS

Just download and install it from [here](https://www.python.org/downloads/windows/)


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

Before run HOME-BIO, the users can manually change the "config_file.txt" or run the "make_config_file.py". 

- In the first scenario, the user will manually modify the correct options for the analysis following the written examples. For some options is mandatory an absolute path, while for others is required a "yes" or "no". Only the "Adapter" option requires a sequence (in capital letters) to use as adapter during the trimming step.
For the "k-mers" option, HOME-BIO is set to "auto". In this way, SPADES will choose the best k-mer length depending on the reads length. Only a comma-separated list of k-mer sizes can be used (all values must be odd, less than 128 and listed in ascending order. e.g., 21,33,55). For more details see [SPADES manual](http://cab.spbu.ru/files/release3.13.0/manual.html).
If the contaminant filtering is not requested the "Path contaminant genome" is not used but it is still mandatory to write something (e.g., "/your_path/" or "none").


- It is also possible to create the config file in an automatic way using "make_config_file.py".
Just type
```
python2.7 make_config_file.py
```
The script will prompt the questions on the screen. Please, use complete pathways and answer all the questions.
It will create a "config.txt" file.

**If the automatic way is used, remember to run the pipeline using "-c config.txt"**

### 3.2 - Download the databases

To run a metagenomic analysis, you should download genome/protein reference databases for bacteria, protozoa, or virus.
We provide a previously indexed version of all of them [on Zenodo](https://doi.org/10.5281/zenodo.4055180) (DB_KKRAKEN2_KAIJU.zip). It is possible to download the zipped folder, but, before running the analysis, please unzip it on your machine.

### 3.3 - Run the pipeline

HOME-BIO uses fastq (or fastq.gz) files as input. Remember to put all your data in one unique folder (e.g., "Fastq_folder", "Data", "Input").
It is possible, now, to run HOME-BIO just typing in the console (for Linux\Debian\Centos\RHEL\Fedora):

```
python2.7 HOME_Bio.py -c config_file.txt
```
or (if the config was generated automatically)
```
python2.7 HOME_Bio.py -c config.txt
```

For Windows users, type
```
C:\python27\python.exe HOME_Bio.py -c config_file.txt
```

Running this command in your console, it will automatically call the Docker container, read the paths from the config_file.txt, and it will launch the analysis.

It is possible to test the pipeline using our test dataset (just two paired-end .fastq files) freely available [here](https://doi.org/10.5281/zenodo.4061297) on Zenodo.


## 4 - OUTPUT

HOME-BIO modules (Shotgun metagenomic and Assembly de novo) will produce output files in tab-delimited (.txt) and graphic (.png) format. In the output folder, it is possible to find, numbered folders for each tool used (e.g. “3_Bowtie2_Alignment” with all the files generated from the alignment; “5_SPADES_Assembly_DeNovo” folder contains per-sample Assembly de novo module output). HOME-BIO generates also a log.txt file containing all the messages prompted by the script during the analysis.
The **"9_Results"** folder has the final output of each module performed in HOME-BIO with the results of the taxonomic classification for the analyzed samples.

## 5 - REQUIREMENTS

It is mandatory to have at least 400 GB of available space on the disk to extract the compressed folder of the provided databases. It is NOT mandatory to use our previously indexed databases, but the user will have to create his own databases and, then, index them with Kraken2 and Kaiju. See the instruction on Kraken2 and Kaiju user manuals.
Other resource requirements will vary from the size of input data being processed.

Minimum requirements: 1 CPU, 4 GB RAM.

Recommended requirements: 4+ CPU, 16+ GB RAM.
