FROM biohaz/basic_ubuntu:latest

# Locale
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# Anaconda
RUN wget -q https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O /tmp/Anaconda3-2020.02-Linux-x86_64.sh
RUN /bin/bash /tmp/Anaconda3-2020.02-Linux-x86_64.sh -bp /opt/Anaconda3
RUN rm /tmp/Anaconda3-2020.02-Linux-x86_64.sh

# Anaconda environment setting
RUN echo "export PATH=/opt/Anaconda3/bin:$PATH" > /etc/profile
ENV PATH /opt/Anaconda3/bin:$PATH

# Software
RUN apt-get update 
RUN pip3 install multiqc
RUN pip install configparser
RUN source activate
RUN conda update -n base -c defaults conda
RUN conda install -c bioconda -c conda-forge fastqc cutadapt bowtie bowtie2 star kraken2 bracken spades kaiju krona

ENTRYPOINT ["/usr/bin/python3", "./Script.py", "-c", "config_file.txt"]
#END
