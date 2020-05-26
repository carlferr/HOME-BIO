#!/usr/bin/python

import sys, getopt
import os
import string
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

def main(argv):
	# Choose of inputfile e configfile	
	configfile = ''
	inputfile_a = ''
        inputfile_b = ''
	try:
		opts, args = getopt.getopt(argv,"hc:",["cfile="])
	except getopt.GetoptError:
		print 'test.py -c <configfile> '
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -c <configfile> '
			sys.exit()
		elif opt in ("-c", "--cfile"):
			configfile = arg
	print 'Config file is : ', configfile
	

	# Reading Config File
	file_config=open(configfile,"r")
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_fastqc=line_config.split()
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_fastqc2=line_config.split()
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_n_thread=line_config.split()
	n_thread=split_n_thread[0]
	flag=n_thread.isdigit()
	n_thread="1"
	if flag :
		n_thread=split_n_thread[0]
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_path_fastq=line_config.split()
	path_fastq=split_path_fastq[0]
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_pairend=line_config.split()
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_output_path=line_config.split()
	os.system("mkdir -p "+split_output_path[0])
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_path_genome=line_config.split()
	path_genome=split_path_genome[0]
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_delete=line_config.split()
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_path_delete=line_config.split()
	path_delete=split_path_delete[0]
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_adapter=line_config.split()
	adapter=split_adapter[0]
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_shotgun=line_config.split()
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	split_denovo=line_config.split()
	line_config=file_config.readline()
	line_config=file_config.readline()
	line_config=file_config.readline()
	
	file_config.close()	
	
	split_output_path[0]="/home/Output/"
	path_fastq="/home/Input/"
	path_genome="/home/Genome/genome"

	fastq=os.listdir(path_fastq)
	fastq.sort(lambda x, y: cmp(string.lower(x), string.lower(y)))
	for i in range(0,len(fastq),2):
		sample_only_name_a=fastq[i].split(".fastq")
		sample_only_name_b=fastq[i+1].split(".fastq")
		sample_only_name_a_noR1=sample_only_name_a[0].split("_R1")
		inputfile_a=path_fastq+"/"+fastq[i]
		inputfile_b=path_fastq+"/"+fastq[i+1]
			
	   	# FASTQC before adapter trimming 
		
		if split_fastqc[0] == "yes" :
			print("FASTQC before adapter trimming")
			os.system("mkdir -p "+split_output_path[0]+"/1_Output_FastQC_before_Trimming")
			os.system("fastqc "+inputfile_a+" "+inputfile_b+" -o "+split_output_path[0]+"/1_Output_FastQC_before_Trimming/"+" -t "+n_thread)
			#os.system("fastqc "+inputfile_a+" -o "+split_output_path[0]+"/1_Output_FastQC_before_Trimming/"+" -t "+n_thread)
					
		# Running MultiQC before adapter trimming

		if split_fastqc[0] == "yes" :
			print("MULTIQC before adapter trimming")	
			os.system("mkdir -p "+split_output_path[0]+"/1_Output_multiQC_before_Trimming")
			os.system("multiqc "+split_output_path[0]+"/1_Output_FastQC_before_Trimming/*.zip -o "+split_output_path[0]+"/1_Output_multiQC_before_Trimming/")	

		# Running trimmomatic for adapter trimming

		print("Running trimmomatic for adapter trimming")
		os.system("mkdir -p "+split_output_path[0]+"/2_Fastq_trimmed")
		#os.system("trimmomatic PE -threads "+n_thread+" -trimlog "+split_output_path[0]+"/2_Fastq_trimmed/trimm_log_file.txt -summary "+split_output_path[0]+"/2_Fastq_trimmed/trimm_summary_file.txt "+inputfile_a+" "+inputfile_b+" "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq  "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"unpaired_trimmed.fastq "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq  "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"unpaired_trimmed.fastq ILLUMINACLIP:/opt/Anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 CROP:76")
		#nextereCTGTCTCTTATA 
		os.system("cutadapt -m 20 -a "+adapter+" -A "+adapter+" -j "+n_thread+" -o "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -p "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq "+inputfile_a+" "+inputfile_b)
                                                                                           
		# FASTQC before adapter trimming 
		
		if split_fastqc2[0] == "yes" :
			print("FASTQC after adapter trimming")
			os.system("mkdir -p "+split_output_path[0]+"/1_Output_FastQC_after_Trimming")
			os.system("fastqc "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -o "+split_output_path[0]+"/1_Output_FastQC_after_Trimming/"+" -t "+n_thread)
					
		# Running MultiQC before adapter trimming

		if split_fastqc2[0] == "yes" :
			print("MULTIQC after adapter trimming")	
			os.system("mkdir -p "+split_output_path[0]+"/1_Output_multiQC_after_Trimming")
			os.system("multiqc "+split_output_path[0]+"/1_Output_FastQC_after_Trimming/*.zip -o "+split_output_path[0]+"/1_Output_multiQC_after_Trimming/")

		if split_delete[0] == "yes" :
			print("Running bowtie2 for alignment on reference genome and delete host organism")
			os.system("mkdir -p "+split_output_path[0]+"/3_Bowtie2")
			os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_metrics_file_host.txt --un-conc "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align_host.fastq --al-conc "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_align_host.fastq --un "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_unpaired_reads_not_align_host.fastq --al "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_unpaired_reads_align_host.fastq -x "+path_delete+" -1 "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -2 "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -S "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"host.sam")
			os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_metrics_file.txt --un-conc "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.fastq --al-conc "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_align.fastq --un "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_unpaired_reads_not_align.fastq --al "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_unpaired_reads_align.fastq -x "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align_host_1.fastq -2 "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align_host_2.fastq -S "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+".sam")

		else :
		
			# Running "bowtie2" for alignment on reference genome
			
			print("Running bowtie2 for alignment on reference genome")
			os.system("mkdir -p "+split_output_path[0]+"/3_Bowtie2")
			os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_metrics_file.txt --un-conc "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.fastq --al-conc "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_align.fastq --un "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_unpaired_reads_not_align.fastq --al "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+"_unpaired_reads_align.fastq -x "+path_genome+" -1 "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -2 "+split_output_path[0]+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -S "+split_output_path[0]+"/3_Bowtie2/"+sample_only_name_a_noR1[0]+".sam")

		if split_shotgun[0] == "yes" :
			# Running "kraken" for bacteria
			
			print("Running kraken for bacteria")
			os.system("mkdir -p "+split_output_path[0]+"/4_Kraken_bacteria")
			os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq  --use-names --report "+split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt --output "+split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_output_kraken.txt --confidence 0.5 --paired "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq")   

			# Running "kraken" for viruses
			
			print("Running kraken for viruses")
			os.system("mkdir -p "+split_output_path[0]+"/4_Kraken_viruses")
			os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt --output "+split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_output_kraken.txt --confidence 0.5 --paired "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq")

			# Running "kraken" for protozoa
			
			print("Running kraken for protozoa")
			os.system("mkdir -p "+split_output_path[0]+"/4_Kraken_protozoa")
			os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt --output "+split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_output_kraken.txt --confidence 0.5 --paired "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq") 

			# Running "kaiju"

			print("Running kaiju")
			os.system("mkdir -p "+split_output_path[0]+"/7_kaiju")
			os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq -j "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq -o "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out")
			os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -r species -o "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out") 
			os.system("kaiju2table -t Db_Kaiju_virus/nodes.dmp -n Db_Kaiju_virus/names.dmp -r species -o "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary2.tsv "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out")
			os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju//names.dmp -i "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out -o "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.krona")
			os.system("ktImportText -o "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.krona")

			#Results 4_kraken
			file_o=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred.txt","w")
			countG=0
			countS=0
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="G":
					output_o.write(line_o)
					countG=countG+1
				if split_o_2[0]=="S":
					output_o.write(line_o)
					countS=countS+1
				line_o=file_o.readline()
			output_o.write("Tot G="+str(countG)+" Tot S="+str(countS))	
			file_o.close()
			output_o.close()
				
			file_o=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","w")
			output_o.write("Tax_id"+"\t"+"Read_Count"+"\t"+"Percentage"+"\n")
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="S":
					lunghezza=len(split_o)
					x=5
					linea=""
					while x < lunghezza-1:
						linea=linea+split_o[x]+"_"
						x=x+1
					linea=linea+split_o[lunghezza-1]
					output_o.write(linea+"\t"+split_o[2]+"\t"+split_o[0]+"\n")
				line_o=file_o.readline()
			file_o.close()
			output_o.close()

			file_o=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred.txt","w")
			countG=0
			countS=0
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="G":
					output_o.write(line_o)
					countG=countG+1
				if split_o_2[0]=="S":
					output_o.write(line_o)
					countS=countS+1
				line_o=file_o.readline()
			output_o.write("Tot G="+str(countG)+" Tot S="+str(countS))	
			file_o.close()
			output_o.close()
				
			file_o=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","w")
			output_o.write("Tax_id"+"\t"+"Read_Count"+"\t"+"Percentage"+"\n")
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="S":
					lunghezza=len(split_o)
					x=5
					linea=""
					while x < lunghezza-1:
						linea=linea+split_o[x]+"_"
						x=x+1
					linea=linea+split_o[lunghezza-1]
					output_o.write(linea+"\t"+split_o[2]+"\t"+split_o[0]+"\n")
				line_o=file_o.readline()
			file_o.close()
			output_o.close()

			file_o=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred.txt","w")
			countG=0
			countS=0
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="G":
					output_o.write(line_o)
					countG=countG+1
				if split_o_2[0]=="S":
					output_o.write(line_o)
					countS=countS+1
				line_o=file_o.readline()
			output_o.write("Tot G="+str(countG)+" Tot S="+str(countS))	
			file_o.close()
			output_o.close()
				
			file_o=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","w")
			output_o.write("Tax_id"+"\t"+"Read_Count"+"\t"+"Percentage"+"\n")
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="S":
					lunghezza=len(split_o)
					x=5
					linea=""
					while x < lunghezza-1:
						linea=linea+split_o[x]+"_"
						x=x+1
					linea=linea+split_o[lunghezza-1]
					output_o.write(linea+"\t"+split_o[2]+"\t"+split_o[0]+"\n")
				line_o=file_o.readline()
			file_o.close()
			output_o.close()

			#CROSSValidation

			file_o1=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_crossvalidated.txt","w")

			line_o1=file_o1.readline()		
			split=line_o1.split("\n")
			output_o.write(split[0]+"\t CrossValidated\n")

			while(line_o1!=""):
				split1=line_o1.split()
				file_o2=open(split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
				line_o2=file_o2.readline()
				trovato=0
				while(line_o2!=""):			
					split2=line_o2.split("\t")
					if split2[3]==split1[4]:
						trovato=1
					line_o2=file_o2.readline()
				if trovato ==1:
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t CrossValidated\n")			
				line_o1=file_o1.readline()	
				file_o2.close()
			
			file_o1.close()
			output_o.close()

			file_o1=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_crossvalidated.txt","w")

			line_o1=file_o1.readline()		
			split=line_o1.split("\n")
			output_o.write(split[0]+"\t CrossValidated\n")

			while(line_o1!=""):
				split1=line_o1.split()
				file_o2=open(split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
				line_o2=file_o2.readline()
				trovato=0
				while(line_o2!=""):			
					split2=line_o2.split("\t")
					if split2[3]==split1[4]:
						trovato=1
					line_o2=file_o2.readline()
				if trovato ==1:
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t CrossValidated\n")			
				line_o1=file_o1.readline()	
				file_o2.close()
			
			file_o1.close()
			output_o.close()

			file_o1=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_crossvalidated.txt","w")

			line_o1=file_o1.readline()		
			split=line_o1.split("\n")
			output_o.write(split[0]+"\t CrossValidated\n")

			while(line_o1!=""):
				split1=line_o1.split()
				file_o2=open(split_output_path[0]+"7_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
				line_o2=file_o2.readline()
				trovato=0
				while(line_o2!=""):			
					split2=line_o2.split("\t")
					if split2[3]==split1[4]:
						trovato=1
					line_o2=file_o2.readline()
				if trovato ==1:
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t CrossValidated\n")			
				line_o1=file_o1.readline()	
				file_o2.close()
			
			file_o1.close()
			output_o.close()

			input=open(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","r")
			readline=input.readline()
			readline=input.readline()
			langs2=""
			students2=""
			while readline!="" :  
			    split=readline.split()
			    langs2=langs2+split[0]+"\t"
			    students2=students2+split[1]+"\t"
			    readline=input.readline()  

			langs3=langs2.split()
			students3=students2.split()

			#print(students3)
			n = len(students3)
			for i in range(n-1): 
			    for j in range(0, n-i-1):
				if int(students3[j]) < int(students3[j+1]) : 
				    students3[j], students3[j+1] = students3[j+1], students3[j]
				    langs3[j], langs3[j+1] = langs3[j+1], langs3[j]
			#print(students3)

			total=0
			labels = ["%s" % i for i in langs3[0:9]]
			for x in students3[0:9]:
			    total=total+int(x)
			fig1, ax1 = plt.subplots(figsize=(10, 10))
			fig1.subplots_adjust(0.3, 0, 1, 1)
			 
			theme = plt.get_cmap('Paired')
			ax1.set_prop_cycle("color", [theme(1. * i / len(students3[0:9]))
						     for i in range(len(students3[0:9]))])
			 
			_, _ = ax1.pie(students3[0:9], startangle=90, radius=1800)
			 
			ax1.axis('equal')
			plt.legend(
			    loc='upper left',
			    labels=['%s, %1.1f%%' % (
				l, (float(s) / total) * 100)
				    for l, s in zip(labels, students3[0:9])],
			    prop={'size': 11},
			    bbox_to_anchor=(0.0, 1),
			    bbox_transform=fig1.transFigure
			)
			
			plt.show()
			plt.savefig(split_output_path[0]+"/4_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"piechart.png")
			input.close()

			input=open(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","r")
			readline=input.readline()
			readline=input.readline()
			langs2=""
			students2=""
			while readline!="" :  
			    split=readline.split()
			    langs2=langs2+split[0]+"\t"
			    students2=students2+split[1]+"\t"
			    readline=input.readline()  

			langs3=langs2.split()
			students3=students2.split()

			#print(students3)
			n = len(students3)
			for i in range(n-1): 
			    for j in range(0, n-i-1):
				if int(students3[j]) < int(students3[j+1]) : 
				    students3[j], students3[j+1] = students3[j+1], students3[j]
				    langs3[j], langs3[j+1] = langs3[j+1], langs3[j]
			#print(students3)

			total=0
			labels = ["%s" % i for i in langs3[0:9]]
			for x in students3[0:9]:
			    total=total+int(x)
			fig1, ax1 = plt.subplots(figsize=(10, 10))
			fig1.subplots_adjust(0.3, 0, 1, 1)
			 
			theme = plt.get_cmap('Paired')
			ax1.set_prop_cycle("color", [theme(1. * i / len(students3[0:9]))
						     for i in range(len(students3[0:9]))])
			 
			_, _ = ax1.pie(students3[0:9], startangle=90, radius=1800)
			 
			ax1.axis('equal')
			plt.legend(
			    loc='upper left',
			    labels=['%s, %1.1f%%' % (
				l, (float(s) / total) * 100)
				    for l, s in zip(labels, students3[0:9])],
			    prop={'size': 11},
			    bbox_to_anchor=(0.0, 1),
			    bbox_transform=fig1.transFigure
			)
			
			plt.show()
			plt.savefig(split_output_path[0]+"/4_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"piechart.png")
			input.close()

			input=open(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","r")
			readline=input.readline()
			readline=input.readline()
			langs2=""
			students2=""
			while readline!="" :  
			    split=readline.split()
			    langs2=langs2+split[0]+"\t"
			    students2=students2+split[1]+"\t"
			    readline=input.readline()  

			langs3=langs2.split()
			students3=students2.split()

			#print(students3)
			n = len(students3)
			for i in range(n-1): 
			    for j in range(0, n-i-1):
				if int(students3[j]) < int(students3[j+1]) : 
				    students3[j], students3[j+1] = students3[j+1], students3[j]
				    langs3[j], langs3[j+1] = langs3[j+1], langs3[j]
			#print(students3)

			total=0
			labels = ["%s" % i for i in langs3[0:9]]
			for x in students3[0:9]:
			    total=total+int(x)
			fig1, ax1 = plt.subplots(figsize=(10, 10))
			fig1.subplots_adjust(0.3, 0, 1, 1)
			 
			theme = plt.get_cmap('Paired')
			ax1.set_prop_cycle("color", [theme(1. * i / len(students3[0:9]))
						     for i in range(len(students3[0:9]))])
			 
			_, _ = ax1.pie(students3[0:9], startangle=90, radius=1800)
			 
			ax1.axis('equal')
			plt.legend(
			    loc='upper left',
			    labels=['%s, %1.1f%%' % (
				l, (float(s) / total) * 100)
				    for l, s in zip(labels, students3[0:9])],
			    prop={'size': 11},
			    bbox_to_anchor=(0.0, 1),
			    bbox_transform=fig1.transFigure
			)
			
			plt.show()
			plt.savefig(split_output_path[0]+"/4_Kraken_viruses/"+sample_only_name_a_noR1[0]+"piechart.png")
			input.close()


		if split_denovo[0] == "yes" :		
			# Running "SPADES" 
			
			print("Running SPADES")
			os.system("mkdir -p "+split_output_path[0]+"/5_SPADES")
			os.system("spades.py --meta -t "+n_thread+" -1 "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq -2 "+split_output_path[0]+"3_Bowtie2/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq -o "+split_output_path[0]+"/5_SPADES/"+sample_only_name_a_noR1[0])

			# Running "kraken" for bacteria
			
			print("Running kraken for bacteria")
			os.system("mkdir -p "+split_output_path[0]+"/6_Kraken_bacteria")
			os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt --output "+split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_output_kraken.txt --confidence 0.5 "+split_output_path[0]+"/5_SPADES/"+sample_only_name_a_noR1[0]+"/contigs.fasta")
			

			# Running "kraken" for viruses

			print("Running kraken for viruses")
			os.system("mkdir -p "+split_output_path[0]+"/6_Kraken_viruses")
			os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt --output "+split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_output_kraken.txt --confidence 0.5 "+split_output_path[0]+"/5_SPADES/"+sample_only_name_a_noR1[0]+"contigs.fasta")

			# Running "kraken" for protozoa

			print("Running kraken for protozoa")
			os.system("mkdir -p "+split_output_path[0]+"/6_Kraken_protozoa")
			os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt --output "+split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_output_kraken.txt --confidence 0.5 "+split_output_path[0]+"/5_SPADES/"+sample_only_name_a_noR1[0]+"contigs.fasta")		

			os.system("mkdir -p "+split_output_path[0]+"/8_kaiju")
			os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path[0]+"5_SPADES/"+sample_only_name_a_noR1[0]+"/contigs.fasta -o "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out")
			os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -r species -o "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out") 
			os.system("kaiju2table -t Db_Kaiju_virus/nodes.dmp -n Db_Kaiju_virus/names.dmp -r species -o "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary2.tsv "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out")
			os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.out -o "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.krona")
			os.system("ktImportText -o "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_e_val-3.krona")

			#Results 6_kraken
			file_o=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred.txt","w")
			countG=0
			countS=0
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="G":
					output_o.write(line_o)
					countG=countG+1
				if split_o_2[0]=="S":
					output_o.write(line_o)
					countS=countS+1
				line_o=file_o.readline()
			output_o.write("Tot G="+str(countG)+" Tot S="+str(countS))	
			file_o.close()
			output_o.close()
				
			file_o=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","w")
			output_o.write("Tax_id"+"\t"+"Read_Count"+"\t"+"Percentage"+"\n")
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="S":
					lunghezza=len(split_o)
					x=5
					linea=""
					while x < lunghezza-1:
						linea=linea+split_o[x]+"_"
						x=x+1
					linea=linea+split_o[lunghezza-1]
					output_o.write(linea+"\t"+split_o[2]+"\t"+split_o[0]+"\n")
				line_o=file_o.readline()
			file_o.close()
			output_o.close()

			file_o=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred.txt","w")
			countG=0
			countS=0
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="G":
					output_o.write(line_o)
					countG=countG+1
				if split_o_2[0]=="S":
					output_o.write(line_o)
					countS=countS+1
				line_o=file_o.readline()
			output_o.write("Tot G="+str(countG)+" Tot S="+str(countS))	
			file_o.close()
			output_o.close()
				
			file_o=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","w")
			output_o.write("Tax_id"+"\t"+"Read_Count"+"\t"+"Percentage"+"\n")
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="S":
					lunghezza=len(split_o)
					x=5
					linea=""
					while x < lunghezza-1:
						linea=linea+split_o[x]+"_"
						x=x+1
					linea=linea+split_o[lunghezza-1]
					output_o.write(linea+"\t"+split_o[2]+"\t"+split_o[0]+"\n")
				line_o=file_o.readline()
			file_o.close()
			output_o.close()

			file_o=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred.txt","w")
			countG=0
			countS=0
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="G":
					output_o.write(line_o)
					countG=countG+1
				if split_o_2[0]=="S":
					output_o.write(line_o)
					countS=countS+1
				line_o=file_o.readline()
			output_o.write("Tot G="+str(countG)+" Tot S="+str(countS))	
			file_o.close()
			output_o.close()
				
			file_o=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","w")
			output_o.write("Tax_id"+"\t"+"Read_Count"+"\t"+"Percentage"+"\n")
			line_o=file_o.readline()
			while(line_o!=""):
				split_o=line_o.split()
				split_o_2=split_o[3]
				if split_o_2[0]=="S":
					lunghezza=len(split_o)
					x=5
					linea=""
					while x < lunghezza-1:
						linea=linea+split_o[x]+"_"
						x=x+1
					linea=linea+split_o[lunghezza-1]
					output_o.write(linea+"\t"+split_o[2]+"\t"+split_o[0]+"\n")
				line_o=file_o.readline()
			file_o.close()
			output_o.close()


			#CROSSVALIDAZIONE	

			file_o1=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_crossvalidated.txt","w")

			line_o1=file_o1.readline()		
			split=line_o1.split("\n")
			output_o.write(split[0]+"\t CrossValidated\n")

			while(line_o1!=""):
				split1=line_o1.split()
				file_o2=open(split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
				line_o2=file_o2.readline()
				trovato=0
				while(line_o2!=""):			
					split2=line_o2.split("\t")
					if split2[3]==split1[4]:
						trovato=1
					line_o2=file_o2.readline()
				if trovato ==1:
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t CrossValidated\n")			
				line_o1=file_o1.readline()	
				file_o2.close()
			
			file_o1.close()
			output_o.close()

			file_o1=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_crossvalidated.txt","w")

			line_o1=file_o1.readline()		
			split=line_o1.split("\n")
			output_o.write(split[0]+"\t CrossValidated\n")

			while(line_o1!=""):
				split1=line_o1.split()
				file_o2=open(split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
				line_o2=file_o2.readline()
				trovato=0
				while(line_o2!=""):			
					split2=line_o2.split("\t")
					if split2[3]==split1[4]:
						trovato=1
					line_o2=file_o2.readline()
				if trovato ==1:
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t CrossValidated\n")			
				line_o1=file_o1.readline()	
				file_o2.close()
			
			file_o1.close()
			output_o.close()

			file_o1=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken.txt","r")
			output_o=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_crossvalidated.txt","w")

			line_o1=file_o1.readline()		
			split=line_o1.split("\n")
			output_o.write(split[0]+"\t CrossValidated\n")

			while(line_o1!=""):
				split1=line_o1.split()
				file_o2=open(split_output_path[0]+"8_kaiju/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
				line_o2=file_o2.readline()
				trovato=0
				while(line_o2!=""):			
					split2=line_o2.split("\t")
					if split2[3]==split1[4]:
						trovato=1
					line_o2=file_o2.readline()
				if trovato ==1:
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t CrossValidated\n")			
				line_o1=file_o1.readline()	
				file_o2.close()
			
			file_o1.close()
			output_o.close()

			#piechart

			
			input=open(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","r")
			readline=input.readline()
			readline=input.readline()
			langs2=""
			students2=""
			while readline!="" :  
			    split=readline.split()
			    langs2=langs2+split[0]+"\t"
			    students2=students2+split[1]+"\t"
			    readline=input.readline()  

			langs3=langs2.split()
			students3=students2.split()

			#print(students3)
			n = len(students3)
			for i in range(n-1): 
			    for j in range(0, n-i-1):
				if int(students3[j]) < int(students3[j+1]) : 
				    students3[j], students3[j+1] = students3[j+1], students3[j]
				    langs3[j], langs3[j+1] = langs3[j+1], langs3[j]
			#print(students3)

			total=0
			labels = ["%s" % i for i in langs3[0:9]]
			for x in students3[0:9]:
			    total=total+int(x)
			fig1, ax1 = plt.subplots(figsize=(10, 10))
			fig1.subplots_adjust(0.3, 0, 1, 1)
			 
			theme = plt.get_cmap('Paired')
			ax1.set_prop_cycle("color", [theme(1. * i / len(students3[0:9]))
						     for i in range(len(students3[0:9]))])
			 
			_, _ = ax1.pie(students3[0:9], startangle=90, radius=1800)
			 
			ax1.axis('equal')
			plt.legend(
			    loc='upper left',
			    labels=['%s, %1.1f%%' % (
				l, (float(s) / total) * 100)
				    for l, s in zip(labels, students3[0:9])],
			    prop={'size': 11},
			    bbox_to_anchor=(0.0, 1),
			    bbox_transform=fig1.transFigure
			)
			
			plt.show()
			plt.savefig(split_output_path[0]+"/6_Kraken_protozoa/"+sample_only_name_a_noR1[0]+"piechart.png")
			input.close()

			input=open(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","r")
			readline=input.readline()
			readline=input.readline()
			langs2=""
			students2=""
			while readline!="" :  
			    split=readline.split()
			    langs2=langs2+split[0]+"\t"
			    students2=students2+split[1]+"\t"
			    readline=input.readline()  

			langs3=langs2.split()
			students3=students2.split()

			#print(students3)
			n = len(students3)
			for i in range(n-1): 
			    for j in range(0, n-i-1):
				if int(students3[j]) < int(students3[j+1]) : 
				    students3[j], students3[j+1] = students3[j+1], students3[j]
				    langs3[j], langs3[j+1] = langs3[j+1], langs3[j]
			#print(students3)

			total=0
			labels = ["%s" % i for i in langs3[0:9]]
			for x in students3[0:9]:
			    total=total+int(x)
			fig1, ax1 = plt.subplots(figsize=(10, 10))
			fig1.subplots_adjust(0.3, 0, 1, 1)
			 
			theme = plt.get_cmap('Paired')
			ax1.set_prop_cycle("color", [theme(1. * i / len(students3[0:9]))
						     for i in range(len(students3[0:9]))])
			 
			_, _ = ax1.pie(students3[0:9], startangle=90, radius=1800)
			 
			ax1.axis('equal')
			plt.legend(
			    loc='upper left',
			    labels=['%s, %1.1f%%' % (
				l, (float(s) / total) * 100)
				    for l, s in zip(labels, students3[0:9])],
			    prop={'size': 11},
			    bbox_to_anchor=(0.0, 1),
			    bbox_transform=fig1.transFigure
			)
			
			plt.show()
			plt.savefig(split_output_path[0]+"/6_Kraken_bacteria/"+sample_only_name_a_noR1[0]+"piechart.png")
			input.close()

			input=open(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"_report_kraken_filtred_count.txt","r")
			readline=input.readline()
			readline=input.readline()
			langs2=""
			students2=""
			while readline!="" :  
			    split=readline.split()
			    langs2=langs2+split[0]+"\t"
			    students2=students2+split[1]+"\t"
			    readline=input.readline()  

			langs3=langs2.split()
			students3=students2.split()

			#print(students3)
			n = len(students3)
			for i in range(n-1): 
			    for j in range(0, n-i-1):
				if int(students3[j]) < int(students3[j+1]) : 
				    students3[j], students3[j+1] = students3[j+1], students3[j]
				    langs3[j], langs3[j+1] = langs3[j+1], langs3[j]
			#print(students3)

			total=0
			labels = ["%s" % i for i in langs3[0:9]]
			for x in students3[0:9]:
			    total=total+int(x)
			fig1, ax1 = plt.subplots(figsize=(10, 10))
			fig1.subplots_adjust(0.3, 0, 1, 1)
			 
			theme = plt.get_cmap('Paired')
			ax1.set_prop_cycle("color", [theme(1. * i / len(students3[0:9]))
						     for i in range(len(students3[0:9]))])
			 
			_, _ = ax1.pie(students3[0:9], startangle=90, radius=1800)
			 
			ax1.axis('equal')
			plt.legend(
			    loc='upper left',
			    labels=['%s, %1.1f%%' % (
				l, (float(s) / total) * 100)
				    for l, s in zip(labels, students3[0:9])],
			    prop={'size': 11},
			    bbox_to_anchor=(0.0, 1),
			    bbox_transform=fig1.transFigure
			)
			
			plt.show()
			plt.savefig(split_output_path[0]+"/6_Kraken_viruses/"+sample_only_name_a_noR1[0]+"piechart.png")
			input.close()
			
if __name__ == "__main__":
	main(sys.argv[1:])


