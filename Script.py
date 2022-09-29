#!/usr/bin/python

import sys, getopt
import os
import os.path
import string
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import configparser


def main(argv):
	# Choose of inputfile e configfile	
	configfile = ''
	inputfile_a = ''
	inputfile_b = ''
        
	try:
		opts, args = getopt.getopt(argv,"hc:",["cfile="])
	except getopt.GetoptError:
		print ('test.py -c <configfile> ')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('test.py -c <configfile> ')
			sys.exit()
		elif opt in ("-c", "--cfile"):
			configfile = arg
	print ('Config file is : ', configfile)
	

	config = configparser.ConfigParser()

	config.read(configfile)

	# Reading Config File
	
	split_fastqc=config['CUSTOM']['quality control before trimming [yes/no]']
	if split_fastqc!='yes':
		split_fastqc=config['DEFAULT']['quality control before trimming [yes/no]']
	print("quality control before trimming: "+split_fastqc)
	
	split_fastqc2=config['CUSTOM']['quality control after trimming [yes/no]']
	if split_fastqc2!='yes':
		split_fastqc2=config['DEFAULT']['quality control after trimming [yes/no]']
	print("quality control after trimming: "+split_fastqc2)

	n_thread=config['CUSTOM']['number of threads [n]']
	if n_thread=='':
		n_thread=config['DEFAULT']['number of threads [n]']
	print("number of threads: "+n_thread)

	path_fastq=config['CUSTOM']['path fastq']
	print("path fastq: "+path_fastq)
	
	split_pairend=config['CUSTOM']['paired-end? [yes/no]']
	if split_pairend!='no':
		split_pairend=config['DEFAULT']['paired-end? [yes/no]']
	print("paired-end: "+split_pairend)
	
	split_output_path=config['CUSTOM']['path output']
	print("path output: "+split_output_path)

	split_path_genome=config['CUSTOM']['path host genome']
	path_genome=split_path_genome
	print("path host genome: "+path_genome)

	split_delete=config['CUSTOM']['filter out contaminant [yes/no]']
	if split_delete!='yes':
		split_delete=config['DEFAULT']['filter out contaminant [yes/no]']
	print("filter out contaminant: "+split_delete)

	split_path_delete=config['CUSTOM']['path contaminant genome']
	path_delete=split_path_delete
	print("path contaminant genome: "+path_delete)

	split_adapter=config['CUSTOM']['adapter?']
	if split_adapter=='':
		split_adapter=config['DEFAULT']['adapter?']
	adapter=split_adapter
	print("adapter: "+adapter)

	split_shotgun=config['CUSTOM']['shotgun module [yes/no]']
	if split_shotgun!='yes':
		split_shotgun=config['DEFAULT']['shotgun module [yes/no]']
	print("shotgun module: "+split_shotgun)

	split_denovo=config['CUSTOM']['assembly denovo module [yes/no]']
	if split_denovo!='yes':
		split_denovo=config['DEFAULT']['assembly denovo module [yes/no]']
	print("assembly denovo module: "+split_denovo)

	split_path_kaiju=config['CUSTOM']['path kraken2 & kaiju databases']
	path_kaiju=split_path_kaiju
	print("path kraken2 & kaiju databases: "+path_kaiju)

	
	split_confidence=config['CUSTOM']['kraken2 confidence']
	if split_confidence=='':
		split_confidence=config['DEFAULT']['kraken2 confidence']
	print("kraken2 confidence: "+split_confidence)
	
	split_DNARNA=config['CUSTOM']['nucleic acid type [dna/rna]']
	if split_DNARNA!='RNA':
		split_DNARNA=config['DEFAULT']['nucleic acid type [dna/rna]']
	print("nucleic acid type: "+split_DNARNA)
	
	split_kmer=config['CUSTOM']['k-mer [auto/21,33,55,77]']
	if split_kmer=='':
		split_kmer=config['DEFAULT']['k-mer [auto/21,33,55,77]']
	print("k-mer: "+split_kmer)
	
		

		
	split_output_path="/home/Output/"
	path_fastq="/home/Input/"
	path_genome="/home/Genome/"
	base_name=[fname.rsplit('.', 2)[0] for fname in os.listdir(path_genome)][0]
	#os.seteuid(65534)
	output_log=open(split_output_path+"/log.txt","w")
	flag_ok="1"

	if split_DNARNA == "RNA":
		print("RNA")
		output_log.write("RNA\n")
		#pairend
		if split_pairend == "no":
			print("Single-end mode")
			output_log.write("Single-end mode\n")
			fastq=os.listdir(path_fastq)
			fastq.sort()
			#fastq.sort(lambda x, y: cmp(string.lower(x), string.lower(y)))
			for i in range(0,len(fastq),1):
				sample_only_name_a=fastq[i].split(".fastq")
				#sample_only_name_b=fastq[i+1].split(".fastq")
				sample_only_name_a_noR1=sample_only_name_a[0].split("_R1")
				inputfile_a=path_fastq+"/"+fastq[i]
				#inputfile_b=path_fastq+"/"+fastq[i+1]
					
			   	# FASTQC before adapter trimming 
				
				print("Sample in running: "+sample_only_name_a_noR1[0])
				output_log.write("Sample in running: "+sample_only_name_a_noR1[0]+"\n")
				if split_fastqc == "yes" :
					print("FASTQC before adapter trimming")
					output_log.write("FASTQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control")
					os.system("fastqc "+inputfile_a+" -o "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/"+" -t "+n_thread)
												
				# Running MultiQC before adapter trimming

				if split_fastqc == "yes" :
					print("MULTIQC before adapter trimming")
					output_log.write("MULTIQC before adapter trimming\n")	
					os.system("mkdir -p "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control/")	

				# Running trimmomatic for adapter trimming

				print("Running cutadapt for adapter trimming")
				output_log.write("Running cutadapt for adapter trimming\n")
				os.system("mkdir -p "+split_output_path+"/2_Fastq_trimmed")
				os.system("mkdir -p "+split_output_path+"/9_Results")

				#os.system("trimmomatic PE -threads "+n_thread+" -trimlog "+split_output_path+"/2_Fastq_trimmed/trimm_log_file.txt -summary "+split_output_path+"/2_Fastq_trimmed/trimm_summary_file.txt "+inputfile_a+" "+inputfile_b+" "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"unpaired_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"unpaired_trimmed.fastq ILLUMINACLIP:/opt/Anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 CROP:76")
				#nextera CTGTCTCTTATA 
				os.system("cutadapt -m 20 -a "+adapter+" -j "+n_thread+" -o "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+inputfile_a)
			
				path=os.path.isfile(split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+inputfile_a)
				if path:
					print("Cutadapt finished without error\n")
					output_log.write("Running cutadapt for adapter trimming\n")
				else:
					print("Error: cannot create adapter trimmed file\n")
					output_log.write("Running cutadapt for adapter trimming\n")
					flag_ok="0"
				                                                                           
				# FASTQC before adapter trimming 
				
				if split_fastqc2 == "yes" :
					print("FASTQC after adapter trimming")
					output_log.write("FASTQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control")
					os.system("fastqc "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -o "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/"+" -t "+n_thread)
							
				# Running MultiQC before adapter trimming

				if split_fastqc2 == "yes" :
					print("MULTIQC after adapter trimming")	
					output_log.write("MULTIQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control/")

				if split_delete == "yes" :
					print("Running STAR for alignment on reference genome and delete host organism")
					output_log.write("Running STAR for alignment on reference genome and delete host organism\n")
					os.system("mkdir -p "+split_output_path+"/3_STAR_Alignment")
					path_genome="/home/Genome/"
					os.system("STAR --runThreadN "+n_thread+" --genomeDir "+path_genome+" --readFilesIn "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq --outFileNamePrefix "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+" --outReadsUnmapped Fastx")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1_host.fastq")
					os.system("STAR --runThreadN "+n_thread+" --genomeDir "+path_genome+" --readFilesIn "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1_host.fastq --outFileNamePrefix "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+" --outReadsUnmapped Fastx")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					path=os.path.isfile(split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					if path:
						print("STAR finished without error\n")
						output_log.write("STAR finished without error\n")
					else:
						print("Error: cannot create STAR file\n")
						output_log.write("Error: cannot create STAR file\n")
						flag_ok="0"

				else :
				
					# Running "STAR" for alignment on reference genome
					
					print("Running STAR for alignment on reference genome")
					output_log.write("Running STAR for alignment on reference genome\n")
					os.system("mkdir -p "+split_output_path+"/3_STAR_Alignment")
					path_genome="/home/Genome/"
					os.system("STAR --runThreadN "+n_thread+" --genomeDir "+path_genome+" --readFilesIn "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq --outFileNamePrefix "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+" --outReadsUnmapped Fastx")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					path=os.path.isfile(split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					if path:
						print("STAR finished without error\n")
						output_log.write("STAR finished without error\n")
					else:
						print("Error: cannot create STAR file\n")
						output_log.write("Error: cannot create STAR file\n")
						flag_ok="0"

				if split_shotgun == "yes" :

					print("Running Shotgun Analysis")
					output_log.write("Running Shotgun Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Shotgun_Module_Results")
					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					else:
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq --confidence "+split_confidence)   
		
					path=os.path.isfile(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for viruses
					
					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa
					
					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq") 
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq --confidence "+split_confidence) 

					path=os.path.isfile(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/7_Kaiju_Shotgun_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_Shotgun.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create Kaiju file\n")	
						output_log.write("Error: cannot create Kaiju file\n")
						flag_ok="0"

					#Results 4_Kraken2
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_Shotgun.txt") 
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_Shotgun.txt")
					
						
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_Shotgun.txt")

					file_o1=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_Shotgun.txt")

					

					input=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
							if (int(students3[j]) < int(students3[j+1])): 
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
					plt.savefig(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_Shotgun.png")
					input.close()


				if split_denovo == "yes" :		
					# Running "SPADES" 
					
					assembly_de_novo_module_results
					print("Running DeNovo Analysis")
					output_log.write("Running DeNovo Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results")					

					print("Running SPADES")
					output_log.write("Running SPADES\n")
					os.system("mkdir -p "+split_output_path+"/5_SPADES_Assembly_DeNovo")
					if split_kmer == "auto" :
						os.system("spades.py -t "+n_thread+" --rna -s "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq  -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					else :
						os.system("spades.py -t "+n_thread+" -k "+split_kmer+" --rna -s "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq  -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])

					path=os.path.isfile(split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					if path:
						print("Spades finished without error\n")
						output_log.write("Spades finished without error\n")
					#else:
					#	print("Error: cannot create Spades file\n")
					#	output_log.write("Error: cannot create Spades file\n")
					#	flag_ok="0"

					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"
					

					# Running "kraken" for viruses

					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta")
					else:
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta --confidence "+split_confidence) 

					path=os.path.isfile(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa

					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta")	
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta --confidence "+split_confidence)	

					path=os.path.isfile(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/8_Kaiju_DeNovo_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_DeNovo.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create Kaiju file\n")
						output_log.write("Error: cannot create Kaiju file\n")
						flag_ok="0"

					#Results 6_Kraken2
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_DeNovo.txt")

					file_o1=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_DeNovo.txt")

					#piechart

					
					input=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_DeNovo.png")
					input.close()
		else:
			print("Pair-end mode")
			output_log.write("Pair-end mode\n")
			fastq=os.listdir(path_fastq)
			fastq.sort()
			#fastq.sort(lambda x, y: cmp(string.lower(x), string.lower(y)))
			for i in range(0,len(fastq),2):
				sample_only_name_a=fastq[i].split(".fastq")
				sample_only_name_b=fastq[i+1].split(".fastq")
				sample_only_name_a_noR1=sample_only_name_a[0].split("_R1")
				inputfile_a=path_fastq+"/"+fastq[i]
				inputfile_b=path_fastq+"/"+fastq[i+1]
			
				
				print("Sample in running: "+sample_only_name_a_noR1[0])
				output_log.write("Sample in running: "+sample_only_name_a_noR1[0]+"\n")
					
			   	# FASTQC before adapter trimming 
				
				if split_fastqc == "yes" :
					print("FASTQC before adapter trimming")
					output_log.write("FASTQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control")
					os.system("fastqc "+inputfile_a+" "+inputfile_b+" -o "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/"+" -t "+n_thread)
												
				# Running MultiQC before adapter trimming

				if split_fastqc == "yes" :
					print("MULTIQC before adapter trimming")	
					output_log.write("MULTIQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control/")	

				# Running trimmomatic for adapter trimming

				print("Running cutadapt for adapter trimming")
				output_log.write("Running cutadapt for adapter trimming\n")
				os.system("mkdir -p "+split_output_path+"/2_Fastq_trimmed")
				os.system("mkdir -p "+split_output_path+"/9_Results")
				#os.system("trimmomatic PE -threads "+n_thread+" -trimlog "+split_output_path+"/2_Fastq_trimmed/trimm_log_file.txt -summary "+split_output_path+"/2_Fastq_trimmed/trimm_summary_file.txt "+inputfile_a+" "+inputfile_b+" "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"unpaired_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"unpaired_trimmed.fastq ILLUMINACLIP:/opt/Anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 CROP:76")
				#nextereCTGTCTCTTATA 
				os.system("cutadapt -m 20 -a "+adapter+" -A "+adapter+" -j "+n_thread+" -o "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -p "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq "+inputfile_a+" "+inputfile_b)

				path=os.path.isfile(split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq")
				if path:
					print("Adapter trimming finished without error\n")
					output_log.write("Adapter trimming finished without error\n")
				else:
					print("Error: cannot create adapter trimmed file\n")
					output_log.write("Error: cannot create adapter trimmed file\n")
					flag_ok="0"
				                                                                           
				# FASTQC before adapter trimming 
				
				if split_fastqc2 == "yes" :
					print("FASTQC after adapter trimming")
					output_log.write("FASTQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control")
					os.system("fastqc "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -o "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/"+" -t "+n_thread)
							
				# Running MultiQC before adapter trimming

				if split_fastqc2 == "yes" :
					print("MULTIQC after adapter trimming")	
					output_log.write("MULTIQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control/")

				if split_delete == "yes" :
					print("Running STAR for alignment on reference genome and delete host organism")
					output_log.write("Running STAR for alignment on reference genome and delete host organism\n")
					os.system("mkdir -p "+split_output_path+"/3_STAR_Alignment")
					path_genome="/home/Genome/"
					os.system("STAR --runThreadN "+n_thread+" --genomeDir "+path_genome+" --readFilesIn "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq --outFileNamePrefix "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+" --outReadsUnmapped Fastx")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1_host.fastq")
					os.system("STAR --runThreadN "+n_thread+" --genomeDir "+path_genome+" --readFilesIn "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1_host.fastq --outFileNamePrefix "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+" --outReadsUnmapped Fastx")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq")

				else :
				
					# Running "STAR" for alignment on reference genome
					
					print("Running STAR for alignment on reference genome")
					output_log.write("Running STAR for alignment on reference genome\n")
					os.system("mkdir -p "+split_output_path+"/3_STAR_Alignment")
					path_genome="/home/Genome/"
					os.system("STAR --runThreadN "+n_thread+" --genomeDir "+path_genome+" --readFilesIn "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq --outFileNamePrefix "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+" --outReadsUnmapped Fastx")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq")
					os.system("mv "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq")

					path=os.path.isfile(split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1")
					if path:
						print("STAR finished without error\n")
						output_log.write("STAR finished without error\n")
					else:
						print("Error: cannot create STAR file\n")
						output_log.write("Error: cannot create STAR file\n")
						flag_ok="0"

				if split_shotgun == "yes" :

					print("Running Shotgun Analysis")
					output_log.write("Running Shotgun Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Shotgun_Module_Results")
										
					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq")   
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq --confidence "+split_confidence) 

					path=os.path.isfile(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"


					# Running "kraken" for viruses
					
					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa
					
					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq") 
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/7_Kaiju_Shotgun_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq -j "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_Shotgun.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create kaiju file\n")
						output_log.write("Error: cannot create kaiju file\n")
						flag_ok="0"

					#Results 4_Kraken2
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_Shotgun.txt") 
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_Shotgun.txt")
					
						
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_Shotgun.txt")

					file_o1=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")	
						else:
							output_o.write(line_o1)		
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_Shotgun.txt")

					

					input=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_Shotgun.png")
					input.close()


				if split_denovo == "yes" :		
					# Running "SPADES" 
					
					print("Running DeNovo Analysis")
					output_log.write("Running DeNovo Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results")

					print("Running SPADES")
					output_log.write("Running SPADES\n")
					os.system("mkdir -p "+split_output_path+"/5_SPADES_Assembly_DeNovo")
					if split_kmer == "auto" :
						os.system("spades.py --rna -t "+n_thread+" -1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq -2 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					else :
						os.system("spades.py -k "+split_kmer+" --rna -t "+n_thread+" -1 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate1.fastq -2 "+split_output_path+"/3_STAR_Alignment/"+sample_only_name_a_noR1[0]+"Unmapped.out.mate2.fastq -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])

					path=os.path.isfile(split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					if path:
						print("SPADES finished without error\n")
						output_log.write("Running cutadapt for adapter trimming\n")
					#else:
					#	print("Error: cannot create SPADES file\n")
					#	output_log.write("Running cutadapt for adapter trimming\n")
					#	flag_ok="0"


					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"
					

					# Running "kraken" for viruses

					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa

					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta")	
					else: 	
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/8_Kaiju_DeNovo_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/transcripts.fasta -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_DeNovo.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create kaiju file\n")
						output_log.write("Error: cannot create kaiju file\n")
						flag_ok="0"

					#Results 6_Kraken2
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")	
						else:
							output_o.write(line_o1)		
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_DeNovo.txt")

					file_o1=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_DeNovo.txt")

					#piechart

					
					input=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_DeNovo.png")
					input.close()
	else :
		print("DNA")
		output_log.write("DNA\n")
		#pairend
		if split_pairend == "no":
			print("Single-end mode")
			output_log.write("Single-end mode\n")
			fastq=os.listdir(path_fastq)
			fastq.sort()
			#fastq.sort(lambda x, y: cmp(string.lower(x), string.lower(y)))
			for i in range(0,len(fastq),1):
				sample_only_name_a=fastq[i].split(".fastq")
				#sample_only_name_b=fastq[i+1].split(".fastq")
				sample_only_name_a_noR1=sample_only_name_a[0].split("_R1")
				inputfile_a=path_fastq+"/"+fastq[i]
				#inputfile_b=path_fastq+"/"+fastq[i+1]
				
				print("Sample in running: "+sample_only_name_a_noR1[0])	
				output_log.write("Sample in running: "+sample_only_name_a_noR1[0]+"\n")

			   	# FASTQC before adapter trimming 
				
				if split_fastqc == "yes" :
					print("FASTQC before adapter trimming")
					output_log.write("FASTQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control")
					os.system("fastqc "+inputfile_a+" -o "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/"+" -t "+n_thread)
												
				# Running MultiQC before adapter trimming

				if split_fastqc == "yes" :
					print("MULTIQC before adapter trimming")
					output_log.write("MULTIQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control/")	

				# Running trimmomatic for adapter trimming

				print("Running cutadapt for adapter trimming")
				output_log.write("Running cutadapt for adapter trimming\n")
				os.system("mkdir -p "+split_output_path+"/2_Fastq_trimmed")
				os.system("mkdir -p "+split_output_path+"/9_Results")
				#os.system("trimmomatic PE -threads "+n_thread+" -trimlog "+split_output_path+"/2_Fastq_trimmed/trimm_log_file.txt -summary "+split_output_path+"/2_Fastq_trimmed/trimm_summary_file.txt "+inputfile_a+" "+inputfile_b+" "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"unpaired_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"unpaired_trimmed.fastq ILLUMINACLIP:/opt/Anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 CROP:76")
				#nextereCTGTCTCTTATA 
				os.system("cutadapt -m 20 -a "+adapter+" -j "+n_thread+" -o "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+inputfile_a)

				path=os.path.isfile(split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq")
				if path:
					print("Adapter trimming finished without error\n")
					output_log.write("Adapter trimming finished without error\n")
				else:
					print("Error: cannot create adapter trimmed file\n")
					output_log.write("Error: cannot create adapter trimmed file\n")
					flag_ok="0"
				                                                                           
				# FASTQC before adapter trimming 
				
				if split_fastqc2 == "yes" :
					print("FASTQC after adapter trimming")
					output_log.write("FASTQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control")
					os.system("fastqc "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -o "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/"+" -t "+n_thread)
							
				# Running MultiQC before adapter trimming

				if split_fastqc2 == "yes" :
					print("MULTIQC after adapter trimming")	
					output_log.write("MULTIQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control/")

				if split_delete == "yes" :
					print("Running bowtie2 for alignment on reference genome and delete host organism")
					output_log.write("Running bowtie2 for alignment on reference genome and delete host organism\n")
					os.system("mkdir -p "+split_output_path+"/3_Bowtie2_Alignment")
					os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_metrics_file_host.txt --un "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align_host.fastq --al "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_align_host.fastq -x "+path_delete+""+base_name+" -U "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -S "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"host.sam")
					os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_metrics_file.txt --un "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq --al "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_align.fastq -x "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align_host.fastq -S "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+".sam")

				else :
				
					# Running "bowtie2" for alignment on reference genome
					
					print("Running bowtie2 for alignment on reference genome")
					output_log.write("Running bowtie2 for alignment on reference genome\n")
					os.system("mkdir -p "+split_output_path+"/3_Bowtie2_Alignment")
					os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_metrics_file.txt --un "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq --al "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_align.fastq -x "+path_genome+""+base_name+" -U "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -S "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+".sam")

				path=os.path.isfile(split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+".sam")
				if path:
					print("Bowtie2 finished without error\n")
					output_log.write("Bowtie2 finished without error\n")
				else:
					print("Error: cannot create bowtie2 file\n")
					output_log.write("Error: cannot create bowtie2 file\n")
					flag_ok="0"

				if split_shotgun == "yes" :

					print("Running Shotgun Analysis")
					output_log.write("Running Shotgun Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Shotgun_Module_Results")
					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq")   
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for viruses
					
					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq")
					else:
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq --confidence "+split_confidence) 

					path=os.path.isfile(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa
					
					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq") 
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq --confidence "+split_confidence) 

					path=os.path.isfile(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/7_Kaiju_Shotgun_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_Shotgun.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create kaiju file\n")
						output_log.write("Error: cannot create kaiju file\n")
						flag_ok="0"

					#Results 4_Kraken2
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_Shotgun.txt") 
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_Shotgun.txt")
					
						
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_Shotgun.txt")

					file_o1=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_Shotgun.txt")

					

					input=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_Shotgun.png")
					input.close()


				if split_denovo == "yes" :		
					# Running "SPADES" 
					
					print("Running DeNovo Analysis")
					output_log.write("Running DeNovo Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results")

					print("Running SPADES")
					output_log.write("Running SPADES\n")
					os.system("mkdir -p "+split_output_path+"/5_SPADES_Assembly_DeNovo")
					if split_kmer == "auto" :
						os.system("spades.py -t "+n_thread+" -s "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					else :
						os.system("spades.py -k "+split_kmer+" -t "+n_thread+" -s "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_reads_not_align.fastq -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])

					path=os.path.isfile(split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					if path:
						print("Spades finished without error\n")
						output_log.write("Spades finished without error\n")
					#else:
					#	print("Error: cannot create Spades file\n")
					#	output_log.write("Error: cannot create Spades file\n")
					#	flag_ok="0"


					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					#else:
					#	print("Error: cannot create Kraken file\n")
					#	output_log.write("Error: cannot create Kraken file\n")
					#	flag_ok="0"
					

					# Running "kraken" for viruses

					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa

					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta")		
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/8_Kaiju_DeNovo_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_DeNovo.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile("config.txt")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create kaiju file\n")
						output_log.write("Error: cannot create kaiju file\n")
						flag_ok="0"

					#Results 6_Kraken2
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_DeNovo.txt")

					file_o1=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_DeNovo.txt")

					#piechart

					
					input=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_DeNovo.png")
					input.close()
		else:
			print("Pair-end mode")
			output_log.write("Pair-end mode\n")
			fastq=os.listdir(path_fastq)
			fastq.sort()
			#fastq.sort(lambda x, y: cmp(string.lower(x), string.lower(y)))
			for i in range(0,len(fastq),2):
				print(fastq)
				sample_only_name_a=fastq[i].split(".fastq")
				print(sample_only_name_a)
				sample_only_name_b=fastq[i+1].split(".fastq")
				print(sample_only_name_b)
				sample_only_name_a_noR1=sample_only_name_a[0].split("_R1")
				print(sample_only_name_a[0])
				inputfile_a=path_fastq+"/"+fastq[i]
				print(inputfile_a)
				inputfile_b=path_fastq+"/"+fastq[i+1]
				print(inputfile_b)
				
				print("Sample in running: "+sample_only_name_a_noR1[0])
				output_log.write("Sample in running: "+sample_only_name_a_noR1[0]+"\n")
	
			   	# FASTQC before adapter trimming 
				
				if split_fastqc == "yes" :
					print("FASTQC before adapter trimming")
					output_log.write("FASTQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control")
					os.system("fastqc "+inputfile_a+" "+inputfile_b+" -o "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/"+" -t "+n_thread)
												
				# Running MultiQC before adapter trimming

				if split_fastqc == "yes" :
					print("MULTIQC before adapter trimming")
					output_log.write("MULTIQC before adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_before_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_before_Trimming_Quality_Control/")	

				# Running trimmomatic for adapter trimming

				print("Running cutadapt for adapter trimming")
				output_log.write("Running cutadapt for adapter trimming\n")
				os.system("mkdir -p "+split_output_path+"/2_Fastq_trimmed")
				os.system("mkdir -p "+split_output_path+"/9_Results")
				#os.system("trimmomatic PE -threads "+n_thread+" -trimlog "+split_output_path+"/2_Fastq_trimmed/trimm_log_file.txt -summary "+split_output_path+"/2_Fastq_trimmed/trimm_summary_file.txt "+inputfile_a+" "+inputfile_b+" "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"unpaired_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq  "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"unpaired_trimmed.fastq ILLUMINACLIP:/opt/Anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 CROP:76")
				#nextereCTGTCTCTTATA 
				os.system("cutadapt -m 20 -a "+adapter+" -A "+adapter+" -j "+n_thread+" -o "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -p "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq "+inputfile_a+" "+inputfile_b)

				path=os.path.isfile(split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq")
				if path:
					print("Adapter trimming finished without error\n")
					output_log.write("Adapter trimming finished without error\n")
				else:
					print("Error: cannot create adapter trimmed file\n")
					output_log.write("Error: cannot create adapter trimmed file\n")
					flag_ok="0"
				                                                                           
				# FASTQC before adapter trimming 
				
				if split_fastqc2 == "yes" :
					print("FASTQC after adapter trimming")
					output_log.write("FASTQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control")
					os.system("fastqc "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -o "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/"+" -t "+n_thread)
							
				# Running MultiQC before adapter trimming

				if split_fastqc2 == "yes" :
					print("MULTIQC after adapter trimming")	
					output_log.write("MULTIQC after adapter trimming\n")
					os.system("mkdir -p "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control")
					os.system("multiqc "+split_output_path+"/1_FastQC_after_Trimming_Quality_Control/*.zip -o "+split_output_path+"/1_multiQC_after_Trimming_Quality_Control/")

				if split_delete == "yes" :
					print("Running bowtie2 for alignment on reference genome and delete host organism")
					output_log.write("Running bowtie2 for alignment on reference genome and delete host organism\n")
					os.system("mkdir -p "+split_output_path+"/3_Bowtie2_Alignment")
					os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_metrics_file_host.txt --un-conc "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align_host.fastq --al-conc "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_align_host.fastq --un "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_unpaired_reads_not_align_host.fastq --al "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_unpaired_reads_align_host.fastq -x "+path_delete+" -1 "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -2 "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -S "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"host.sam")
					os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_metrics_file.txt --un-conc "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.fastq --al-conc "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_align.fastq --un "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_unpaired_reads_not_align.fastq --al "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_unpaired_reads_align.fastq -x "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align_host_1.fastq -2 "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align_host_2.fastq -S "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+".sam")

				else :
				
					# Running "bowtie2" for alignment on reference genome
					
					print("Running bowtie2 for alignment on reference genome")
					output_log.write("Running bowtie2 for alignment on reference genome\n")
					os.system("mkdir -p "+split_output_path+"/3_Bowtie2_Alignment")
					os.system("bowtie2 -p "+n_thread+" --end-to-end --very-sensitive --met-file "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_metrics_file.txt --un-conc "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.fastq --al-conc "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_align.fastq --un "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_unpaired_reads_not_align.fastq --al "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_unpaired_reads_align.fastq -x "+path_genome+"host_genome_idx -1 "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_a[0]+"_trimmed.fastq -2 "+split_output_path+"/2_Fastq_trimmed/"+sample_only_name_b[0]+"_trimmed.fastq -S "+split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+".sam")

					path=os.path.isfile(split_output_path+"/3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+".sam")
					if path:
						print("Bowtie2 finished without error\n")
						output_log.write("Bowtie2 finished without error\n")
					else:
						print("Error: cannot create bowtie2 file\n")
						output_log.write("Error: cannot create bowtie2 file\n")
						flag_ok="0"

				if split_shotgun == "yes" :

					print("Running Shotgun Analysis")
					output_log.write("Running Shotgun Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Shotgun_Module_Results")
					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq")   
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq  --use-names --report "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for viruses
					
					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa
					
					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq") 
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq#.fastq --classified-out "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq#.fastq --use-names --report "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 --paired "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq --confidence "+split_confidence) 

					path=os.path.isfile(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/7_Kaiju_Shotgun_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq -j "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju//names.dmp -i "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_Shotgun.html "+split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create kaiju file\n")
						output_log.write("Error: cannot create kaiju file\n")
						flag_ok="0"

					#Results 4_Kraken2
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_Shotgun.txt") 
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_Shotgun.txt")
					
						
					file_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_Shotgun.txt")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_Shotgun.txt")
						
					file_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_Shotgun.txt")

					file_o1=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"7_Kaiju_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")	
						else:
							output_o.write(line_o1)		
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_Shotgun.txt")

					

					input=open(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_protozoa_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_bacteria_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_Shotgun.png")
					input.close()

					input=open(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/4_Kraken2_viruses_Shotgun_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Shotgun_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_Shotgun.png")
					input.close()


				if split_denovo == "yes" :		
					# Running "SPADES" 
					
					print("Running DeNovo Analysis")
					output_log.write("Running DeNovo Analysis\n")
					os.system("mkdir -p "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results")

					print("Running SPADES")
					output_log.write("Running SPADES\n")
					os.system("mkdir -p "+split_output_path+"/5_SPADES_Assembly_DeNovo")
					if split_kmer == "auto" :
						os.system("spades.py --meta -t "+n_thread+" -1 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq -2 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					else :
						os.system("spades.py --meta -k "+split_kmer+" -t "+n_thread+" -1 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.1.fastq -2 "+split_output_path+"3_Bowtie2_Alignment/"+sample_only_name_a_noR1[0]+"_paired_reads_not_align.2.fastq -o "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])

					path=os.path.isfile(split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0])
					if path:
						print("SPADES finished without error\n")
						output_log.write("SPADES finished without error\n")
					#else:
					#	print("Error: cannot create SPADES file\n")
					#	output_log.write("Error: cannot create SPADES file\n")
					#	flag_ok="0"

					# Running "kraken" for bacteria
					
					print("Running kraken for bacteria")
					output_log.write("Running kraken for bacteria\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_bacteria --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"
					

					# Running "kraken" for viruses

					print("Running kraken for viruses")
					output_log.write("Running kraken for viruses\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta")
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_viruses --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta --confidence "+split_confidence)

					path=os.path.isfile(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: ccannot create Kraken file\n")
						output_log.write("Error: ccannot create Kraken file\n")
						flag_ok="0"

					# Running "kraken" for protozoa

					print("Running kraken for protozoa")
					output_log.write("Running kraken for protozoa\n")
					os.system("mkdir -p "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation")
					if (split_confidence == "default") or (split_confidence == "0.5"):
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta")		
					else: 
						os.system("kraken2 --db Db_Kraken2_Kaiju_protozoa --threads "+n_thread+" --unclassified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_unclass_seq.fastq --classified-out "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_class_seq.fastq --use-names --report "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt --output "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt --confidence 0.5 "+split_output_path+"/5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta --confidence "+split_confidence)
					
					path=os.path.isfile(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_output_kraken2.txt")
					if path:
						print("Kraken finished without error\n")
						output_log.write("Kraken finished without error\n")
					else:
						print("Error: cannot create Kraken file\n")
						output_log.write("Error: cannot create Kraken file\n")
						flag_ok="0"	

					# Running "kaiju"

					print("Running kaiju")
					output_log.write("Running kaiju\n")
					os.system("mkdir -p "+split_output_path+"/8_Kaiju_DeNovo_Annotation")
					os.system("kaiju -z "+n_thread+" -v -t Db_Kaiju/nodes.dmp -f Db_Kaiju/nr/kaiju_db_nr.fmi -E 0.001 -i "+split_output_path+"5_SPADES_Assembly_DeNovo/"+sample_only_name_a_noR1[0]+"/contigs.fasta -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					os.system("kaiju2table -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -e -r species -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out") 
					os.system("kaiju2krona -t Db_Kaiju/nodes.dmp -n Db_Kaiju/names.dmp -i "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"krona.out.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")
					os.system("ktImportText -o "+split_output_path+"9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_protein_DeNovo.html "+split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.krona")

					path=os.path.isfile(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_kaiju_e_val-3.out")
					if path:
						print("Kaiju finished without error\n")
						output_log.write("Kaiju finished without error\n")
					else:
						print("Error: cannot create kaiju file\n")
						output_log.write("Error: cannot create kaiju file\n")
						flag_ok="0"

					#Results 6_Kraken2
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_bacteria_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_bacteria_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_virus_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_virus_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_annotation_protozoa_DeNovo.txt")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt","w")
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
					os.system("cp "+split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_tax_count_protozoa_DeNovo.txt")
						
					file_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","w")
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

					file_o1=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")	
						else:
							output_o.write(line_o1)		
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_bacteria_DeNovo.txt")

					file_o1=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2.txt","r")
					output_o=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt","w")

					line_o1=file_o1.readline()		
					split=line_o1.split("\n")
					output_o.write(split[0]+"\t Protein_validated\n")

					while(line_o1!=""):
						split1=line_o1.split()
						file_o2=open(split_output_path+"8_Kaiju_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"kaiju_summary.tsv","r")
						line_o2=file_o2.readline()
						trovato=0
						while(line_o2!=""):			
							split2=line_o2.split("\t")
							if split2[3]==split1[4]:
								trovato=1
							line_o2=file_o2.readline()
						if trovato ==1:
							split=line_o1.split("\n")
							output_o.write(split[0]+"\t Protein_validated\n")
						else:
							output_o.write(line_o1)			
						line_o1=file_o1.readline()	
						file_o2.close()
					
					file_o1.close()
					output_o.close()
					os.system("cp "+split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_proteinvalidated.txt "+split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_result_protein_validated_virus_DeNovo.txt")

					#piechart

					
					input=open(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_protozoa_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_protozoa_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_bacteria_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_bacteria_DeNovo.png")
					input.close()

					input=open(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_report_kraken2_filtered_count.txt","r")
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
					plt.savefig(split_output_path+"/6_Kraken2_viruses_DeNovo_Annotation/"+sample_only_name_a_noR1[0]+"_piechart.png")
					plt.savefig(split_output_path+"/9_Results/Assembly_De_Novo_Module_Results/"+sample_only_name_a_noR1[0]+"_HOMEBIO_abundance_piechart_virus_DeNovo.png")
					input.close()
	
	if flag_ok=="1":
		print("Analysis finished successfully \n")
		output_log.write("Analysis finished successfully \n")
	else:
		print("Analysis finished with Error: check log file for more information\n")
		output_log.write("Analysis finished with Error: check log file for more information\n")
			
			
	output_log.close()	
		
			
if __name__ == "__main__":
	main(sys.argv[1:])

