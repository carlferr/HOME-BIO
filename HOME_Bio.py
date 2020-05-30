#!/usr/bin/python

import sys, getopt
import os
import string

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
	print("HOME-Bio")
	

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
	split_path_kaiju=line_config.split()
	path_kaiju=split_path_kaiju[0]
	
	file_config.close()	

	cwd = os.getcwd()

	os.system("docker run -it --rm -v "+cwd+"/Script.py:/home/Script.py:ro -v "+cwd+"/config_file.txt:/home/config_file.txt:ro -v "+path_fastq+"/:/home/Input:ro -v "+split_output_path[0]+"/:/home/Output:rw -v "+path_genome+":/home/Genome:ro -v "+path_kaiju+"/KRAKENdb_bacteria:/home/Db_Kraken2_Kaiju_bacteria:ro -v "+path_kaiju+"/KRAKENdb_protozoa:/home/Db_Kraken2_Kaiju_protozoa:ro -v "+path_kaiju+"/KRAKENdb_viruses:/home/Db_Kraken2_Kaiju_viruses:ro -v "+path_kaiju+"/KAIJUdb:/home/Db_Kaiju:ro biohaz/home_bio")

			
			
			
if __name__ == "__main__":
	main(sys.argv[1:])


