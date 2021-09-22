#!/usr/bin/python

import sys, getopt
import os
import string
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
	print("HOME-Bio")
	
	config = configparser.ConfigParser()

	config.read(configfile)

	# Reading Config File
	
	split_fastqc=config['CUSTOM']['quality control before trimming [yes/no]']
	if split_fastqc!='yes':
		split_fastqc=config['DEFAULT']['quality control before trimming [yes/no]']
	
	split_fastqc2=config['CUSTOM']['quality control after trimming [yes/no]']
	if split_fastqc2!='yes':
		split_fastqc2=config['DEFAULT']['quality control after trimming [yes/no]']

	n_thread=config['CUSTOM']['number of threads [n]']
	if n_thread=='':
		n_thread=config['DEFAULT']['number of threads [n]']

	path_fastq=config['CUSTOM']['path fastq']
	
	split_pairend=config['CUSTOM']['paired-end? [yes/no]']
	if split_pairend!='no':
		split_pairend=config['DEFAULT']['paired-end? [yes/no]']
	
	split_output_path=config['CUSTOM']['path output']
	os.system("mkdir "+split_output_path)

	split_path_genome=config['CUSTOM']['path host genome']
	path_genome=split_path_genome

	split_delete=config['CUSTOM']['filter out contaminant [yes/no]']
	if split_delete!='yes':
		split_delete=config['DEFAULT']['filter out contaminant [yes/no]']

	split_path_delete=config['CUSTOM']['path contaminant genome']
	path_delete=split_path_delete

	split_adapter=config['CUSTOM']['adapter?']
	if split_adapter=='':
		split_adapter=config['DEFAULT']['adapter?']
	adapter=split_adapter

	split_shotgun=config['CUSTOM']['shotgun module [yes/no]']
	if split_shotgun!='yes':
		split_shotgun=config['DEFAULT']['shotgun module [yes/no]']

	split_denovo=config['CUSTOM']['assembly denovo module [yes/no]']
	if split_denovo!='yes':
		split_denovo=config['DEFAULT']['assembly denovo module [yes/no]']

	split_path_kaiju=config['CUSTOM']['path kraken2 & kaiju databases']
	path_kaiju=split_path_kaiju
	



	cwd = os.getcwd()
	#os.seteuid(1000)
	os.system("""docker run -it --rm -v "{}\Script.py:/home/Script.py:ro" -v "{}\config_file.txt:/home/config_file.txt:ro" -v "{}\:/home/Input:ro" -v "{}\:/home/Output:rw" -v "{}\:/home/Genome:ro" -v "{}\KRAKENdb_bacteria:/home/Db_Kraken2_Kaiju_bacteria:ro" -v "{}\KRAKENdb_protozoa:/home/Db_Kraken2_Kaiju_protozoa:ro" -v "{}\KRAKENdb_viruses:/home/Db_Kraken2_Kaiju_viruses:ro" -v "{}\KAIJUdb:/home/Db_Kaiju:ro" biohaz/home_bio""" .format(cwd, cwd, path_fastq, split_output_path, path_genome, path_kaiju, path_kaiju, path_kaiju, path_kaiju))

		
			
			
if __name__ == "__main__":
	main(sys.argv[1:])


