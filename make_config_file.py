#!/usr/bin/python
import sys
import os.path

output=open("config.txt","w")

print("Quality Control before Trimming [yes/no] default -> no")
output.write("Quality Control before Trimming [yes/no] default -> no\n")
output.write(sys.stdin.readline())
output.write("\n")


print("Quality Control after Trimming [yes/no] default -> no")
output.write("Quality Control after Trimming [yes/no] default -> no\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Number of threads [n] default -> 1")
output.write("Number of threads [n] default -> 1\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Path Fastq")
output.write("Path Fastq\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Paired-end? [yes/no] default -> yes")
output.write("Paired-end? [yes/no] default -> yes\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Path output")
output.write("Path output\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Path host genome")
output.write("Path host genome\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Filter out contaminant [yes/no] default -> no")
output.write("Filter out contaminant [yes/no] default -> no\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Path contaminant genome (Not used if the previous answer is no)")
output.write("Path contaminant genome (Not used if the previous answer is no\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Adapter?")
output.write("Adapter?\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Shotgun module [yes/no] default -> no")
output.write("Shotgun module [yes/no] default -> no\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Assembly DeNovo module [yes/no] default -> no")
output.write("Assembly DeNovo module [yes/no] default -> no\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Path Kraken2 & Kaiju databases")
output.write("Path Kraken2 & Kaiju databases\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Kraken2 confidence [FLOAT 0 to 1] default -> 0.5")
output.write("Kraken2 confidence\n")
output.write(sys.stdin.readline())
output.write("\n")

print("Nucleic acid type [DNA/RNA] default -> DNA")
output.write("Nucleic acid type [DNA/RNA] default -> DNA\n")
output.write(sys.stdin.readline())
output.write("\n")

print("k-mer size used in Assembly DeNovo module (must be odd and less than 128) default -> auto")
output.write("k-mer [auto/21,33,55,77]\n")
output.write(sys.stdin.readline())
output.write("\n")

output.close()
path=os.path.isfile("config.txt")
if path:
	print("Config file created correctly\n")
else:
	print("Error: cannot create config file\n")

