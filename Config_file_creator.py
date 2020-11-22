import configparser
config = configparser.ConfigParser()
config['DEFAULT'] = {'Quality Control before Trimming [yes/no]': 'no',
                     'Quality Control after Trimming [yes/no]': 'no',
                     'Number of threads [n]': '1',
                     'Paired-end? [yes/no]': 'yes',
                     'Filter out contaminant [yes/no]': 'no',
                     'Adapter?': 'CTGTCTCTTATA',
                     'Shotgun module [yes/no]': 'no',
                     'Assembly DeNovo module [yes/no]': 'no',
                     'Kraken2 confidence': '0.5',
                     'Nucleic acid type [DNA/RNA]': 'DNA',
                     'k-mer [auto/21,33,55,77]': 'auto'}
config['CUSTOM'] = {'Quality Control before Trimming [yes/no]': 'no',
                    'Quality Control after Trimming [yes/no]': 'no',
                    'Number of threads [n]': '16',
                    'Path Fastq': '/root/Scrivania/HOME-BIO-master/fastq/',
                    'Paired-end? [yes/no]': 'yes',
                    'Path output': '/root/Scrivania/HOME-BIO-master/output/',
                    'Path host genome': '/root/Scrivania/HOME-BIO-master/Genome/',
                    'Filter out contaminant [yes/no]': 'no',
                    'Path contaminant genome': '/root/Scrivania/HOME-BIO-master/Genome/',
                    'Adapter?': 'CTGTCTCTTATA',
                    'Shotgun module [yes/no]': 'yes',
                    'Assembly DeNovo module [yes/no]': 'yes',
                    'Path Kraken2 & Kaiju databases': '/root/Scrivania/HOME-BIO-master/Db_Kraken2_Kaiju',
                    'Kraken2 confidence': '0.5',
                    'Nucleic acid type [DNA/RNA]': 'DNA',
                    'k-mer [auto/21,33,55,77]': 'auto'}
                    

with open('config_file.txt', 'w') as configfile:
   config.write(configfile)

