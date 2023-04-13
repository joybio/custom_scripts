#!/share/data3/yangjunbo/miniconda3/bin/python

"""
extract virus genome by provide acc
"""
__date__ = "2022-6-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import optparse
import pickle
from optparse import OptionParser
import os
import re
import time
from time import strftime

parser = OptionParser('Usage: %prog  -i [input] -m [mapped_id] -v [virus_dict] -o [output] \n \
			Options: -n [number]')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: acc, which should provide by user, one accession number per line.')

parser.add_option('-v','--virus',
                dest='virus',
		default = "/share/data3/yangjunbo/database/virus_dict",
                help='raw data from: /share/data3/yangjunbo/database/virus_dict;\
		virus hash: keys = virus name; values = virus sequence.')

parser.add_option('-n','--num',
                dest='num',
                default = 1,
                help='column of NCBI_ID, which column is used to search nt database. Default: [1]')

parser.add_option('-o','--out',
                dest='out',
                help='Out file: Total virus genome.fa')

(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
print("INFO {} Start: Load virus genome......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
###viruses need to dowload###
data = open(options.input,"r")
###viruses genomes in NT database###
virus_dict = open(options.virus,"rb")
###virus genome hash###
seq_dict = pickle.load(virus_dict)
###Total genome seq###
outpath = str(options.out)
genome_virus = outpath + ".mapped.log.txt"
log_genome_virus = open(genome_virus,"w")
###export non genome virus###
non_genome_virus = outpath + ".unmapped.log.txt"
non_genome_virus = open(non_genome_virus,"w")
number = int(options.num)

out = open(options.out,"w")
print("INFO {} Thread 1: Start analysis......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
#number of mapped NCBI ID
n = 0
for i in data:
	line = i.strip()
	i = i.strip().split("\t")
	NCBI_ID = i[number]
	if NCBI_ID in seq_dict.keys():
		log_genome_virus.write(NCBI_ID +"\n")
		n += 1
		###Total genome in specific virus###
		out.write(seq_dict[NCBI_ID])
		###specifice sequence of specific virus### 
	#	NCBI_out = outpath + "_" +  NCBI_ID +".fa"
	#	with open(NCBI_out,'w') as f:
	#		f.write(">" + NCBI_ID + "\n" + seq_dict[NCBI_ID])
	else:
		non_genome_virus.write(NCBI_ID +"\n")
print("INFO {} Thread 1: END......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
print("#################\nmap to NCBI_number: %s\n#################"%n)
data.close()
virus_dict.close()
out.close()
log_genome_virus.close()
non_genome_virus.close()

