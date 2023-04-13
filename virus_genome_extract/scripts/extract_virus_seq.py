#!/share/data1/lisong/Software/anaminiminiconda3/bin/python


"""
extract virus genome by provide virus name
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

parser = OptionParser('Usage: %prog  -i [input] -m [mapped_id] -v [virus_dict] -o [output]')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: virus name, which should provide by user, one virus per line.')

parser.add_option('-m','--mapped',
                dest='mapped',
                help='mapped_id, output of the first step. see run.sh')

parser.add_option('-v','--virus',
                dest='virus',
		default = "virus_dict",
                help='raw data from: /share/data/public/DataBase/update_version/NT/viruses/viruses.nt.fa;\
		virus hash: keys = virus name; values = virus sequence.')

parser.add_option('-o','--out',
                dest='out',
                help='Out file: Total virus genome.fa')

(options,args) = parser.parse_args()
print("INFO {} Start: Load virus genome......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
###viruses need to dowload###
data = open(options.input,"r")
###viruses genomes in NT database###
mapped = open(options.mapped,"r")
mapped_list = []
for m in mapped:
	mapped_list.append(m)
virus_dict = open(options.virus,"rb")
###virus genome hash###
seq_dict = pickle.load(virus_dict)
###Total genome seq###
out = open(options.out,"w")
###number of viruses need to download###
n = 1
###Result file###
D = "Results"
os.system("mkdir {}".format(D))
###logfile###
log_genome_virus = open("genome_virus_log.txt","w")
###export non genome virus###
non_genome_virus = open("non_genome_virus.txt","w")
print("INFO {} Thread 1: Start analysis......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))

for i in data:
	line = i.strip()
	virus = line.replace(" ","_").replace("/","_").replace("(","").replace(")","")
	os.system("mkdir -p {}".format(D+ "/" +virus))
	log_genome_virus.write("####################\n"+line + ":\n####################\n\t")
	m = 0
	virus_id = D + "/" +virus+"/Total_" +virus + ".fa"
	with open(virus_id,"w") as v:
		for k in mapped_list:
			ID = k.strip()
			if  n == 1:
				out.write(k+seq_dict[ID] + "\n")
			mapped_ID = k.lstrip(">").split(".")
			NM_id = D + "/" +virus + "/" + mapped_ID[0] + ".fa"
			if re.search(line,k):
				log_genome_virus.write(mapped_ID[0] + "|")
				m += 1
			###Total genome in specific virus###
				v.write(k+seq_dict[ID] + "\n")
			###specifice sequence of specific virus### 
				with open(NM_id,'w') as f:
					f.write(k+seq_dict[ID])
			else:
				pass
	if m == 0:
		non_genome_virus.write(i)
	log_genome_virus.write("\n")
	n+=1
print("INFO {} Thread 1: END......".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))
print("#################\nVirus number: %s\n#################"%n)
data.close()
mapped.close()
virus_dict.close()
out.close()
log_genome_virus.close()
non_genome_virus.close()

