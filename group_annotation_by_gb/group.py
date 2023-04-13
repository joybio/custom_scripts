#!/share/data3/yangjunbo/miniminiconda/bin/python
__date__ = "2022-7-19"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import pickle
import optparse
from optparse import OptionParser
import re

parser = OptionParser('Usage: %prog -i [input] -o [output] \n extract annotation information by gb file \n \
			return a binary dict, which can loaded by pickle\n \
			key = accession number, value = group information')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: gb file.')

parser.add_option('-o','--out',
                dest='out',
                help='Out file: Dictionary: key = NCBI ID/accession number. value = subgroup')


(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

from collections import defaultdict

data = open(options.input,"r")
out = open(options.out,"wb")
gb_dict = {}
group_dict = defaultdict(list)
for i in data:
	line = i
	i = i.strip()
	if i.startswith("LOCUS"):
		i = re.split(r"[ ]+",i)
		key = i[1]
		#print(key)
		gb_dict[key] = ''
	else:
		gb_dict[key] += line
#print(gb_dict)
data.close()
for i in gb_dict.keys():
	value = gb_dict[i].split("\n")
	for k in value:
		if k.startswith("VERSION"):
			key = re.split(r"[ ]+",k)[1]
		k = k.split("/note=")
		if len(k) > 1:
			group_dict[key].append(k[1].replace("\"",""))
			break
			#print(key,k[1].replace("\"",""))
pickle.dump(group_dict,out)
out.close()











