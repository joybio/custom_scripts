#!/share/data3/yangjunbo/miniminiconda/bin/python
"""This is used to annotate the input file by the gb_decode dict."""

__date__ = "2022-7-19"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import pickle
import optparse
from optparse import OptionParser
import re

parser = OptionParser('Usage: %prog -i [input] -o [output] -v [virus name] -g [group dict]\n append subtype information to the raw data. \n')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: accession and title file.')

parser.add_option('-g','--group',
                dest='group',
                help='Input file: group file.')

parser.add_option('-o','--out',
                dest='out',
                help='Out file: Dictionary: key = NCBI ID/accession number. value = subgroup')


(options,args) = parser.parse_args()
import sys
from sys import argv
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

data = open(options.input,"r")
group = open(options.group,"rb")
out = open(options.out,"w")
group_dict = pickle.load(group)
for i in data:
	line = i.strip()
	i = i.split("\t")
	key = i[0]
	out.write(line)
	#print(key)
	if key in group_dict.keys():
		if len(group_dict[key]) > 2:
			out.write("\t" + group_dict[key][0] + "\t" + group_dict[key][1] +"\t" + group_dict[key][2] + "\n")
		else:
			out.write("\t" + group_dict[key][0] + "\t" + group_dict[key][1] + "\n")
	else:
		out.write("\t" + "undetected" + "\n")

data.close()
group.close()
out.close()





