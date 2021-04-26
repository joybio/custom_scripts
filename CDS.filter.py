#!/root/miniconda3/bin/python

"""ref format: startswith ATG; endswith TAG TGA TAA """

__date__ = "2019-10-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#imports
import re
import os
import optparse
#import bisect
#sort  package
import linecache
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-i','--in',dest='input',
			help='input file')

parser.add_option('-o','--out',	dest='out',
			help='count and percentage file')
(options,args) = parser.parse_args()

ref = open(options.input,"r")
out = open(options.out,"w")

ref_format = open("ref_format","w")
for i in ref:
	i = i.strip()
	if i.startswith(">"):
		ref_format.write("\n" + i + "\n")
	else:
		ref_format.write(i)
ref_format.close()
ref.close()

ref_dict = {}
with open("ref_format","r") as f:
	for i in f:
		i = i.strip()
		if i.startswith(">"):
			i = i.strip().split("\t")
			name = i[0]
			ref_dict[name] = ''
		else:
			num = len(i)%3
			if num == 0:
				if i.startswith("ATG"):
					if i.endswith("TAG"):
						ref_dict[name] += i
					elif i.endswith("TAA"):
						ref_dict[name] += i
					elif i.endswith("TAG"):
						ref_dict[name] += i
					else:
						pass
				else:
					pass
f.close()
os.system("rm ref_format")
for i in ref_dict.keys():
	if ref_dict[i] != '':
		out.write(i + '\n' + ref_dict[i] + '\n')
out.close()

