#!/share/data3/yangjunbo/miniminiconda/bin/python
"""This is used to check files in the input directory. move those files with zero size to the output directory, and record these file names into test."""

__date__ = "2022-6-24"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import os
import time
from time import strftime
import optparse
from optparse import OptionParser
import pandas as pd

parser = OptionParser('Usage: %prog -i [input] -o [output] \n This is used to check files in the input directory. \
			Move those files with zero size to the output directory, and record these file names into test.\n')
parser.add_option('-i','--input',
		dest='input',
		help='Input directory')
parser.add_option('-o','--out',
		dest='out',
		help='Output directory')

(options,args) = parser.parse_args()
import sys
from sys import argv
import math
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

Input = options.input
output = options.out
if os.path.isdir(Input):
	print("Input directory is OK!!!\n")
else:
	print("Input directory not exist!!! please check!!!\n")

if os.path.isdir(output):
	print("Output directory is OK!!!\n")
else:
	os.system("mkdir {}".format(output))

file_list = os.listdir(Input)
#print(file_list)
out = open(options.out+"/"+options.out,"w")
for i in file_list:
	size = os.path.getsize(options.input+"/"+i)
	if size:
		pass
	else:
		out.write(i+"\n")
		os.system("mv {}/{} {}".format(Input,i,output))
out.close()










