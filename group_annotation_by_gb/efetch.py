#!/share/data3/yangjunbo/miniminiconda/bin/python
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

parser = OptionParser('Usage: %prog -i [input] -o [output] \n return a binary dict, which can loaded by pickle')
parser.add_option('-i','--input',
                dest='input',
                help='Input file: list of accession number.')
parser.add_option('-n','--num',
                dest='num',
		default=40,
		type="int",
                help='number')
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

data = open(options.input,"r")
output = options.out

if os.path.isdir(output):
	print("Output directory is OK!!!\n")
else:
	os.system("mkdir {}".format(output))
print("INFO {} Start: Loading accession number......\n \n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))

accession = pd.read_csv(data)
row_number = accession.shape[0]

n = options.num
for idx,row in accession.iterrows():
	if idx < n:
		accession_id = row[0]
		gb_file = str(accession_id) + ".gb"
		os.system("efetch -db nuccore -format gb -id {} > {}/{} &".format(accession_id,output,gb_file))
	else:
		os.system("efetch -db nuccore -format gb -id {} > {}/{} &".format(accession_id,output,gb_file))
		time.sleep(20)
		n += 10

data.close()
print("INFO {} Done !!!\n \n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time()))))


