#!/root/miniconda3/bin/python

"""extract.protein.py """

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
parser.add_option('-i','--in',
				dest='input',
				help='input file')
parser.add_option('-g','--gtf',
				dest='gtf',
				help='/home/l/backup1/refgenome/Arabidopsis/protein/Arabidopsis_thaliana.TAIR10.pep.all.fa')
parser.add_option('-t','--temp',
				dest='temp',
				help='temp file')
parser.add_option('-o','--out',
				dest='out',
				help='count and percentage file')
(options,args) = parser.parse_args()

data = open(options.input,'r')
ref = open(options.gtf,'r')
count= open(options.out,'w')
temp = open(options.temp,'w')
ref_dict={}
for i in ref:
	i = i.strip()
	if i.startswith(">"):
		i = i.lstrip('>').split('.')
		name = i[0]
		ref_dict[name] = ''
	else:
		ref_dict[name] += i 
ref.close()
total=0
for i in data:
	i = i.strip().split(',')
	if i[0] in ref_dict.keys():
		temp.write(ref_dict[i[0]])
		total = total + 1
	else:
		pass
print("protein number:",total)
data.close()
temp.close()

num_dict = {}
aa_num = 0
count.write("aa\taa_num\taa_percentage\n")
with open(options.temp,'r') as temp:
	for i in temp:
		i = list(i)
		i_set = set(i)
		for  j in i_set:
			num_dict[j] = i.count(j)
			aa_num += num_dict[j]
		for k in i_set:	
			count.write(k + '\t' + str(num_dict[k]) + "\t" + str(round(num_dict[k],4)/aa_num) + "\n")


temp.close()
count.close()


