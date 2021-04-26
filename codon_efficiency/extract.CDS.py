#!/root/miniconda3/bin/python

"""extract.cds.py """

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
		help='/home/l/backup1/refgenome/homo_sapiens/CDS/gene.human.CDS.format.fa')
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

#dict
for i in ref:
	i = i.strip()
	if i.startswith(">"):
		i = i.lstrip(">")
		i = i.lstrip("hg38_knownGene_").split('.')
		name = i[0]
		ref_dict[name] = ''
	else:
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
ref.close()
#search

total=0
for i in data:
	i = i.strip().split(',')
	if i[1] in ref_dict.keys():
		a =len(str(ref_dict[i[1]]))%3
		if a == 0:
			temp.write(ref_dict[i[1]])
			total = total + 1
	else:
		pass
print("CDS number:",total)
data.close()
temp.close()

num_dict = {}

def cut_text(text,lenth):
    textArr = re.findall('.{'+str(lenth)+'}', text)
    textArr.append(text[(len(textArr)*lenth):])
    return textArr
#def cut_text(text,length):
#	textarr = re.findall(r'.{3}', text)
#	return textarr
codon_num = 0
with open(options.temp,'r') as temp:
	for i in temp:
		i=i.strip()
		num_dict = {}
		code=cut_text(i,3)
		for j in code:
			if re.match("[ACGT]",j):
				num_dict[j] = code.count(j)
		for m in num_dict.keys():
			codon_num += num_dict[m]
		for k in num_dict.keys():
			if re.match("[ACGT]",k):
				count.write(k + '\t' + str(num_dict[k]) + "\t" +\
					str(round(num_dict[k],4)/codon_num) + "\n")

temp.close()
count.close()



