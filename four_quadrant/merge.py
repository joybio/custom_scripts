#!/root/miniconda3/bin/python

"""Create a summary that shows log2FC of RNAseq and m6Aseq"""

__date__ = "2020-11-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#imports
import re
import os
import optparse
import bisect
#sort  package
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-r','--rna',
                dest='rna',
                help='log2FC of RNAseq')

parser.add_option('-m','--m6A',
                dest='m6A',
                help='log2FC of m6Aseq')

parser.add_option('-o','--out',
                dest='out',
                help='out file')
(options,args) = parser.parse_args()


RNA_dict = {}
with open(options.rna,"r") as d:
	d.readline()
	for i in d:
		i = i.strip().split("\t")
		key = i[0]
		value = float(i[1].replace("NA","0"))
		RNA_dict[key] = value
m6A_dict = {}
with open(options.m6A,"r") as f:
	f.readline()
	with open(options.out,"w") as o:
		o.write("Geneid\tRNA_fc\tCLIP_fc\n")
		for i in f:
			i = i.strip().split("\t")
			if i[0] not in m6A_dict.keys():
				m6A_dict[i[0]] = float(i[1])
			else:
				m6A_dict[i[0]] += float(i[1])
with open(options.out,"w") as o:
	o.write("Geneid\tRNA_fc\tCLIP_fc\n")
	for i in m6A_dict.keys():
		if i in RNA_dict.keys():
			o.write(i+"\t"+str(RNA_dict[i])+"\t"+str(m6A_dict[i])+"\n")

