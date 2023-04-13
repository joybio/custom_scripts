#!/share/data3/yangjunbo/miniminiconda/bin/python
#utf:coding-8

"""This is used for decode genbank file.make a dict of transcript ID (key) and subgroup (value)."""

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
				help='Input directory, which contain gb files.')
parser.add_option('-o','--out',
				dest='out',
				help='Output dict. key = VERSION (transcript_id); value = subtype')

(options,args) = parser.parse_args()
import sys
from sys import argv
import math
import re
from Bio import SeqIO
import pickle
if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

Input = options.input
if os.path.isdir(Input):
	print("Input directory is OK!!!\n")
else:
	print("Input directory not exist!!! please check!!!\n")

output = open(options.out,"wb")
file_list = os.listdir(Input)
group_dict = {}
for i in file_list:
	f = Input + "/" + i
	print(f)
	for record in SeqIO.parse(f, "genbank"):
		for feature in record.features:
			'''
			if feature.type == "CDS":
				symbol = feature.qualifiers.get("gene", ["???"])[0]
				gene_id = feature.qualifiers.get("db_xref", ["???"])[0]
				gene_id = re.sub('GeneID:', '', gene_id)
				transcript_id = feature.qualifiers.get("coded_by", ["???"])[0]
				transcript_id = re.sub(':.*', '', transcript_id)
			'''
			if feature.type == "source":
				species_name = feature.qualifiers.get("organism", ["???"])[0]
				species_id = feature.qualifiers.get("db_xref", ["???"])[0]
				species_id = re.sub('taxon:', '', species_id)
				subtype = feature.qualifiers.get("note",["NA"])[0]
				serotype = feature.qualifiers.get("serotype",["NA"])[0]
				if subtype != "NA":
					subtype = re.sub("note:","",subtype)
				else:
					if serotype != "NA":
						serotype = re.sub("serotype=","",serotype)
					else:
						subtype = feature.qualifiers.get("strain",["NA"])[0]
						subtype = re.sub("strain=","",subtype)
			'''	
			if feature.type == "Region":
				cdd_id = feature.qualifiers.get("db_xref", ["???"])[0]
				cdd_id = re.sub('CDD:', '', cdd_id)
			'''
		if serotype != "NA":
			group_dict[record.id] = [species_name,subtype,serotype]
		else:
			group_dict[record.id] = [species_name,subtype]
		
pickle.dump(group_dict,output)
output.close()

#.qualifiers
#– 存储feature附加信息（Python字典）。键（key）为值（value）所存信息的单字简要描述，值为实际信息。比如，键为 “evidence” ，而值为 “computational (non-experimental)”。 这只是为了提醒人们注意，该feature没有被实验所证实（湿实验）。Note：为与GenBank/EMBL文件中的feature tables对应，规定.qualifiers 中值为字符串数组（即使只有一个字符串）。



