#!/root/miniconda3/bin/python

"""Create a summary that shows the assignment of enrichment peak to annotation features"""

__date__ = "2019-10-2"
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
parser.add_option('-c','--csv',
                dest='csv',
                help='ribosome.csv')

parser.add_option('-g','--tran',
                dest='tran',
                help='/home/l/backup1/refgenome/homo_sapiens/gene_id2trans_id')
parser.add_option('-o','--out',
                dest='out',
                help='out.csv')
(options,args) = parser.parse_args()

data = open(options.csv,"r")
out = open(options.out,'w')
ref = open(options.tran,"r")

geneiddict = {}
for i in ref:
	i = i.strip().split("\t")
	geneiddict[i[0]] = i[2]
ref.close()

for i in data:
	i = i.strip().split(",")
	if i[0] in geneiddict.keys():
		out.write(i[0] + "," + geneiddict[i[0]] + "\n")
data.close()
out.close()



