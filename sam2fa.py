#!/root/miniconda3/bin/python

"""extract unmapped reads from sam"""

__date__ = "2019-10-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
import bisect
#sort  package
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-s','--sam',
                dest='sam',
                help='sam file')

parser.add_option('-o','--out',
                dest='out',
                help='out fa')
(options,args) = parser.parse_args()

sam = open(options.sam,'r')
out = open(options.out,'w')

for i in sam:
	i=i.strip().split('\t')
	out.write(i[0] + '\n' + i[9] + '\n' + '+' + '\n' + i[10] + '\n')



