#!/share/data3/yangjunbo/miniminiconda/bin/python
__date__ = "2022-6-24"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "yangjunbo"

import re
import os
import time
from time import strftime
import optparse
from optparse import OptionParser
import pandas as pd

parser = OptionParser('Usage: %prog -i [input] -o [output] \n return a binary dict, which can loaded by pickle')
parser.add_option('-i','--input',
        dest='input',
        help='Input file: file with group information.')
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
output = open(options.out,"w")
for i in data:
    pre_line = i
    line = i.strip()
    i = i.strip().split("\t")
    group_information = i[-1]
    #print(group_information)
    subgroup = "none"
    general_pattern1 = re.compile('\s([1-9])\n')
    general_pattern2 = re.compile("T([1-9])")
    general_pattern3 = re.compile("type-([1-9])")
    general_pattern4 = re.compile('MRV-?_?([1-9])')
    if general_pattern1.search(pre_line):
        subgroup = general_pattern1.search(pre_line).group(1)
    elif general_pattern2.search(group_information):
        subgroup = general_pattern2.search(group_information).group(1)
    elif general_pattern3.search(group_information):
        subgroup = general_pattern3.search(group_information).group(1)
    elif general_pattern4.search(pre_line):
                subgroup = general_pattern4.search(pre_line).group(1)
    if subgroup == "1":
        output.write(line + "\t" + subgroup + "\t538120\n")
    elif subgroup == "2":
        output.write(line + "\t" + subgroup + "\t538121\n")
    elif subgroup == "3":
        output.write(line + "\t" + subgroup + "\t538123\n")
    elif subgroup == "4":
                output.write(line + "\t" + subgroup + "\t538122\n")
    elif subgroup == "5":
                output.write(line + "\t" + subgroup + "\t2730540\n")
    else:
        output.write(line + "\t" + subgroup + "\tnone\n")
data.close()
output.close()




