#!/bin/python

"""export a dict: geneID：[bed format]
It is used for the result of prokka (proksryotes or virus)
"""
import json

"""bed info: chrom  start   stop gene_ID:geneName  type  strand"""

__date__ = "2023-4-12"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "MIT"

"""
The MIT License (MIT)

Copyright (c) 2022 Junbo Yang <yang_junbo_hi@126.com> <1806389316@pku.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager
from optparse import OptionParser
import random
from collections import defaultdict
import re
import os
import sys
from pathlib import Path


def argsParse():
    parser = OptionParser('Usage: %prog -i genome.fa -o GC')

    parser.add_option('-i', '--input',
                      dest='input',
                      help='Input file: gff file of prokka.')

    parser.add_option('-o', '--out',
                      dest='out',
                      help='prefix of output file.')

    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Input file must be specified !!!")
        sys.exit(1)
    elif options.out is None:
        parser.print_help()
        print("Output file must be specified !!!")
        sys.exit(1)
    return parser.parse_args()


def build_dict(Input):
    gene_ID_pattern = re.compile('ID=(\w*);')
    ##  .*?; stop until find ";"
    gene_name_pattern = re.compile('Name=(\w*);')
    Input_dict = {}
    with open(Input, "r") as f:
        for line in f:
            if line.startswith("#"):
                pass
            elif line.startswith(">"):
                break
            else:
                i = line.strip().split("\t")
                chrom = i[0]
                start = i[3]
                stop = i[4]
                gene_type = i[2]
                strand = i[6]
                information = ''.join(i[8:])
                if gene_ID_pattern.search(information):
                    gene_ID = gene_ID_pattern.search(information).group(1)
                else:
                    gene_ID = "NA"
                if gene_name_pattern.search(information):
                    gene_Name = gene_name_pattern.search(information).group(1)
                else:
                    gene_Name = "NA"
                Input_dict[gene_ID] = [chrom, start, stop, gene_ID + ":" + gene_Name, gene_type, strand]
    return Input_dict


if __name__ == "__main__":
    options, args = argsParse()
    out_dict = build_dict(options.input)
    with open(options.out + '.json', "w") as fj:
        json.dump(dict(out_dict), fj, indent=4)
    with open(options.out + '.bed', "w") as fbed:
        for i in out_dict.keys():
            fbed.write('\t'.join(map(str, out_dict[i])) + "\n")
