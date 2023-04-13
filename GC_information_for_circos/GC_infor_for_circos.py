#!/bin/python

"""Stastics of GC information, GCskew / GCratio / Nratio"""

__date__ = "2023-4-11"
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
                      help='Input file: genome.fa')

    parser.add_option('-w', '--window',
                      dest='window',
                      default=200,
                      type="int",
                      help='Window length. Default: 200')

    parser.add_option('-o', '--out',
                      dest='out',
                      default="GC",
                      type="str",
                      help='Prefix of output. Default: GC.')

    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif options.input is None:
        parser.print_help()
        print("Input file must be specified !!!")
        sys.exit(1)
    count = 3  # times for the packages install
    while count:
        try:
            import pyfaidx  #
            print('Dependent package pyfaidx is OK.\nDpendent module pyfaidx is OK.')
            break
        except:
            print('Dependent package pyfaidx is not found!!! \n Start intalling ....')
            os.system('pip install pyfaidx')
            count -= 1
            continue
    return parser.parse_args()

def windows_dict(Input):
    window_dict = defaultdict(str)
    with open(Input, 'r') as f:
        for i in f:
            if i.startswith(">"):
                i = i.lstrip(">").strip()
                key = i.replace(":", "-")
            else:
                value = i.strip()
                window_dict[key] += value
    return window_dict

def windows_GC(Input, Output):
    GC_ratio_file = Output + ".GC_ratio"
    GC_skew_file = Output + ".GC_skew"
    N_file = Output + ".N_ratio"
    with open(GC_ratio_file, "w") as out:
        with open(GC_skew_file, "w") as out2:
            with open(N_file, "w") as Nout:
                for i in Input.keys():
                    bed = i.split("-")
                    GC_ratio = (list(Input[i]).count("G") + list(Input[i]).count("C"))/len(list(Input[i]))
                    GC_skew = (list(Input[i]).count("G") - list(Input[i]).count("C"))/\
                              (list(Input[i]).count("G") + list(Input[i]).count("C"))
                    N_ratio = list(Input[i]).count("N")/len(list(Input[i]))
                    out.write("\t".join(map(str, bed)) + "\t" + str(GC_ratio) + "\n")
                    out2.write("\t".join(map(str, bed)) + "\t" + str(GC_skew) + "\n")
                    Nout.write("\t".join(map(str, bed)) + "\t" + str(N_ratio) + "\n")

if __name__ == "__main__":
    options, args = argsParse()
    e1 = time.time()
    genome_size = options.input.split("/")[-1] + ".size"
    os.system("faidx {} -i chromsizes > {}".format(options.input, genome_size))
    genome_window = options.input.split("/")[-1] + ".window"
    os.system("bedtools makewindows -g {} -w {}  > {}".format(genome_size, options.window, genome_window))
    genome_window_fa = options.input.split("/")[-1] + ".window.fa"
    os.system("bedtools getfasta -fi {} -bed {} -fo {}".format(options.input, genome_window, genome_window_fa))
    genome_window_dict = windows_dict(genome_window_fa)
    windows_GC(genome_window_dict, options.out)
    e2 = time.time()
    print("INFO {} Total times: {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
                                           round(float(e2 - e1), 2)))





















