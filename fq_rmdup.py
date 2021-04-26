#!/bin/python
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
detected = []
unique = []
for rec in SeqIO.parse(open('inputfile.fastq', 'rU'), 'fastq'):
   cksum = seguid(rec.seq)
   if cksum not in detected:
       unique.append(rec)
       detected.append(cksum)
SeqIO.write(unique, open('deduplicated.fastq','w'), 'fastq')


