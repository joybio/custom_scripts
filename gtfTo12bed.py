#!/bin/python
#coding: utf-8 -*-
"""
Created on Wed Dec 16 11:03:05 2015
@author: biochen
"""
import sys, getopt, csv, re
opts, args = getopt.getopt(sys.argv[1:], "hi:o:")
for op, value in opts:
	if op == "-i":
		GTF_in_File_Name = value
	elif op == "-o":
		BED_in_File_Name = value
	elif op == "-h":
		print("Usage: python3 GTF2BED.py -i input.gtf -o output.bed")
		sys.exit()
GTF_in_file = open(GTF_in_File_Name, "r")
BED_out_file = open(BED_in_File_Name, "w")
GTFs = csv.reader(GTF_in_file, dialect = "excel-tab")
BEDs = csv.writer(BED_out_file, dialect = "excel-tab")
transcripts = {}
for GTF in GTFs:
	if GTF[0].startswith("#"):
		pass
#\w 匹配字母数字或者下划线；\s匹配任意空白字符；\d匹配数字
#\W 匹配非字母，非数字，非下划线；\S匹配非空白字符；\D匹配非数字
#[...] 匹配字符组中的字符
#[^...] 匹配除了字符组中字符的所有字符
	elif GTF[2] == "exon":
		transcript_id = re.findall('transcript_id "([^;]*)"', GTF[8])[0]
		if transcript_id in transcripts:
			transcripts[transcript_id][6].append([int(GTF[3]), int(GTF[4])])
		else:
			transcripts[transcript_id] = []
			transcripts[transcript_id].append(GTF[0])
			transcripts[transcript_id].append(GTF[3])
			transcripts[transcript_id].append(GTF[4])
			transcripts[transcript_id].append(GTF[5])
			transcripts[transcript_id].append(GTF[6])
			transcripts[transcript_id].append(transcript_id)
			transcripts[transcript_id].append([[int(GTF[3]), int(GTF[4])]])
for transcript in transcripts:
	transcripts[transcript][6].sort()
	transcript_start = transcripts[transcript][6][0][0]
	transcript_end = transcripts[transcript][6][len(transcripts[transcript][6]) - 1][1]
	exon_sizes = ""
	exon_starts = ""
	for exon in transcripts[transcript][6]:
		exon_size = exon[1] - exon[0] + 1
		exon_start = exon[0] - transcript_start
		exon_sizes = exon_sizes + str(exon_size) + ","
		exon_starts = exon_starts + str(exon_start) + ","
	BED = []
	BED.append(transcripts[transcript][0])
	BED.append(transcript_start - 1)
	BED.append(transcript_end)
	BED.append(transcripts[transcript][5])
	BED.append(transcripts[transcript][3])
	BED.append(transcripts[transcript][4])
	BED.append(transcript_start - 1)
	BED.append(transcript_end)
	BED.append("0, 0, 0")
	BED.append(len(transcripts[transcript][6]))
	BED.append(exon_sizes)
	BED.append(exon_starts)
	BEDs.writerow(BED)
GTF_in_file.close()
BED_out_file.close()


