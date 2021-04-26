#!/bin/bash
awk '{print $1"\t"$3"\t"$1".count""\t""2"}' mapfile1.txt > mapfile2.txt



