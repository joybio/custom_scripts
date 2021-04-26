#!/bin/bash
ls *.UCSC.bowtie2.sam | while read id;do(grep "AS:" $id | grep -v "XS:" > $(basename $id .UCSC.bowtie2.sam).uniq.sam);done
ls *.uniq.sam | while read id;do(grep -e "YT:Z:CP" $id > $(basename $id '.uniq.sam').pair-end_unique.sam);done


