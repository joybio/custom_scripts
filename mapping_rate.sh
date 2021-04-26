#!/bin/bash
ls *.txt | while read id;do(echo $id >> mapping_rate.csv;cat $id | grep -E "in total" | awk -F ' ' '{print $1}'>> mapping_rate.csv;cat $id | grep "0 mapped" | awk -F ' ' '{print $1"\n"$5}'>> mapping_rate.csv;sed -i 's/(//g'  mapping_rate.csv);done
sed -i 'N;s/\n/,/g' mapping_rate.csv
sed -i 'N;s/\n/,/g' mapping_rate.csv
echo "sample,total_reads,mapped_reads,map_rate" > head.csv
cat head.csv mapping_rate.csv > mapping_stast.csv

