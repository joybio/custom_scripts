#!/bin/bash
cat ../addition_list | while read id;do(pre=${id//\(/\\(};name=${pre//\)/\\)};echo $name;cat /share/data3/yangjunbo/database/viruses.format.nt.fa| grep -E "`echo $name`" >> mapped_id);done
'''
#####################################################################################################################
#只运行一次
#python scripts/prepare_pickle.py -i /share/data/public/DataBase/update_version/NT/viruses/viruses.nt.fa -o virus_dict
#####################################################################################################################
'''
python scripts/extract_virus_seq.py -i ../addition_list  -m mapped_id -o Total_sequence.fa


