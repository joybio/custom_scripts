#!/bin/bash
# two steps for sequence extraction
#/share/data3/yangjunbo/database/accession_taxid.dict
#/share/data3/yangjunbo/database/virus_dict

cat virus_taxid | grep Enterovirus > Enterovirus
cat Enterovirus | grep -E  "Enterovirus A|Enterovirus B|Enterovirus C|Enterovirus D" > EnterovirusA_D.acc
cat EnterovirusA_D.acc | grep -E  "Enterovirus A" > Enterovirus_A.acc
cat EnterovirusA_D.acc | grep -E  "Enterovirus B" > Enterovirus_B.acc
cat EnterovirusA_D.acc | grep -E  "Enterovirus C" > Enterovirus_C.acc
cat EnterovirusA_D.acc | grep -E  "Enterovirus D" > Enterovirus_D.acc

ls *.acc | while read id;do(extract_value_from_dict.py -i $id -d /share/data3/yangjunbo/database/accession_taxid.dict -o ${id}.ID &);done

ls *.ID |while read id;do(extract_virus_seq_by_acc.py -i $id -o ${id}.fa &);done

