acc=$1
out=$2
mkdir -p $out
#/data/weikf/bin/miniconda3/envs/quast/bin/efetch -db nuccore -format fasta  -id $acc  >$out/$acc.fa
/share/data/weikf/bin/efetch -db nuccore -format fasta  -id $acc  >$out/$acc.fa
