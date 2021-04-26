#!/bin/bash
mkdir stat
ls *bam | while read id
do
samtools stats -@ 16 --reference /home/l/backup1/refgenome/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa ${id} > ./stat/$(basename ${id} .bam).stat
done

ls stat/*stat | while read id
do
plot-bamstats -p stat/$(basename $id '.stat')/ ${id}
done

ls *.bam | while read id; do(file=$(basename ${id} .bam);qualimap bamqc --java-mem-size=40G -gff /home/l/backup1/refgenome/homo_sapiens/Homo_sapiens.GRCh38.101.gtf -nr 100000 -nw 500 -nt 16 -bam $id -outdir ./qualimap/${file});done

ls *.bam | while read id; do(qualimap rnaseq --java-mem-size=40G -bam $id -gtf /home/l/backup1/refgenome/homo_sapiens/Homo_sapiens.GRCh38.101.gtf -oc ./qualimap/count.matrix -outdir ./qualimap/rnaseq/$(basename ${id} .bam) -pe -outformat PDF:HTML);done

qualimap multi-bamqc --java-mem-size=40G -r -d mapfile1.txt -gff /home/l/backup1/refgenome/homo_sapiens/Homo_sapiens.GRCh38.101.gtf -outdir multi_bamqc -outformat PDF:HTML

awk ‘{print $1,$2}’ mapfile1.txt | while read a b; do qualimap comp-counts --java-mem-size=40G -bam $b -gtf /home/l/backup1/refgenome/homo_sapiens/Homo_sapiens.GRCh38.101.gtf -out $a.count -pe; done

qualimap counts --java-mem-size=40G -c -d mapfile2.txt -outdir counts_result -outformat PDF:HTML

