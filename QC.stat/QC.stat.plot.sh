#!/bin/bash
ls *bam | while read id
do
samtools stats -@ 16 --reference /home/l/backup1/refgenome/mm/Mus_musculus.GRCm38.dna.toplevel.fa ${id} > ./stat/$(basename ${id} .bam).stat
done

ls ./*stat | while read id
do
plot-bamstats -p $(basename ${id} .stat) ${id}
done

ls *.bam | while read id; do(file=$(basename ${id} .bam);qualimap bamqc --java-mem-size=40G -gff /home/l/backup1/refgenome/mm/Mus_musculus.GRCm38.98.chr.gtf -nr 100000 -nw 500 -nt 16 -bam $id -outdir ./qualimap/${file});done

ls *.bam | while read id; do(qualimap rnaseq --java-mem-size=40G -bam $id -gtf /home/l/backup1/refgenome/mm/Mus_musculus.GRCm38.98.chr.gtf -oc ./qualimap/count.matrix -outdir ./qualimap/rnaseq/$(basename ${id} .bam) -pe -outformat PDF:HTML);done

qualimap multi-bamqc --java-mem-size=40G -r -d mapfile1.txt -gff /home/l/backup1/refgenome/mm/Mus_musculus.GRCm38.98.chr.gtf -outdir multi_bamqc -outformat PDF:HTML

awk ‘{print $1,$2}’ mapfile1.txt | while read a b; do qualimap comp-counts --java-mem-size=40G -bam $b -gtf /home/l/backup1/refgenome/mm/Mus_musculus.GRCm38.98.chr.gtf -out $a.count -pe; done

qualimap counts --java-mem-size=40G -c -d mapfile2.txt -outdir counts_result -outformat PDF:HTML

