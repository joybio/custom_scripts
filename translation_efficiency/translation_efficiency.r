library(tidyverse)
library(DESeq2)
library(edgeR)

#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF("/home/l/backup1/refgenome/homo_sapiens/Homo_sapiens.GRCh38.95.gtf",format="gtf")
#exons_gene <- exonsBy(txdb, by = "gene")
#exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
#write.csv(exons_gene_lens,"exons_gene_lens.csv")
exons_gene_lens <- read.csv("exons_gene_lens.csv",sep = "\t",stringsAsFactors = F)
#将，替换成制表符，并去掉引号
texons_gene_lens <- t(exons_gene_lens)
texons_gene_lens <- texons_gene_lens[-1,]
write.csv(texons_gene_lens,"texons_gene_lens.csv")
#mut_IP_24 <- read.table('24.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#mut_IP_38 <- read.table('38.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#mut_IP_WT1 <- read.table('Col-0.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#mut_IP_WT2 <- read.table('Col-0_2.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#reads_matrix <- cbind(mut_IP_24,mut_IP_38,mut_IP_WT1,mut_IP_WT2)
reads_matrix <- read.csv("riboCntTable.table",header = T,sep = "\t",row.names=1,stringsAsFactors = F)
#reads_matrix <- reads_matrix[,-1]
d <- DGEList(counts = reads_matrix)
d$genes$Length <- c(texons_gene_lens)
mycounts <- rpkm(d)
write.table(mycounts,"ribo_readscounts.RPKM.xls",quote= F,sep="\t",row.names = T)

reads_matrix2 <- read.csv("rnaCntTable.table",header = T,sep = "\t",row.names=1,stringsAsFactors = F)
#reads_matrix <- reads_matrix[,-1]
d <- DGEList(counts = reads_matrix2)
d$genes$Length <- c(texons_gene_lens)
mycounts <- rpkm(d)
write.table(mycounts,"rna_readscounts.RPKM.xls",quote= F,sep="\t",row.names = T)


library(riborex)
rna <- read.csv("rnaCntTable.table",head=T,sep="\t",row.names=1)
ribo <- read.csv("riboCntTable.table",head=T,sep="\t",row.names=1)
rnacond <- c("control", "control", "treated", "treated")
ribocond <- c("control", "control","treated", "treated")
res.deseq2 <- riborex(rna, ribo, rnacond, ribocond)
write.table(res.deseq2,"translation_efficiency.xls",quote= F,sep="\t",row.names = T)


ribo_rpkm <- read.csv("rna_readscounts.RPKM.xls",header = T,sep = "\t",stringsAsFactors = F)
rna_rpkm <- read.csv("ribo_readscounts.RPKM.xls",header = T,sep = "\t",stringsAsFactors = F)
TE <- read.csv("translation_efficiency.xls",header = T,sep = "\t",stringsAsFactors = F)
head(rpkm)
head(TE)

ribo_rpkm$gene_id <- row.names(ribo_rpkm)
rna_rpkm$gene_id <- row.names(rna_rpkm)
TE$gene_id <- row.names(TE)
rpkm <- merge(ribo_rpkm,rna_rpkm,by = c("gene_id"))
merge <- merge(rpkm,TE,by = c("gene_id"))

write.table(merge,"merge.RPKM.TE.xls",quote= F,sep="\t",row.names = F)


