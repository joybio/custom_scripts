#coding:utf-8
rm(list=ls())
setwd("/Users/apple/Downloads/WCL/2020.7.2/DEseq/")
library(tidyverse)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(GenomicFeatures)
#txdb <- makeTxDbFromGFF("/Users/apple/Downloads/refgenome/Arabidopsis_thaliana.TAIR10.41.gtf",format="gtf")
#exons_gene <- exonsBy(txdb, by = "gene")
#exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
#write.csv(exons_gene_lens,"exons_gene_lens.csv")
exons_gene_lens <- read.csv("/Users/apple/Downloads/refgenome/exons_gene_lens.csv",stringsAsFactors = F)
texons_gene_lens <- t(exons_gene_lens)
texons_gene_lens <- texons_gene_lens[-1,]
write.csv(texons_gene_lens,"texons_gene_lens.csv")
#mut_IP_24 <- read.table('24.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#mut_IP_38 <- read.table('38.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#mut_IP_WT1 <- read.table('Col-0.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#mut_IP_WT2 <- read.table('Col-0_2.dup.count.txt',header = FALSE, sep = "\t", quote = "",stringsAsFactors =F)
#reads_matrix <- cbind(mut_IP_24,mut_IP_38,mut_IP_WT1,mut_IP_WT2)
reads_matrix <- read.csv("total_readscounts.format.txt",header = T,sep = "\t",stringsAsFactors = F)
rownames(reads_matrix) <- reads_matrix$Geneid
reads_matrix <- reads_matrix[,-1]
d <- DGEList(counts = reads_matrix)
d$genes$Length <- c(texons_gene_lens)
mycounts <- rpkm(d)
write.csv(mycounts,"total_readscounts.RPKM.txt")
condition <- factor(c(rep("treat",2),rep("control",2)), levels = c("control","treat"))
colData <- data.frame(row.names=colnames(mycounts), condition)
write.csv(colData,"condition.csv")
dds <- DESeqDataSetFromMatrix(reads_matrix, colData, design= ~ condition)
dds <- DESeq(dds)
res= results(dds)
#图中x轴为baseMean，y轴为logFoldChange，若padj<0.1点标为红色。
pdf("plotMA.pdf")
plotMA(res)
dev.off()

res = res[order(res$pvalue),]
write.csv(res,file="All_results.csv")
table(res$padj<0.05)
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file="diff_gene_deseq2.csv")
#对单个基因在各组中的读数进行可视化。计数经归一化。
d = plotCounts(dds,gene = which.min(res$padj),intgroup = "condition",returnData = TRUE)
ggplot(d,aes(x = condition,y = count)) +
  geom_point(position = position_jitter(w = 0.1,h = 0)) +
  scale_y_log10()
ggsave("plotCounts.pdf")
dev.off()

# using rlog transformed data:
pdf("plotPCA.rlog.pdf")
rld <- rlog(dds)
plotPCA(rld)
dev.off()

# also possible to perform custom transformation:
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
pdf("plotPCA.custom.pdf")
plotPCA(DESeqTransform(se))
dev.off()

library(pheatmap)
pdf(file="sample.heatmap.pdf")
vsd = vst(dds,blind = FALSE)
sampleDist = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleDist)
rownames(sampleDistMatrix) = paste(vsd$condition,vsd$type,sep='_')
colnames(sampleDistMatrix) = NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist)
dev.off()

library(clusterProfiler)
library(topGO)
library(pathview)
library(org.At.tair.db)
library(biomaRt)
ALL = as.data.frame(subset(res, padj < 0.05))
DEG_ALL = as.character(row.names(ALL))
up = as.data.frame(subset(res, padj < 0.05 & log2FoldChange > 1))
DEG_up = as.character(row.names(up))
down = as.data.frame(subset(res, padj < 0.05 & log2FoldChange < -1))
DEG_down = as.character(row.names(up))

DEG_ALL.entrez_id = mapIds(x = org.At.tair.db,keys = DEG_ALL,keytype = "TAIR",column = "ENTREZID")
DEG_ALL.entrez_id = na.omit(DEG_ALL.entrez_id)

DEG_up.entrez_id = mapIds(x = org.At.tair.db,keys = DEG_up,keytype = "TAIR",column = "ENTREZID")
DEG_up.entrez_id = na.omit(DEG_up.entrez_id)

DEG_down.entrez_id = mapIds(x = org.At.tair.db,keys = DEG_down,keytype = "TAIR",column = "ENTREZID")
DEG_down.entrez_id = na.omit(DEG_down.entrez_id)

dir.create("ALL")
setwd("ALL/")
enrich.go_ALL.ALL = enrichGO(gene = DEG_ALL.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_ALL.BP = enrichGO(gene = DEG_ALL.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_ALL.CC = enrichGO(gene = DEG_ALL.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_ALL.MF = enrichGO(gene = DEG_ALL.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)

ALL.ALL.df = as.data.frame(enrich.go_ALL.ALL)
ALL.BP.df = as.data.frame(enrich.go_ALL.BP)
ALL.CC.df = as.data.frame(enrich.go_ALL.CC)
ALL.MF.df = as.data.frame(enrich.go_ALL.MF)
ALL.BP.df$Class <- rep("biological_process")
ALL.CC.df$Class <- rep("cellular_component")
ALL.MF.df$Class <- rep("molecular_function")
ALL.BP.df <- ALL.BP.df[order(-ALL.BP.df$Count),]
ALL.CC.df <- ALL.CC.df[order(-ALL.CC.df$Count),]
ALL.MF.df <- ALL.MF.df[order(-ALL.MF.df$Count),]
ALL.BP.df <- head(ALL.BP.df,n=10)
ALL.CC.df <- head(ALL.CC.df,n=10)
ALL.MF.df <- head(ALL.MF.df,n=10)
GO <- rbind(ALL.BP.df,ALL.CC.df,ALL.MF.df)
CPCOLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
#CPCOLS<-c("#CC6666", "#9999CC", "#66CC99")
#CPCOLS<-c("#999999", "#E69F00", "#56B4E9")
dorder = factor(as.character(GO$Description),levels = rev(as.character(GO$Description))) #对数据按Description排序
pdf("GO_annotation_all.pdf")
p <- ggplot(GO,aes(x=Description,y=Count,fill=Class)) #定义X轴，Y轴的数据和颜色填充 
p + geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder)) +  #定义柱形图的宽度和间距
  coord_flip() + #转换横纵坐标 
  scale_y_log10(breaks = c(1,10,100,1000)) + #Y轴按log10排序，展示1，10，100，1000刻度
  scale_fill_manual(values = CPCOLS)
#scale_fill_discrete(name="Ontology") + #修改legend的tittle
theme(panel.background = element_rect(fill = "transparent",colour = NA)) + #清楚背景颜色
  xlab("Term") #修改X轴标签
dev.off()

ALL.BP.df$Rich_factor <- (ALL.BP.df$Count)/sum(ALL.BP.df$Count)
ALL.BP.df <- ALL.BP.df[order(-ALL.BP.df$Count),]
ALL.BP.df.ten <- head(ALL.BP.df,n=10)
B=as.data.frame(ALL.BP.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_BP")
p+theme_bw()
ggsave("enrichmentgo_all_bp.bubble.pdf")

ALL.BP.df.fifteen <- head(ALL.BP.df,n=15)
ALL.BP.df.fifteen$padjust = -log10(ALL.BP.df.fifteen$pvalue)
B=as.data.frame(ALL.BP.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_all_bp.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#66C3A5",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

ALL.CC.df$Rich_factor <- (ALL.CC.df$Count)/sum(ALL.CC.df$Count)
ALL.CC.df <- ALL.CC.df[order(-ALL.CC.df$Count),]
ALL.CC.df.ten <- head(ALL.CC.df,n=10)
B=as.data.frame(ALL.CC.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_CC")
p+theme_bw()
ggsave("enrichmentgo_all_CC.bubble.pdf")

ALL.CC.df.fifteen <- head(ALL.CC.df,n=15)
ALL.CC.df.fifteen$padjust = -log10(ALL.CC.df.fifteen$pvalue)
B=as.data.frame(ALL.CC.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_all_CC.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#8DA1CB",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()


ALL.MF.df$Rich_factor <- (ALL.MF.df$Count)/sum(ALL.MF.df$Count)
ALL.MF.df <- ALL.MF.df[order(-ALL.MF.df$Count),]
ALL.MF.df.ten <- head(ALL.MF.df,n=10)
B=as.data.frame(ALL.MF.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_MF")
p+theme_bw()
ggsave("enrichmentgo_all_MF.bubble.pdf")

ALL.MF.df.fifteen <- head(ALL.MF.df,n=15)
ALL.MF.df.fifteen$padjust = -log10(ALL.MF.df.fifteen$pvalue)
B=as.data.frame(ALL.MF.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_all_MF.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

ALL_KEGG = enrichKEGG(gene = DEG_ALL,organism = "ath",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
ALL.KEGG.df = as.data.frame(ALL_KEGG)
ALL.KEGG.df$Rich_factor <- (ALL.KEGG.df$Count)/sum(ALL.KEGG.df$Count)
ALL.KEGG.df <- ALL.KEGG.df[order(-ALL.KEGG.df$Count),]
ALL.KEGG.df.ten <- head(ALL.KEGG.df,n=10)
B=as.data.frame(ALL.KEGG.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_KEGG")
p+theme_bw()
ggsave("enrichmentgo_all_KEGG.bubble.pdf")

ALL.KEGG.df.fifteen <- head(ALL.KEGG.df,n=15)
ALL.KEGG.df.fifteen$padjust = -log10(ALL.KEGG.df.fifteen$pvalue)
B=as.data.frame(ALL.KEGG.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_all_KEGG.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

dir.create("GO_raw_analysis")
setwd("GO_raw_analysis/")
pdf(file = "enrich.go_ALL.All_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_ALL.ALL, showCategor=20)
dev.off()
pdf(file = "enrich.go_ALL.BP_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_ALL.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_ALL.CC_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_ALL.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_ALL.MF_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_ALL.MF, showCategory=20)
dev.off()
pdf(file = "enrich.go_ALL.ALL_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_ALL.ALL, showCategory=20)
dev.off()
pdf(file = "enrich.go_ALL.BP_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_ALL.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_ALL.CC_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_ALL.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_ALL.MF_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_ALL.MF, showCategory=20)
dev.off()

setwd("../")
pdf(file = "enrich.go_ALL.BP_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_ALL.BP)
dev.off()
pdf(file = "enrich.go_ALL.CC_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_ALL.CC)
dev.off()
pdf(file = "enrich.go_ALL.MF_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_ALL.MF)
dev.off()
#cnetplot通路-基因网络图对于富集到的GO terms之间的基因重叠关系进行展示，如果两个GO terms系的差异基因存在重叠，说明这两个节点存在overlap关系
pdf("ALL_cnetplot.pdf")
cnetplot(enrich.go_ALL.BP, categorySize="pvalue",showCategory = 3,vertex.label.cex=0.1,colorEdge = TRUE)
dev.off()
pdf("ALL_cnetplot.circular.pdf")
cnetplot(enrich.go_ALL.BP, circular = TRUE, colorEdge = TRUE,showCategory = 3,vertex.label.cex=0.1)
dev.off()
#热图
pdf("GO_ALL.heatmap.pdf",width = 30,height = 8)
heatplot(enrich.go_ALL.ALL)
dev.off()
#pathway-pathway Network
pdf("ALL_emapplot.pdf")
emapplot(enrich.go_ALL.ALL)
dev.off()
#GO通路网络图
pdf("GO_ALL.netpathway.pdf")
goplot(enrich.go_ALL.BP)
dev.off()


#KEGG的可视化
ALL_KEGG = enrichKEGG(gene = DEG_ALL,organism = "ath",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
pdf(file = "ALL_KEGG.pdf")
dotplot(ALL_KEGG, showCategory=30)
dev.off()

dir.create("ALL_pathview")
setwd("ALL_pathview/")
DEG.entrez_id <- DEG_ALL.entrez_id
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath04075',species="ath")
#Ribosome biogenesis in eukaryotes
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03008',species="ath")
#Ribosome
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03010',species="ath")
#Zeatin biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00908',species="ath")
#03015mRNA surveillance pathway
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03015',species="ath")
#04141Protein processing in endoplasmic reticulum
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath04141',species="ath")
#01230Biosynthesis of amino acids
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath01230',species="ath")
#00910Nitrogen metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00910',species="ath")
#01070Biosynthesis of plant hormones
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath01070',species="ath")
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00250',species="ath")
#00260  Glycine, serine and threonine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00260',species="ath")
#00270  Cysteine and methionine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00270',species="ath")
#00280  Valine, leucine and isoleucine degradation
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00280',species="ath")
#00290  Valine, leucine and isoleucine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00290',species="ath")
#00300  Lysine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00300',species="ath")
#00310  Lysine degradation
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00310',species="ath")
#00220  Arginine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00220',species="ath")
#00330  Arginine and proline metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00330',species="ath")
#00340  Histidine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00340',species="ath")
#00350  Tyrosine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00350',species="ath")
#00360  Phenylalanine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00360',species="ath")
#00380  Tryptophan metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00380',species="ath")
#00400  Phenylalanine, tyrosine and tryptophan biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00400',species="ath")
#00410  beta-Alanine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00410',species="ath")
#00430  Taurine and hypotaurine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00430',species="ath")
#00440  Phosphonate and phosphinate metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00440',species="ath")
#00450  Selenocompound metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00450',species="ath")
#00460  Cyanoamino acid metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00460',species="ath")
#00480  Glutathione metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00480',species="ath")

setwd("../../")
dir.create("up")
setwd("up/")
enrich.go_up.ALL = enrichGO(gene = DEG_up.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_up.BP = enrichGO(gene = DEG_up.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_up.CC = enrichGO(gene = DEG_up.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_up.MF = enrichGO(gene = DEG_up.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)

up.ALL.df = as.data.frame(enrich.go_up.ALL)
up.BP.df = as.data.frame(enrich.go_up.BP)
up.CC.df = as.data.frame(enrich.go_up.CC)
up.MF.df = as.data.frame(enrich.go_up.MF)
up.BP.df$Class <- rep("biological_process")
up.CC.df$Class <- rep("cellular_component")
up.MF.df$Class <- rep("molecular_function")
up.BP.df <- up.BP.df[order(-up.BP.df$Count),]
up.CC.df <- up.CC.df[order(-up.CC.df$Count),]
up.MF.df <- up.MF.df[order(-up.MF.df$Count),]
up.BP.df <- head(up.BP.df,n=10)
up.CC.df <- head(up.CC.df,n=10)
up.MF.df <- head(up.MF.df,n=10)
GO <- rbind(up.BP.df,up.CC.df,up.MF.df)
CPCOLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
#CPCOLS<-c("#CC6666", "#9999CC", "#66CC99")
#CPCOLS<-c("#999999", "#E69F00", "#56B4E9")
dorder = factor(as.character(GO$Description),levels = rev(as.character(GO$Description))) #对数据按Description排序
pdf("GO_annotation_up_all.pdf",width = 30,height = 8)
p <- ggplot(GO,aes(x=Description,y=Count,fill=Class)) #定义X轴，Y轴的数据和颜色填充 
p + geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder)) +  #定义柱形图的宽度和间距
  coord_flip() + #转换横纵坐标 
  scale_y_log10(breaks = c(1,10,100,1000)) + #Y轴按log10排序，展示1，10，100，1000刻度
  scale_fill_manual(values = CPCOLS)
#scale_fill_discrete(name="Ontology") + #修改legend的tittle
theme(panel.background = element_rect(fill = "transparent",colour = NA)) + #清楚背景颜色
  xlab("Term") #修改X轴标签
dev.off()

up.BP.df$Rich_factor <- (up.BP.df$Count)/sum(up.BP.df$Count)
up.BP.df <- up.BP.df[order(-up.BP.df$Count),]
up.BP.df.ten <- head(up.BP.df,n=10)
B=as.data.frame(up.BP.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_BP")
p+theme_bw()
ggsave("enrichmentgo_up_bp.bubble.pdf",width = 30,height = 8)

up.BP.df.fifteen <- head(up.BP.df,n=15)
up.BP.df.fifteen$padjust = -log10(up.BP.df.fifteen$pvalue)
B=as.data.frame(up.BP.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_up_bp.bar.pdf",width = 30, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#66C3A5",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

up.CC.df$Rich_factor <- (up.CC.df$Count)/sum(up.CC.df$Count)
up.CC.df <- up.CC.df[order(-up.CC.df$Count),]
up.CC.df.ten <- head(up.CC.df,n=10)
B=as.data.frame(up.CC.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_CC")
p+theme_bw()
ggsave("enrichmentgo_up_CC.bubble.pdf",width = 30,height = 8)

up.CC.df.fifteen <- head(up.CC.df,n=15)
up.CC.df.fifteen$padjust = -log10(up.CC.df.fifteen$pvalue)
B=as.data.frame(up.CC.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_up_CC.bar.pdf",width = 30, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#8DA1CB",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()


up.MF.df$Rich_factor <- (up.MF.df$Count)/sum(up.MF.df$Count)
up.MF.df <- up.MF.df[order(-up.MF.df$Count),]
up.MF.df.ten <- head(up.MF.df,n=10)
B=as.data.frame(up.MF.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_MF")
p+theme_bw()
ggsave("enrichmentgo_up_MF.bubble.pdf",width = 30,height = 8)

up.MF.df.fifteen <- head(up.MF.df,n=15)
up.MF.df.fifteen$padjust = -log10(up.MF.df.fifteen$pvalue)
B=as.data.frame(up.MF.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_up_MF.bar.pdf",width = 30, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

up_KEGG = enrichKEGG(gene = DEG_up,organism = "ath",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
up.KEGG.df = as.data.frame(up_KEGG)
up.KEGG.df$Rich_factor <- (up.KEGG.df$Count)/sum(up.KEGG.df$Count)
up.KEGG.df <- up.KEGG.df[order(-up.KEGG.df$Count),]
up.KEGG.df.ten <- head(up.KEGG.df,n=10)
B=as.data.frame(up.KEGG.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_KEGG")
p+theme_bw()
ggsave("enrichmentgo_up_KEGG.bubble.pdf")

up.KEGG.df.fifteen <- head(up.KEGG.df,n=15)
up.KEGG.df.fifteen$padjust = -log10(up.KEGG.df.fifteen$pvalue)
B=as.data.frame(up.KEGG.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_up_KEGG.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

dir.create("GO_raw_analysis")
setwd("GO_raw_analysis/")
pdf(file = "enrich.go_up.ALL_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_up.up, showCategor=20)
dev.off()
pdf(file = "enrich.go_up.BP_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_up.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_up.CC_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_up.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_up.MF_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_up.MF, showCategory=20)
dev.off()
pdf(file = "enrich.go_up.ALL_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_up.ALL, showCategory=20)
dev.off()
pdf(file = "enrich.go_up.BP_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_up.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_up.CC_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_up.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_up.MF_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_up.MF, showCategory=20)
dev.off()
setwd("../")
pdf(file = "enrich.go_up.BP_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_up.BP)
dev.off()
pdf(file = "enrich.go_up.CC_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_up.CC)
dev.off()
pdf(file = "enrich.go_up.MF_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_up.MF)
dev.off()
#cnetplot通路-基因网络图对于富集到的GO terms之间的基因重叠关系进行展示，如果两个GO terms系的差异基因存在重叠，说明这两个节点存在overlap关系
pdf("up_cnetplot.pdf")
cnetplot(enrich.go_up.BP, categorySize="pvalue",showCategory = 3,vertex.label.cex=0.1,colorEdge = TRUE)
dev.off()
pdf("up_cnetplot.circular.pdf")
cnetplot(enrich.go_up.BP, circular = TRUE, colorEdge = TRUE,showCategory = 3,vertex.label.cex=0.1)
dev.off()
#热图
pdf("GO_up.heatmap.pdf",width = 30,height = 8)
heatplot(enrich.go_up.ALL)
dev.off()
#pathway-pathway Network
pdf("up_emapplot.pdf")
emapplot(enrich.go_up.ALL)
dev.off()
#GO通路网络图
pdf("GO_up.netpathway.pdf")
goplot(enrich.go_up.BP)
dev.off()


#KEGG的可视化
up_KEGG = enrichKEGG(gene = DEG_up,organism = "ath",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
pdf(file = "up_KEGG.pdf")
dotplot(up_KEGG, showCategory=30)
dev.off()

dir.create("up_pathview")
setwd("up_pathview/")
DEG.entrez_id <- DEG_up.entrez_id
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath04075',species="ath")
#Ribosome biogenesis in eukaryotes
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03008',species="ath")
#Ribosome
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03010',species="ath")
#Zeatin biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00908',species="ath")
#03015mRNA surveillance pathway
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03015',species="ath")
#04141Protein processing in endoplasmic reticulum
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath04141',species="ath")
#01230Biosynthesis of amino acids
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath01230',species="ath")
#00910Nitrogen metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00910',species="ath")
#01070Biosynthesis of plant hormones
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath01070',species="ath")
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00250',species="ath")
#00260  Glycine, serine and threonine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00260',species="ath")
#00270  Cysteine and methionine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00270',species="ath")
#00280  Valine, leucine and isoleucine degradation
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00280',species="ath")
#00290  Valine, leucine and isoleucine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00290',species="ath")
#00300  Lysine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00300',species="ath")
#00310  Lysine degradation
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00310',species="ath")
#00220  Arginine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00220',species="ath")
#00330  Arginine and proline metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00330',species="ath")
#00340  Histidine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00340',species="ath")
#00350  Tyrosine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00350',species="ath")
#00360  Phenylalanine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00360',species="ath")
#00380  Tryptophan metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00380',species="ath")
#00400  Phenylalanine, tyrosine and tryptophan biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00400',species="ath")
#00410  beta-Alanine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00410',species="ath")
#00430  Taurine and hypotaurine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00430',species="ath")
#00440  Phosphonate and phosphinate metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00440',species="ath")
#00450  Selenocompound metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00450',species="ath")
#00460  Cyanoamino acid metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00460',species="ath")
#00480  Glutathione metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00480',species="ath")

setwd("../../")
dir.create("down")
setwd("down/")
enrich.go_down.ALL = enrichGO(gene = DEG_down.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont =
                                "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_down.BP = enrichGO(gene = DEG_down.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = 
                               "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_down.CC = enrichGO(gene = DEG_down.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
enrich.go_down.MF = enrichGO(gene = DEG_down.entrez_id,OrgDb = org.At.tair.db,keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)

down.ALL.df = as.data.frame(enrich.go_down.ALL)
down.BP.df = as.data.frame(enrich.go_down.BP)
down.CC.df = as.data.frame(enrich.go_down.CC)
down.MF.df = as.data.frame(enrich.go_down.MF)
down.BP.df$Class <- rep("biological_process")
down.CC.df$Class <- rep("cellular_component")
down.MF.df$Class <- rep("molecular_function")
down.BP.df <- down.BP.df[order(-down.BP.df$Count),]
down.CC.df <- down.CC.df[order(-down.CC.df$Count),]
down.MF.df <- down.MF.df[order(-down.MF.df$Count),]
down.BP.df <- head(down.BP.df,n=10)
down.CC.df <- head(down.CC.df,n=10)
down.MF.df <- head(down.MF.df,n=10)
GO <- rbind(down.BP.df,down.CC.df,down.MF.df)
CPCOLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
#CPCOLS<-c("#CC6666", "#9999CC", "#66CC99")
#CPCOLS<-c("#999999", "#E69F00", "#56B4E9")
dorder = factor(as.character(GO$Description),levels = rev(as.character(GO$Description))) #对数据按Description排序
pdf("GO_annotation_down_all.pdf",width = 30,height = 8)
p <- ggplot(GO,aes(x=Description,y=Count,fill=Class)) #定义X轴，Y轴的数据和颜色填充 
p + geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,aes(x=dorder)) +  #定义柱形图的宽度和间距
  coord_flip() + #转换横纵坐标 
  scale_y_log10(breaks = c(1,10,100,1000)) + #Y轴按log10排序，展示1，10，100，1000刻度
  scale_fill_manual(values = CPCOLS)
#scale_fill_discrete(name="Ontology") + #修改legend的tittle
theme(panel.background = element_rect(fill = "transparent",colour = NA)) + #清楚背景颜色
  xlab("Term") #修改X轴标签
dev.off()

down.BP.df$Rich_factor <- (down.BP.df$Count)/sum(down.BP.df$Count)
down.BP.df <- down.BP.df[order(-down.BP.df$Count),]
down.BP.df.ten <- head(down.BP.df,n=10)
B=as.data.frame(down.BP.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_BP")
p+theme_bw()
ggsave("enrichmentgo_down_bp.bubble.pdf",width = 30,height = 8)

down.BP.df.fifteen <- head(down.BP.df,n=15)
down.BP.df.fifteen$padjust = -log10(down.BP.df.fifteen$pvalue)
B=as.data.frame(down.BP.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_down_bp.bar.pdf",width = 30, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#66C3A5",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

down.CC.df$Rich_factor <- (down.CC.df$Count)/sum(down.CC.df$Count)
down.CC.df <- down.CC.df[order(-down.CC.df$Count),]
down.CC.df.ten <- head(down.CC.df,n=10)
B=as.data.frame(down.CC.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_CC")
p+theme_bw()
ggsave("enrichmentgo_down_CC.bubble.pdf",width = 30,height = 8)

down.CC.df.fifteen <- head(down.CC.df,n=15)
down.CC.df.fifteen$padjust = -log10(down.CC.df.fifteen$pvalue)
B=as.data.frame(down.CC.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_down_CC.bar.pdf",width = 30, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#8DA1CB",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

down.MF.df$Rich_factor <- (down.MF.df$Count)/sum(down.MF.df$Count)
down.MF.df <- down.MF.df[order(-down.MF.df$Count),]
down.MF.df.ten <- head(down.MF.df,n=10)
B=as.data.frame(down.MF.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_MF")
p+theme_bw()
ggsave("enrichmentgo_down_MF.bubble.pdf",width = 30,height = 8)

down.MF.df.fifteen <- head(down.MF.df,n=15)
down.MF.df.fifteen$padjust = -log10(down.MF.df.fifteen$pvalue)
B=as.data.frame(down.MF.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_down_MF.bar.pdf",width = 30, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

down_KEGG = enrichKEGG(gene = DEG_down,organism = "ath",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
down.KEGG.df = as.data.frame(down_KEGG)
down.KEGG.df$Rich_factor <- (down.KEGG.df$Count)/sum(down.KEGG.df$Count)
down.KEGG.df <- down.KEGG.df[order(-down.KEGG.df$Count),]
down.KEGG.df.ten <- head(down.KEGG.df,n=10)
B=as.data.frame(down.KEGG.df.ten)
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
p=ggplot(B,aes(B$Rich_factor,B$Description))
p=p+geom_point(aes(size=B$Count,color=B$pvalue))
p=p+scale_color_gradient(low="red",high="green")
p=p+labs(color="p.value",size="gene number",x="Rich Factor",y="GO term",title="GO_KEGG")
p+theme_bw()
ggsave("enrichmentgo_down_KEGG.bubble.pdf")

down.KEGG.df.fifteen <- head(down.KEGG.df,n=15)
down.KEGG.df.fifteen$padjust = -log10(down.KEGG.df.fifteen$pvalue)
B=as.data.frame(down.KEGG.df.fifteen)
B <- B[order(B$pvalue),]
B$Description = factor(as.character(B$Description),levels = rev(as.character(B$Description)))
pdf(file="enrichmentgo_down_KEGG.bar.pdf",width = 10, height = 8)
p<-ggplot(B,aes(x=B$Description,y=B$padjust))+geom_bar(stat='identity',fill="#FD8D62",position=position_dodge(0.6),width = 0.5) + labs(x="",y="-log10(pvalue)",title="GO Term") + theme_minimal()
p+coord_flip()+ theme(plot.margin = unit(rep(4,4),"lines"),
                      axis.text.x=element_text(hjust = 1,colour="black",family="Times",size=8), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
                      axis.text.y=element_text(family="Times",size=8,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
                      axis.title.y=element_text(family="Times",size = 8,face="plain"), #设置y轴标题的字体属性
                      panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
                      legend.position = "none"
)
#axis.line=element_line(colour="black"),
dev.off()

dir.create("GO_raw_analysis")
setwd("GO_raw_analysis/")
pdf(file = "enrich.go_down.All_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_down.ALL, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.BP_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_down.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.CC_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_down.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.MF_barplot.pdf",width = 30,height = 8)
barplot(enrich.go_down.MF, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.ALL_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_down.ALL, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.BP_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_down.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.CC_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_down.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.MF_dotplot.pdf",width = 30,height = 8)
dotplot(enrich.go_down.MF, showCategory=20)
dev.off()

setwd("../")
pdf(file = "enrich.go_down.BP_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_down.BP, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.CC_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_down.CC, showCategory=20)
dev.off()
pdf(file = "enrich.go_down.MF_plotGOgraph.pdf",width = 10,height = 8)
plotGOgraph(enrich.go_down.MF, showCategory=20)
dev.off()
#cnetplot通路-基因网络图对于富集到的GO terms之间的基因重叠关系进行展示，如果两个GO terms系的差异基因存在
#重叠，说明这两个节点存在overlap关系
pdf("down_cnetplot.pdf")
cnetplot(enrich.go_down.BP, categorySize="pvalue",showCategory = 3,vertex.label.cex=0.1,colorEdge = TRUE)
dev.off()
pdf("down_cnetplot.circular.pdf")
cnetplot(enrich.go_down.BP, circular = TRUE, colorEdge = TRUE,showCategory = 3,vertex.label.cex=0.1)
dev.off()
#热图
pdf("GO_down.heatmap.pdf",width = 30,height = 8)
heatplot(enrich.go_down.ALL)
dev.off()
#pathway-pathway Network
pdf("down_emapplot.pdf")
emapplot(enrich.go_down.ALL)
dev.off()
#GO通路网络图
pdf("GO_down.netpathway.pdf")
goplot(enrich.go_down.BP)
dev.off()


#KEGG的可视化
down_KEGG = enrichKEGG(gene = DEG_down,organism = "ath",keyType = "kegg",pvalueCutoff = 0.05,pAdjustMethod = 'BH',minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
pdf(file = "down_KEGG.pdf")
dotplot(down_KEGG, showCategory=30)
dev.off()


dir.create("down_pathview")
setwd("down_pathview/")
DEG.entrez_id <- DEG_down.entrez_id
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath04075',species="ath")
#Ribosome biogenesis in eukaryotes
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03008',species="ath")
#Ribosome
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03010',species="ath")
#Zeatin biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00908',species="ath")
#03015mRNA surveillance pathway
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath03015',species="ath")
#04141Protein processing in endoplasmic reticulum
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath04141',species="ath")
#01230Biosynthesis of amino acids
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath01230',species="ath")
#00910Nitrogen metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00910',species="ath")
#01070Biosynthesis of plant hormones
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath01070',species="ath")
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00250',species="ath")
#00260  Glycine, serine and threonine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00260',species="ath")
#00270  Cysteine and methionine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00270',species="ath")
#00280  Valine, leucine and isoleucine degradation
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00280',species="ath")
#00290  Valine, leucine and isoleucine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00290',species="ath")
#00300  Lysine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00300',species="ath")
#00310  Lysine degradation
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00310',species="ath")
#00220  Arginine biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00220',species="ath")
#00330  Arginine and proline metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00330',species="ath")
#00340  Histidine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00340',species="ath")
#00350  Tyrosine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00350',species="ath")
#00360  Phenylalanine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00360',species="ath")
#00380  Tryptophan metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00380',species="ath")
#00400  Phenylalanine, tyrosine and tryptophan biosynthesis
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00400',species="ath")
#00410  beta-Alanine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00410',species="ath")
#00430  Taurine and hypotaurine metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00430',species="ath")
#00440  Phosphonate and phosphinate metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00440',species="ath")
#00450  Selenocompound metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00450',species="ath")
#00460  Cyanoamino acid metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00460',species="ath")
#00480  Glutathione metabolism
pathview(gene.data = DEG.entrez_id, pathway.id = 'ath00480',species="ath")

