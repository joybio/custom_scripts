library(Vennerable)
library(VennDiagram)
#library(reshape2)
ribo_2020 <- read.csv("../translation_efficiency.csv",head=T,sep = ",",row.names=1)
ribo_2020 <- na.omit(ribo_2020)

TE_up_ribo_2020 <- subset(ribo_2020,ribo_2020$log2FoldChange >= 0.585 & ribo_2020$padj <= 0.05)
gene_TE_up_ribo_2020 <- row.names(TE_up_ribo_2020)
TE_down_ribo_2020 <- subset(ribo_2020,ribo_2020$log2FoldChange <= -0.585 & ribo_2020$padj <= 0.05)
gene_TE_down_ribo_2020 <- row.names(TE_down_ribo_2020)
ribo_2019_5 <- read.csv("/home/l/backup3/ZW/20190515/eclip/ECLIP/DESeq/translation_efficiency.csv",head=T,sep = ",",row.names=1)
ribo_2019_5 <- na.omit(ribo_2019_5)
TE_up_ribo_2019_5 <- subset(ribo_2019_5,ribo_2019_5$log2FoldChange >= 0.585 & ribo_2019_5$padj <= 0.05)
TE_down_ribo_2019_5 <- subset(ribo_2019_5,ribo_2019_5$log2FoldChange <= -0.585 & ribo_2019_5$padj <= 0.05)
gene_TE_up_ribo_2019_5 <- row.names(TE_up_ribo_2019_5)
gene_TE_down_ribo_2019_5 <- row.names(TE_down_ribo_2019_5)
TE_up_venn_plot <- venn.diagram(x=list(TE_2020=gene_TE_up_ribo_2020,TE_2019_5=gene_TE_up_ribo_2019_5),filename ="TE_up_venn_plot.png",imagetype ="png",lwd = 0.1,fill = c("#66C3A5", "#8DA1CB"),alpha = 0.6,label.col = "black",cex = 1,fontfamily = "Times",cat.col = c("#66C3A5", "#8DA1CB"),cat.cex = 1,cat.fontfamily = "Times",margin = 0.05,cat.dist = c(0.03, 0.03), cat.pos = c(-20, 20))

TE_dwon_venn_plot <- venn.diagram(x=list(TE_2020=gene_TE_down_ribo_2020,TE_2019_5=gene_TE_down_ribo_2019_5),filename ="TE_down_venn_plot.png",imagetype ="png",lwd = 0.1,fill = c("#FD8D62","#56B4E9"),alpha = 0.6,label.col = "black",cex = 1,fontfamily = "Times",cat.col = c("#FD8D62","#56B4E9"),cat.cex = 1,cat.fontfamily = "Times",margin = 0.05,cat.dist = c(0.03, 0.03), cat.pos = c(-20, 20))

RNA_2020 <- read.csv("../DEGs/All_results.csv",head=T,sep = ",",row.names=1)
RNA_2020 <- na.omit(RNA_2020)

RNA_up_RNA_2020 <- subset(RNA_2020,RNA_2020$log2FoldChange >= 0.585 & RNA_2020$padj <= 0.05)
gene_RNA_up_RNA_2020 <- row.names(RNA_up_RNA_2020)
RNA_down_RNA_2020 <- subset(RNA_2020,RNA_2020$log2FoldChange <= -0.585 & RNA_2020$padj <= 0.05)
gene_RNA_down_RNA_2020 <- row.names(RNA_down_RNA_2020)
RNA_2019_5 <- read.csv("/home/l/backup3/ZW/20190515/eclip/ECLIP/DESeq/DEGs/All_results.csv",head=T,sep = ",",row.names=1)
RNA_2019_5 <- na.omit(RNA_2019_5)
RNA_up_RNA_2019_5 <- subset(RNA_2019_5,RNA_2019_5$log2FoldChange >= 0.585 & RNA_2019_5$padj <= 0.05)
RNA_down_RNA_2019_5 <- subset(RNA_2019_5,RNA_2019_5$log2FoldChange <= -0.585 & RNA_2019_5$padj <= 0.05)
gene_RNA_up_RNA_2019_5 <- row.names(RNA_up_RNA_2019_5)
gene_RNA_down_RNA_2019_5 <- row.names(RNA_down_RNA_2019_5)

TE_RNA_2020_up_venn_plot <- venn.diagram(x=list(TE_up_2020=gene_TE_up_ribo_2020,RNA_up_2020=gene_RNA_up_RNA_2020),filename ="TE_RNA_2020_up_venn_plot.png",imagetype ="png",lwd = 0.1,fill = c("#66C3A5", "#8DA1CB"),alpha = 0.6,label.col = "black",cex = 1,fontfamily = "Times",cat.col = c("#66C3A5", "#8DA1CB"),cat.cex = 1,cat.fontfamily = "Times",margin = 0.05,cat.dist = c(0.03, 0.03), cat.pos = c(-20, 20))


TE_RNA_2020_down_venn_plot <- venn.diagram(x=list(TE_down_2020=gene_TE_down_ribo_2020,RNA_down_2020=gene_RNA_down_RNA_2020),filename ="TE_RNA_2020_down_venn_plot.png",imagetype ="png",lwd = 0.1,fill = c("#66C3A5", "#8DA1CB"),alpha = 0.6,label.col = "black",cex = 1,fontfamily = "Times",cat.col = c("#66C3A5", "#8DA1CB"),cat.cex = 1,cat.fontfamily = "Times",margin = 0.05,cat.dist = c(0.03, 0.03), cat.pos = c(-20, 20))


TE_RNA_2019_up_venn_plot <- venn.diagram(x=list(TE_up_2019=gene_TE_up_ribo_2019_5,RNA_up_2019=gene_RNA_up_RNA_2019_5),filename ="TE_RNA_2019_up_venn_plot.png",imagetype ="png",lwd = 0.1,fill = c("#66C3A5", "#8DA1CB"),alpha = 0.6,label.col = "black",cex = 1,fontfamily = "Times",cat.col = c("#66C3A5", "#8DA1CB"),cat.cex = 1,cat.fontfamily = "Times",margin = 0.05,cat.dist = c(0.03, 0.03), cat.pos = c(-20, 20))


TE_RNA_2019_down_venn_plot <- venn.diagram(x=list(TE_down_2019=gene_TE_down_ribo_2019_5,RNA_down_2019=gene_RNA_down_RNA_2019_5),filename ="TE_RNA_2019_down_venn_plot.png",imagetype ="png",lwd = 0.1,fill = c("#66C3A5", "#8DA1CB"),alpha = 0.6,label.col = "black",cex = 1,fontfamily = "Times",cat.col = c("#66C3A5", "#8DA1CB"),cat.cex = 1,cat.fontfamily = "Times",margin = 0.05,cat.dist = c(0.03, 0.03), cat.pos = c(-20, 20))







