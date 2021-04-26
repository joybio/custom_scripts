#For ribosome occupancy analysis, the translation efficiency (TE) was calculated by dividing the CDS footprinting abundancy
#(excluding the starting 15 codons) by its mRNA abundancy. The genes with 1.5 fold change of TE were considered as the differential
#translation genes (Zinshteyn and Gilbert, 2013).
rm(list=ls())
library(ggplot2)
library(dplyr)
DATA=read.csv("merge.xls",header = T,sep='\t')
data=na.omit(DATA)
data_1 <- subset(data,data$RNA_fc >= 1 & data$CLIP_fc >= 1)
data_2 <- subset(data,data$RNA_fc >= 1 & data$CLIP_fc <= -1)
data_3 <- subset(data,data$RNA_fc <= -1 & data$CLIP_fc >= 1)
data_4 <- subset(data,data$RNA_fc <= -1 & data$CLIP_fc <= -1)
label <- rbind(data_1,data_2,data_3,data_4)
label$label <- label$Geneid
data_label <- data.frame(DATA=NULL)
row.names(data) <- data$Geneid
for(i in as.character(row.names(data))){
  if(i %in% as.character(label$Geneid)){
    data_label <- data_label
  }
  else{
    data_label <- rbind(data_label,data[i,])
    }
}
data_label$label <- rep("")
data <- rbind(label,data_label)

color.vec=rep("gray",length(data$Geneid))
color.vec[data$RNA_fc >= 1 & data$CLIP_fc >= 1] = "#E69F00"
color.vec[data$RNA_fc >= 1 & data$CLIP_fc <= -1] = "#8DA1CB"
color.vec[data$RNA_fc <= -1 & data$CLIP_fc >= 1] = "#FD8D62"
color.vec[data$RNA_fc <= -1 & data$CLIP_fc <= -1] = "#66C3A5"

v=ggplot(data,aes(x=RNA_fc,y=CLIP_fc))
v+geom_point(colour=color.vec,size=1)+labs(title="Volcanoplot",x="expression difference(log2FC)", y="methylation difference(log2FC)")+
  xlim(-5,5)+ylim(-10,10)+
  geom_vline(xintercept =0,linetype=4)+
  geom_hline(yintercept=0,linetype=4)+
  geom_text(label=data$label,colour="black",size=0.5,vjust=0.5)
ggsave("diff_vocanno.pdf")
