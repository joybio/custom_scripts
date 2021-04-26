library(Vennerable)
pdf("m6A_venn.CJZ.KYG.pdf")
venn <- Venn(SetNames=c("m6A-seal-CJZ","m6A-seal-KYG"),Weight=c(0,5252,4785,1669))
p <- compute.Venn(venn)
gp <- VennThemes(p)
gp[["Face"]][["11"]]$fill <-"#66C3A5"
gp[["Face"]][["01"]]$fill <-"#8DA1CB"
gp[["Face"]][["10"]]$fill <-"#FD8D62"

gp[["SetText"]][["Set1"]]$col <- "black"
gp[["SetText"]][["Set2"]]$col <- "black"

gp[["SetText"]][["Set1"]]$fontsize <- 10
gp[["SetText"]][["Set2"]]$fontsize <- 10

gp[["Set"]][["Set1"]]$col <- "NA"
gp[["Set"]][["Set2"]]$col <- "NA"
gp[["Set"]][["Set3"]]$col <- "NA"

gp[["FaceText"]][["DarkMatter"]]$fontsize <- 0
gp[["FaceText"]][["11"]]$fontsize <- 10
gp[["FaceText"]][["10"]]$fontsize <- 10
gp[["FaceText"]][["01"]]$fontsize <- 10

plot(venn,show = list(FaceText = "weight", SetLabels = T, Faces = T,Universe=F,DarkMatter = T),doWeight=T,gp=gp)
dev.off()

pdf("m6A_venn.CJZ.SJZ.pdf")
venn <- Venn(SetNames=c("m6A-seal-CJZ","m6A-seal-SJZ"),Weight=c(0,5099,6881,1822))
plot(venn,show = list(FaceText = "weight", SetLabels = T, Faces = T,Universe=F,DarkMatter = T),doWeight=T,gp=gp)
dev.off()

pdf("m6A_venn.KYG.SJZ.pdf")
venn <- Venn(SetNames=c("m6A-seal-KYG","m6A-seal-SJZ"),Weight=c(0,4736,6985,1718))
plot(venn,show = list(FaceText = "weight", SetLabels = T, Faces = T,Universe=F,DarkMatter = T),doWeight=T,gp=gp)
dev.off()

pdf("m6A_venn.CJZ.KYG.SJZ.pdf")
venn <- Venn(SetNames=c("m6A-seal-CJZ","m6A-seal-KYG","m6A-seal-SJZ"),Weight=c(0,4074,3558,5807,1025,1178,1074,644))
p <- compute.Venn(venn)
gp <- VennThemes(p)

gp[["SetText"]][["Set1"]]$col <- "black"
gp[["SetText"]][["Set2"]]$col <- "black"
gp[["SetText"]][["Set3"]]$col <- "black"

gp[["SetText"]][["Set1"]]$fontsize <- 10
gp[["SetText"]][["Set2"]]$fontsize <- 10
gp[["SetText"]][["Set3"]]$fontsize <- 10

gp[["Set"]][["Set1"]]$col <- "NA"
gp[["Set"]][["Set2"]]$col <- "NA"
gp[["Set"]][["Set3"]]$col <- "NA"

gp[["FaceText"]][["DarkMatter"]]$fontsize <- 0
gp[["FaceText"]][["011"]]$fontsize <- 10
gp[["FaceText"]][["010"]]$fontsize <- 10
gp[["FaceText"]][["001"]]$fontsize <- 10
gp[["FaceText"]][["101"]]$fontsize <- 10
gp[["FaceText"]][["100"]]$fontsize <- 10
gp[["FaceText"]][["111"]]$fontsize <- 10
gp[["FaceText"]][["110"]]$fontsize <- 10

plot(venn,show = list(FaceText = "weight", SetLabels = T, Faces = T,Universe=F,DarkMatter = T),doWeight=T,gp=gp)
dev.off()
