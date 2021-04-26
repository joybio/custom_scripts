#http://www.bioconductor.org/packages/release/bioc/vignettes/Guitar/inst/doc/Guitar-Overview.R
library(Guitar)
stBedFiles=list(system.file("extdata","ipcom_col_vs_ipcom_ww.hyper.bed.bed", package="Guitar"))
txdb=makeTxDbFromGFF("/home/l/backup1/refgenome/Arabidopsis/Arabidopsis_thaliana.TAIR10.42.gff3")
#Ara_transcript_txdb <- makeGuitarTxdb(txdb,txPrimaryOnly = FALSE)
GuitarPlot(txTxdb = txdb, 
	stBedFiles = stBedFiles,
	headOrtail = TRUE,
	txpromoterLength = 0,
	txtailLength = 0,
	enableCI = FALSE,
	mapFilterTranscript = TRUE,
	pltTxType =  c("mrna"),
	stGroupName = c("hyper"),
	miscOutFilePrefix = "consistant hyper in ipcom_col and ipcom_ww")
