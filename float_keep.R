args <- commandArgs(trailingOnly = TRUE)


#file <- read.delim("~/sunyd/identify/oryza_rnaseq/SRRpi/methylation/lncrna.methylation",header=F,stringsAsFactors=F)
file <- read.delim(args[1],header=T,row.names = 1,stringsAsFactors=F)
options("scipen"=100)
file1 = round(file,4)
write.table(file1,args[2],quote=F,sep="\t",col.names=NA,row.names=TRUE)


