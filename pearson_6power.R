args <- commandArgs(trailingOnly = TRUE)


#file <- read.delim("~/sunyd/identify/oryza_rnaseq/SRRpi/methylation/lncrna.methylation",header=F,stringsAsFactors=F)
file <- read.delim(args[1],header=F,row.names = NULL,stringsAsFactors=F)
options("scipen"=100)
file[,3] <- round(file[,3]^(1/6),4)

write.table(file,args[2],quote=F,sep="\t",col.names=FALSE,row.names=FALSE)


