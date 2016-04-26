args <- commandArgs(trailingOnly = TRUE)


#file <- read.delim("~/sunyd/identify/oryza_rnaseq/SRRpi/methylation/lncrna.methylation",header=F,stringsAsFactors=F)
file <- read.delim(args[1],header=T,row.names = NULL,stringsAsFactors=F)
file[which(file$strand == "."),9] = rep("unknown",length(file[which(file$strand == "."),9]))
write.table(file,args[2],quote=F,sep="\t",col.names=TRUE,row.names=FALSE)


