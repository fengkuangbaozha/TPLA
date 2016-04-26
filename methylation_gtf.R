args <- commandArgs(trailingOnly = TRUE)
print("gtf starting")

file <- read.delim(file=args[1],header=T,stringsAsFactors=F)
options("scipen"=100)
a = file
#a = file[which(file$ratio!=0),]   ####only keep ratio is not 0
b = cbind(a$chr,rep("BSMAP",nrow(a)),rep("exon",nrow(a)),a$pos,a$pos,c("."),a$strand,c("."),paste(paste(paste(a$context,rownames(a),sep=""),a$C_count,sep="; "),a$CT_count,sep="; "))
write.table(b,args[2],quote=F,sep="\t",col.names=FALSE,row.names=FALSE)
