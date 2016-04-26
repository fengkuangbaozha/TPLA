args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
#allcount.avs.high <-read.delim(file=args[1],row.names=1,header=F)
#allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
#allgroup <- allgroup[,2]
filename <- args[5]

allcount <-read.delim(file=args[1],row.names=1,header=T)    ###read in the cuffnorm gene count info
combgroup <- read.delim(file=args[2],header=F)     ###combine the same group or tissue together info
allcount1 = allcount
colnames(allcount1) <- combgroup[,1]   ##put the group name for the colnames of the data

#########read in the gene names and lncrna names for analysis
genename <- read.delim(file=args[3],header=F,stringsAsFactors = FALSE) 
genenamepure <- unique(as.matrix(genename[,2]))
lncname <- read.delim(file = args[4],header=F)
lncnamepure <- unique(as.matrix(lncname[,2]))
#lncnamedel <- read.delim(file = args[5],header=F,stringsAsFactors = FALSE)
#lncnamepure.del <- unique(as.matrix(lncnamedel[,5]))

library(TCC)

############combine the replicate counts together and use mean value##################
allcount.avs = data.frame(t(apply(allcount1,1,tapply,names(allcount1),mean)))    ########calculate the group mean fpkm value
log2countsorigin <- log2(allcount.avs+1)
write.table(log2countsorigin,paste(filename,"_expression_mean.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)
#explnccpm <- t(apply(log2countsorigin[which(rownames(log2countsorigin) %in% lncnamepure),],1,summary))   ######find fpkm max, min, medium, mean
#explnccpm <- explnccpm[,c(1,3,4,6)]
#write.table(explnccpm,paste(filename,"_expression_summary.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)

##########extract the fpkm value high gene name##############
genenamepure <- as.matrix(rownames(log2countsorigin)[which(rownames(log2countsorigin) %in% genenamepure)])   #####extract the high expressed gene name
lncnamepure <- as.matrix(rownames(log2countsorigin)[which(rownames(log2countsorigin) %in% lncnamepure)])
#lncnamepure.del <- as.matrix(rownames(log2countsorigin)[which(rownames(log2countsorigin) %in% lncnamepure.del)])
gene <- log2countsorigin[which(rownames(log2countsorigin) %in% genenamepure),]  ######extract the high expressed genename expression, log 2 value
write.table(gene,paste(filename,"_expression_gene.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)
lnc <- log2countsorigin[which(rownames(log2countsorigin) %in% lncnamepure),]
write.table(lnc,paste(filename,"_expression_lnc.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)
#lnc.del <- log2countsorigin[which(rownames(log2countsorigin) %in% lncnamepure.del),]
#write.table(lnc.del,paste(filename,"_expression_lnc.del.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)

#ts <- ROKU(log2countsorigin, upper.limit = 0.25, sort = FALSE)                                                                                                          
#tsorigin <- as.matrix(ts$H)
#write.table(tsorigin,paste(filename,"_tissuespecific.txt",sep=""),quote=F,sep="\t",col.names=FALSE,row.names=TRUE)

##########extract the fpkm value high gene name##############
keep <- apply(allcount.avs,1,function(x){any(x > 1)})    #####only extract FPKM > 1 in at least one sample
allcount.avs.high <- allcount.avs[keep,]
log2allcount.avs.high <- log2(allcount.avs.high+1)
genenamepure.high <- as.matrix(rownames(log2allcount.avs.high)[which(rownames(log2allcount.avs.high) %in% genenamepure)])   #####extract the high expressed gene name
lncnamepure.high <- as.matrix(rownames(log2allcount.avs.high)[which(rownames(log2allcount.avs.high) %in% lncnamepure)])
#lncnamepure.high.del <- as.matrix(rownames(log2allcount.avs.high)[which(rownames(log2allcount.avs.high) %in% lncnamepure.del)])
gene.high <- log2allcount.avs.high[which(rownames(log2allcount.avs.high) %in% genenamepure),]  ######extract the high expressed genename expression, log 2 value
write.table(gene.high,paste(filename,"_expression_gene.high.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)
lnc.high <- log2allcount.avs.high[which(rownames(log2allcount.avs.high) %in% lncnamepure),]
write.table(lnc.high,paste(filename,"_expression_lnc.high.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)
#lnc.high.del <- log2allcount.avs.high[which(rownames(log2allcount.avs.high) %in% lncnamepure.del),]
#write.table(lnc.high.del,paste(filename,"_expression_lnc.high.del.txt",sep=""),quote=F,sep="\t",col.names=NA,row.names=TRUE)
