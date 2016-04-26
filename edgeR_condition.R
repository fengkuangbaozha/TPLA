args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
allcount <-read.delim(file=args[1],row.names=1,header=F)
allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
allgroup <- allgroup[,2]
pair1 <- args[3]
pair2 <- args[4]
filename <- args[5]
##########built a DGElist and store your data in it#################
library(edgeR)
DE <- function(all.count,group,num,name){
cds <- DGEList(all.count,group = group)
cds <- calcNormFactors(cds)
keep <- rowSums(cpm(cds)>1)>=num   ####value > 1 and at least $num numbers
cds <-cds[keep,]      
dim(cds)

pdf(paste(name,"pdf",sep="."))
plotMDS(cds, labels=group)  ##### check the replicate distribution
dev.off
############estimate the common distribution of the data###########
cds <- estimateCommonDisp(cds,verbose=TRUE)
#names(cds)
#sqrt(cds$common.dispersion)   ##the common distribution value
cds <- estimateTagwiseDisp(cds)
plotBCV(cds)
abline(h=sqrt(cds$common.dispersion), col="firebrick", lwd=3)
##############get the counts and get the differential expressed genes###########
cpmcounts <- cpm(cds$counts) ###get the count number
log2counts <-log2(cds$counts+1)
#et.RCSC <-exactTest(cds,pair=c(group[1],group[length(group)]))
et.RCSC <-exactTest(cds,pair=c(1,2))
toptag <-as.matrix(all.count)[(rownames(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1))),]
write.table(file=paste(name,"_diff.txt",sep=""),cbind(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1),toptag),quote=F,sep="\t")
return(toptag)
}
count1 <- allcount[,c(which(allgroup==pair1))]     ####group 1 remains to be comparied
count2 <- allcount[,c(which(allgroup==pair2))]
group1 <- allgroup[c(which(allgroup==pair1))]
group2 <- allgroup[c(which(allgroup==pair2))]

diff <- DE(cbind(count1,count2),c(group1,group2),3,filename)
