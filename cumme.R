args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
##############this program should be finished under R graph because read in needs too much time############
library("cummeRbund")
setwd(dir)
cuff <- readCufflinks()
print(samples(isoforms(cuff)))
filename <- args[2] 
c1 <- args[3]
c2 <- args[4]
pdf(paste(filename,".isodiff.pdf",sep=""))
gene.fpkm <- repFpkmMatrix(isoforms(cuff))
fpkm.pca <- prcomp(t(gene.fpkm), retx=TRUE)
c <- round(100*summary(fpkm.pca)$importance[2,1],digits=2)
d <- round(100*summary(fpkm.pca)$importance[2,2],digits=2)
plot(fpkm.pca$x[,1:2], pch=c(15:20), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),col=rainbow(18),main="FPKM PCA Plot of Two Conditions")
legend("bottomright",cex=0.6,border=F, c(row.names(fpkm.pca$x)),pch=c(15:20), col=rainbow(18),bty="n")
csDensity(isoforms(cuff),replicates=T)
csDendro(isoforms(cuff),replicates=T)
dev.off()
sigids<-getSig(cuff,c1,c2,alpha=0.05,level="isoforms") 
myGenes<-getGenes(cuff,sigids,sampleIdList=c(c1,c2))
iso.diff <- diffData(isoforms(myGenes))
iso.diff <- iso.diff[iso.diff$significant == "yes",]
iso.diff <- iso.diff[abs(iso.diff$log2_fold_change) >= 0.6,]
print(dim(iso.diff))
write.table(iso.diff, paste(filename,".isodiff.txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE,quote=FALSE) 


