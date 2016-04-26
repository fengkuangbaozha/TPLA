all.count <-read.table(file="ar63.count1.txt",row.names=1,header=T)

group <- c(rep("h0ck", 4), rep("h0", 4),rep("h1.5ck", 2), rep("h1.5", 3),rep("h2.5ck",2),rep("h2.5",3),rep("h3ck",2),rep("h3",3),rep("pi",3),rep("m45ck",4),rep("m45",4),rep("h2ck",3),rep("h2",3))
library(edgeR)
cds <- DGEList(all.count,group =group)
cds <- calcNormFactors(cds)
keep <- rowSums(cpm(cds)>1)>=2 
cds <-cds[keep,]


log2counts <-log2(cds$counts+1)


pdf(file="DM.pairplot.pdf")
#png(file="ar63.pairplot.png")
panel.cor <- function(x,y, ...)
{ 
  par(usr=c(0,1,0,1))
   txt <- as.character(format(cor(x,y),digits=4))
   text(0.5,0.5,txt,cex =4*abs(cor(x,y)))
}
pairs(log2counts[,1:40],upper.panel=panel.cor,,main="Relationship of 40 samples in ar63 RNA-Seq")


dev.off()

library(edgeR)
group[c(1:8)] ->group0h
group[c(9:13)] ->group1.5h
group[c(14:18)] ->group2.5h
group[c(19:23)] ->group3h
group[c(27:34)] ->group45
group[c(35:39)] ->group2h


all.count[,c(1:8)] ->h0
all.count[,c(9:13)] ->h1.5
all.count[,c(14:18)] ->h2.5
all.count[,c(19:23)] ->h3
all.count[,c(27:34)] ->h45
all.count[,c(35:39)] ->h2


####0h
DE <- function(x,y,a,b){
cds <- DGEList(x,group =y)
cds <- calcNormFactors(cds)
keep <- rowSums(cpm(cds)>1)>=a 
cds <-cds[keep,]
cds <- estimateCommonDisp(cds,verbose=TRUE)
cds <- estimateTagwiseDisp(cds)
cpmcounts <- cpm(cds$counts) 
log2counts <-log2(cds$counts+1)

#################the spearman correlation analysis###############
#pdf(file="SPL.pairplot.pdf")
pdf(file=paste(y[5],".pairplot.pdf",sep=""))
panel.cor <- function(x,y, ...)
{ 
  par(usr=c(0,1,0,1))
   txt <- as.character(format(cor(x,y),digits=4))
   text(0.5,0.5,txt,cex =2*abs(cor(x,y)))
}
pairs(log2counts[,1:(a+b)],upper.panel=panel.cor,main=paste("Relationship of  samples in " ,y[5], " RNA-Seq",sep=""))


dev.off()


####################the PCA analysis#################
tcounts <- t(cpmcounts)

pca.total <- prcomp(log2(tcounts+0.00000001), retx=TRUE)
pdf(paste("pca.",y[5],".pdf",sep=""))
round(100*summary(pca.total)$importance[2,1],digits=2) -> c
round(100*summary(pca.total)$importance[2,2],digits=2) -> d
plot(pca.total$x[,1:2], pch=c(15:18), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),col=c(rep("black",a),rep("red",b)),main="PCA Plot of Samples")

legend("bottomright",cex=0.6,border=F, c(row.names(pca.total$x)),pch=c(15:18), col=c(rep("black",a),rep("red",b)),bty="n")
dev.off()
#################get the top differentially expressed genes################
et.RCSC <-exactTest(cds,pair=c(y[1],y[5]))
rawcount <-as.matrix(x)[(rownames(subset(topTags(et.RCSC,n=30000)$table,FDR<0.05&abs(logFC) >=  1))),]
rawcount2 <-as.matrix(x)[(rownames(topTags(et.RCSC,n=30000)$table,)),]
write.table(file=paste(y[5],".fdr5.txt",sep=""),cbind(subset(topTags(et.RCSC,n=30000)$table,FDR<0.05&abs(logFC) >=  1),rawcount),quote=F,sep="\t")
write.table(file=paste(y[5],".1.5fdr5.txt",sep=""),cbind(subset(topTags(et.RCSC,n=30000)$table,FDR<0.05&abs(logFC) >=  0.58),rawcount),quote=F,sep="\t")
write.table(file=paste(y[5],"_diff.txt",sep=""),cbind(topTags(et.RCSC,n=100000)$table,rawcount2),sep="\t",quote=F)
}
DE(h0,group0h,4,4)
DE(h45,group45,4,4)
DE(h1.5,group1.5h,2,3)
DE(h3,group3h,2,3)
DE(h2.5,group2.5h,2,3)




