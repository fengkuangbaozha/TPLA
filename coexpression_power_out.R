args <- commandArgs(trailingOnly = TRUE)
#########read in the gene names and lncrna names for analysis
#gene.high <- read.delim(file="~/sunyd/identify/tomato_rnaseq/SRRm82/heatmap/normlnc_expression_gene.high.txt",header=T,stringsAsFactors = FALSE)
#lnc.high <- read.delim(file="~/sunyd/identify/tomato_rnaseq/SRRm82/heatmap/normlnc_expression_lnc.high.txt",header=T,stringsAsFactors = FALSE)
#lnc.high.del <- read.delim(file="~/sunyd/identify/tomato_rnaseq/SRRm82/heatmap/normlnc_expression_lnc.high.del.txt",header=T,stringsAsFactors = FALSE)
#genename <- read.delim(file="~/sunyd/identify/tomato_rnaseq/SRRm82/cuff85.combined.gtf-correspond-genename-GO",header=F,stringsAsFactors = FALSE)
#kmean.floral <- read.delim(file="~/sunyd/identify/tomato_rnaseq/SRRm82/heatmap/normlnckcluster.floral",header=F,stringsAsFactors = FALSE)
#filename <- "~/sunyd/identify/tomato_rnaseq/SRRm82/coexpression/floral."
#number <- as.numeric(16)
gene.high <- read.delim(file=args[1],header=T,stringsAsFactors = FALSE) 
genenamepure.high<- as.matrix(rownames(gene.high))
lnc.high <- read.delim(file=args[2],header=T,stringsAsFactors = FALSE) 
lncnamepure.high<- as.matrix(rownames(lnc.high))
lnc.high.del <- read.delim(file=args[3],header=T,stringsAsFactors = FALSE) 
lncnamepure.high.del <- as.matrix(rownames(lnc.high.del))
genename <- read.delim(file=args[4],header=F,stringsAsFactors = FALSE) 
kmean.floral <- read.delim(file=args[5],header=F,stringsAsFactors = FALSE)
kmean.floral.log2counts <- rbind(lnc.high[as.matrix(kmean.floral[,1]),],gene.high)
log2allcount.avs.high <- rbind(gene.high,lnc.high,lnc.high.del)
filename <- args[6]
number <- as.numeric(args[7])
##################start coexpression analysis############
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(32)
library(flashClust)

#########calculate the softpower number for doing TOM similarity##########
coexpression.power <- function(datExpr,filename){
########## Choosing the soft-thresholding power: analysis of network topology###############
gene.names=rownames(t(datExpr))
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
		##### Plot the results
#sizeGrWindow(9, 5)
pdf(paste(filename,"softthreshold.pdf",sep=""))
par(mfrow = c(1,2));
cex1 = 0.9;
		# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
		# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")
		# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
}

kmean.floral <- coexpression.power(t(kmean.floral.log2counts),filename)

############get the correlation value between lncrna and gene and find the top 20 closest gene for lncrna###########
lnc_gene <- rbind(lnc.high,gene.high)      ######lnc.high.del
adj.test = adjacency(t(lnc_gene),selectCols = c(1:nrow(lnc.high)),type = "unsigned", power = number)  #########caculate the correlation value between lncrna and gene
adj.test_lnc <- t(adj.test[as.matrix(rownames(gene.high)),])     #########only extract lncRNA_gene correlation
adj.test_lnc <- adj.test_lnc[,which(as.matrix(colnames(adj.test_lnc)) %in% as.matrix(genename[,2]))]
cor.genename <- genename[which(as.matrix(genename[,2]) %in% as.matrix(colnames(adj.test_lnc))),c(2,5)]     #########correspond TCONS and Csa
cor.genename <- cor.genename[order(cor.genename[,1]),]               ###########order the name
colnames(adj.test_lnc) <- cor.genename[,2]           ######give the cucumber gff name to the correlation matrix
adj.test_lnc.top <- t(apply(adj.test_lnc,1,function(x){y <- colnames(adj.test_lnc)[order(x,decreasing = TRUE)]; y <- y[1:20]; return(y)})) ###get the top 20 large correlation value
write.table(adj.test_lnc.top, paste(filename,"correlation.txt",sep=""), sep="\t", row.names=TRUE, col.names=FALSE,quote=FALSE)

#################find the correlation over 0.95 closest gene for lncrna#########
system.time(adj.list <- lapply(seq_len(nrow(adj.test_lnc)), function(i) adj.test_lnc[i,]))   ###convert matrix to list by row
#system.time(li <- split(adj.test_lnc,row(adj.test_lnc)))
names(adj.list) <- rownames(adj.test_lnc)
adj.keep <- lapply(adj.list,function(x) x[x>= 0.95^number])
adj.keep.convert <- sapply(names(adj.keep),function(x) paste(x,paste(names(adj.keep[[x]]),collapse="\t"),sep="\t"))    #####write to the file
write(adj.keep.convert,file=paste(filename,"correlation0.95.txt",sep=""),append=F)
