args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
#allcount.avs.high <-read.delim(file=args[1],row.names=1,header=F)
#allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
#allgroup <- allgroup[,2]
filename <- args[4]

#########read in the gene names and lncrna names for analysis
gene.high <- read.delim(file=args[1],header=T,stringsAsFactors = FALSE) 
genenamepure.high<- as.matrix(rownames(gene.high))
lnc.high <- read.delim(file=args[2],header=T,stringsAsFactors = FALSE) 
lncnamepure.high<- as.matrix(rownames(lnc.high))
lnc.high.del <- read.delim(file=args[3],header=T,stringsAsFactors = FALSE) 
lncnamepure.high.del <- as.matrix(rownames(lnc.high.del))
log2allcount.avs.high <- rbind(gene.high,lnc.high,lnc.high.del)
legendname <- c("Genes","LncRNAs","Intron","Inter","Exon","Anti","iso","match")
color8 <- c("black","red","blue","green","purple","pink","yellow","grey")
legendname3 <- c("Genes","HC-lncRNAs","Other_lncRNAs")
color3 <- c("blue","red","green")

##########built a DGElist and store your data in it#################
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(TCC)
library(cluster)

############calculate maximum cpm value for gene expression plot analysis###########
expgenehigh <- apply(gene.high,1,max)
explnchigh <- apply(lnc.high,1,max)
explnchigh.del <- apply(lnc.high.del,1,max)
valu <- c(data.matrix(expgenehigh),data.matrix(explnchigh),data.matrix(explnchigh.del))
tim <- c(length(expgenehigh),length(explnchigh),length(explnchigh.del))   
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
png(paste(filename,"_expression.png",sep=""))
boxplot(values ~ vars, data = df, col = color3,main="Transcript expression distribution (Maximum)",ylab="log2FPKM Maximum")
dev.off()

#plot(density(explnchigh),type = "l", xlab = "log2FPKM",ylab = "Density",col="red",main="Expression level distribution (Maximum)")
#lines(density(expgenehigh),col = "black")
#lines(density(iexplnchigh),col = "blue")
#lines(density(uexplnchigh),col = "green")
#lines(density(oexplnchigh),col = "purple")
#lines(density(xexplnchigh),col = "pink")
#lines(density(jexplnchigh),col = "yellow")
#lines(density(eexplnchigh),col = "grey")
#legend("topright",legend=legendname,col=color8,bg="white",lwd=2)

####################the PCA analysis#################
#keep.all <- apply(allcount,1,function(x){any(x > 1)})
#allcounts.high <- allcount[keep.all,]
#tcounts <- t(allcounts.high)
#pca.total <- prcomp(log2(tcounts+1), retx=TRUE)
#png(paste("normhtseq","pca",".png",sep=""))
#	c <- round(100*summary(pca.total)$importance[2,1],digits=2)
#	d <- round(100*summary(pca.total)$importance[2,2],digits=2)
#	plot(pca.total$x[,1:2], pch=c(1:25), col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),main="PCA Plot of Samples")
#	legend("topleft",cex=0.5,border=F, c(row.names(pca.total$x)),pch=c(1:25),col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)),bty="n")
#dev.off()

############calculate tissue specific score###########3
ts <- ROKU(log2allcount.avs.high, upper.limit = 0.25, sort = FALSE)
tsorigin <- as.matrix(ts$H)
#write.table(tsorigin,paste(filename,"_tissuespecific.txt",sep=""),quote=F,sep="\t",col.names=FALSE,row.names=TRUE)
tsgenehigh <- as.matrix(tsorigin[which(rownames(tsorigin) %in% genenamepure.high),])
tslnchigh <- as.matrix(tsorigin[which(rownames(tsorigin) %in% lncnamepure.high),])       
tslnchigh.del <- as.matrix(tsorigin[which(rownames(tsorigin) %in% lncnamepure.high.del),])       

valu <- c(tsgenehigh,tslnchigh,tslnchigh.del)
tim <- c(nrow(tsgenehigh),nrow(tslnchigh),nrow(tslnchigh.del))
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
png(paste(filename,"_tissuespecific.png",sep=""))
boxplot(values ~ vars, data = df, col = color3,main="Transcript tissue specific score distribution",ylab="Tissue specific score")
dev.off()

###########seperate the data in Highlight and lowlight condition###########
#a = seq(1,ncol(lnc.high),by=2)
#b = seq(2,ncol(lnc.high),by=2)
#lr <- lnc.high[,b]
#keep.lr <- apply(lr,1,function(x){any(x > 1)})    #####only extract FPKM > 1 in at least one sample
#lr.high <- lr[keep.lr,]
#hr <- lnc.high[,a]
#keep.hr <- apply(hr,1,function(x){any(x > 1)})
#hr.high <- hr[keep.hr,]

###identigy kmeans value#######
kmidentify <- function(x){
wss <- (nrow(x)-1)*sum(apply(x,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(x,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
	}
kmidentify(lnc.high)
############perform Kmeans, Hcluster heatmap###############
heatmapfunc <- function(counts,num,filename){
	####kmeans clustering####
	set.seed(1)
	heat <- pheatmap(counts,kmeans_k=num,scale="row",color=colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 6)
	#heat <- pheatmap(lnclog2counts,scale="row",color=colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 4)
	heat1 <- ((as.matrix(counts[names((sort(heat$kmeans$cluster))),])))
	write.table(cbind(sort(heat$kmeans$cluster),counts[names((sort(heat$kmeans$cluster))),]),file=paste(filename,"kcluster.txt",sep=""),quote=F,sep="\t" )
	print(table(heat$kmeans$cluster))         ####check the gene number of each cluster
	png(paste(filename,"kheatmap.png",sep=""))
	col1 <- colorRampPalette(c("green","black","red"))(100)  ##set the color for heatmap
	pheatmap(t(heat1),scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 18)
	dev.off()
	####Hierachical clustering#############
#	hclus <- agnes(t(counts),method = "gaverage")
#	png(paste(filename,"hheatmap.png",sep=""))
#	pheatmap(hclus$height,scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 8)
#	dev.off()
}
#heatmapfunc(log2allcount.avs.high,"lnc_gene")
#gene.high <- gene.high[,c(4,7,3,5,9,8,1,6,2)]
heatmapfunc(gene.high,6,paste(filename,"gene",sep=""))
#lnc.high <- lnc.high[,c(4,7,3,5,9,8,1,6,2)]
heatmapfunc(lnc.high,6,filename)

