#args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
#cpmcounts <-read.delim(file=args[1],row.names=1,header=F)
#allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
#allgroup <- allgroup[,2]
#filename <- args[3]

allcount <-read.delim("~/sunyd/identify/cucumber_rnaseq/htseq.result",row.names=1,header=F)    ###read in the htseq count info
allgroup <- read.delim("~/sunyd/identify/cucumber_rnaseq/htseqgroup",header=F)     ####the group info, second col is group 
allgroup <- allgroup[,2]        ##only extract the group info
combgroup <- read.delim("~/sunyd/identify/cucumber_rnaseq/htseqgroupcombine",header=F)     ###combine the same group or tissue together info
name <- combgroup[,3]   ##only extract the group info
lncnameorigin <- rownames(allcount)[grep("XLOC",rownames(allcount))]  ###extract only the lncRNA name position
genenameorigin <- rownames(allcount)[grep("Csa",rownames(allcount))]  ###extract only the lncRNA name position
lncgtf <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.gtf",header=F)   ##lncRNA name and classification
#name <- c(allgroup[1],paste(allgroup[2:4],1:3,sep=""),paste(allgroup[5:7],1:3,sep=""),paste(allgroup[8:20],1:13,sep=""),paste(allgroup[21:22],1:2,sep=""),paste(allgroup[23:24],1:2,sep=""),paste(allgroup[25:34],1:10,sep=""),paste(allgroup[35],1,sep=""),paste(allgroup[36:37],1:2,sep=""),paste(allgroup[38],1,sep=""),paste(allgroup[39],1,sep=""),paste(allgroup[40:42],1:3,sep=""),paste(allgroup[43:75],1:33,sep=""))
#name <- c(allgroup[1],paste(allgroup[2:4],1:3,sep=""),paste(allgroup[5:7],1:3,sep=""),paste(allgroup[8:20],1:13,sep=""),paste(allgroup[21:28],1:8,sep=""),paste(allgroup[29],1,sep=""),paste(allgroup[30:31],1:2,sep=""),paste(allgroup[32],1,sep=""),paste(allgroup[33],1,sep=""),paste(allgroup[34:36],1:3,sep=""),paste(allgroup[37:69],1:33,sep=""))

############extract the lncRNA name info and classify them by position############
nameext <- function(x){
 y <- strsplit(as.character(x[,9]),";") #####split the last item, and extract geneid, transcriptid
 z <- data.frame(NA,NA,NA,NA)
 for (i in 1:length(y)){
         a <- as.matrix(y[[i]])
         z[i,1] <- a[apply(a, 1, function (x) grepl('gene_id', x))]
         z[i,2] <- a[apply(a, 1, function (x) grepl('transcript_id', x))]
         z[i,3] <- a[apply(a, 1, function (x) grepl('class_code', x))]
         z[i,4] <- a[apply(a, 1, function (x) grepl('oId', x))]
                }
 z <- cbind(z,(x[,5]-x[,4]),x[,c(4,5,7)])   ###add the exon position and length info####
 colnames(z) <- c("gene_id","transcript_id","class_code","oId","exonlength","start","end","strand")
 return(z)
}
lncgenename <- nameext(lncgtf)  ######calculate the lncRNA info, extract name info
codeext <- function(x,y){
    class <- x[which(x[,3]== (paste(" class_code ",y,sep=""))),]
	class1 <- strsplit(as.character(class[,1])," ")        
	class2 <- strsplit(as.character(class[,2])," ")
	z <- data.frame(NA,NA)
	 for (i in 1:length(class1)){	
		a <- class1[[i]]
		b <- class2[[i]]
		z[i,1] <- a[grepl('XLOC', a)]
		z[i,2] <- b[grepl('TCONS', b)]
		}
	colnames(z) <- c("gene_id","transcript_id")
	return(z)   
        }
jclass <- unique(codeext(lncgenename,"j")[,1,drop=FALSE])
iclass <- unique(codeext(lncgenename,"i")[,1,drop=FALSE])
oclass <- unique(codeext(lncgenename,"o")[,1,drop=FALSE])
uclass <- unique(codeext(lncgenename,"u")[,1,drop=FALSE])
xclass <- unique(codeext(lncgenename,"x")[,1,drop=FALSE])
eclass <- unique(codeext(lncgenename,"=")[,1,drop=FALSE])
cclass <- unique(codeext(lncgenename,"c")[,1,drop=FALSE])
##########built a DGElist and store your data in it#################
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(TCC)
library(cluster)

cds <- DGEList(allcount,group = allgroup)
cds <- calcNormFactors(cds)
##############get the counts and log2counts for all lncRNA and genes ###########
cpmorigin <- cpm(cds$counts)
log2countsorigin <- log2(cpmorigin+1)
##############get the counts and log2counts for high expressed lncRNA and genes ########### 
keep <- rowSums(cpm(cds)>1)>=1   ####value > 1 and at least $num numbers, high expression value ones
cds <-cds[keep,]      
print(dim(cds))
cpmcountshigh <- cpm(cds$counts) ###get the count number
log2countshigh <- log2(cpmcountshigh+1)
lncnamehigh <- rownames(cpmcountshigh)[grep("XLOC",rownames(cpmcountshigh))]   #####only extract lncRNA name rows and delete gene name rows
genenamehigh <- rownames(cpmcountshigh)[grep("Csa",rownames(cpmcountshigh))]   #####only extract lncRNA name rows and delete gene name rows
############calculate maximum cpm value for gene expression plot analysis###########
expgeneorigin <- apply(log2countsorigin[genenameorigin,],1,max)
explncorigin <- apply(log2countsorigin[lncnameorigin,],1,max)
expgenehigh <- apply(log2countshigh[genenamehigh,],1,max)
explnchigh <- apply(log2countshigh[lncnamehigh,],1,max)

pdf("lnc_expression.pdf")
plot(density(explnchigh),type = "l", xlab = "cpm",ylab = "Density",col="red",main="Expression level distribution (Maximum)")
lines(density(expgenehigh),col = "blue")
legend("topright",legend=c("Long noncoding RNAs","Genes"),col=c("red","blue"),bg="white",lwd=2)
dev.off()
############calculate tissue specific score###########3
ts <- ROKU(log2countsorigin, upper.limit = 0.25, sort = FALSE)
tsgeneorigin <- as.matrix(as.matrix(ts$H)[genenameorigin,])           ####ROKU(log2countsorigin[-lncnameorigin,], upper.limit = 0.25, sort = TRUE)
tslncorigin <- as.matrix(as.matrix(ts$H)[lncnameorigin,])           ####ROKU(log2countsorigin[lncnameorigin,], upper.limit = 0.25, sort = TRUE)
tsgenehigh <- as.matrix(tsgeneorigin[genenamehigh,])
tslnchigh <- as.matrix(tslncorigin[lncnamehigh,])                     ####ROKU(log2countshigh[lncnamehigh,], upper.limit = 0.25, sort = TRUE)

pdf("lnc_tissuespecific.pdf")
df <- data.frame(values = c(tslnchigh,tsgenehigh),vars = rep(c("Long noncoding RNAs","Genes"), times = c(nrow(tslnchigh),nrow(tsgenehigh))))
boxplot(values ~ vars, data = df)
dev.off()

####################the PCA analysis#################
tcounts <- t(cpmcountshigh)
pca.total <- prcomp(log2(tcounts+1), retx=TRUE)
pdf(paste("htseq","pca",".pdf",sep=""))
	c <- round(100*summary(pca.total)$importance[2,1],digits=2)
	d <- round(100*summary(pca.total)$importance[2,2],digits=2)
	plot(pca.total$x[,1:2], pch=c(1:25), col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),main="PCA Plot of Samples")
	legend("topleft",cex=0.5,border=F, c(row.names(pca.total$x)),pch=c(1:25),col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)),bty="n")
dev.off()

############combine the replicate counts together and use mean value##################
cpmcounts <- matrix(NA, nrow=nrow(cpmcountshigh),ncol=28)   ###combine the replicate together
rownames(cpmcounts) <- rownames(cpmcountshigh)
cpmcounts[,1] <- apply(cpmcountshigh[,1:2],1,mean)
cpmcounts[,2] <- apply(cpmcountshigh[,3:4],1,mean)
cpmcounts[,3] <- apply(cpmcountshigh[,5:6],1,mean)
cpmcounts[,4] <- apply(cpmcountshigh[,7:9],1,mean)
cpmcounts[,5] <- apply(cpmcountshigh[,10:12],1,mean)
cpmcounts[,6] <- cpmcountshigh[,13]
cpmcounts[,7] <- cpmcountshigh[,14]
cpmcounts[,8] <- apply(cpmcountshigh[,15:16],1,mean) 
cpmcounts[,10] <- apply(cpmcountshigh[,17:29],1,mean)
cpmcounts[,28] <- apply(cpmcountshigh[,30:32],1,mean)
cpmcounts[,11] <- apply(cpmcountshigh[,33:34],1,mean)
cpmcounts[,12] <- apply(cpmcountshigh[,35:36],1,mean)
cpmcounts[,13] <- apply(cpmcountshigh[,37:38],1,mean)
cpmcounts[,14] <- apply(cpmcountshigh[,39:40],1,mean)
cpmcounts[,16] <- apply(cpmcountshigh[,41:43],1,mean)
cpmcounts[,17] <- apply(cpmcountshigh[,44:46],1,mean)
cpmcounts[,19] <- apply(cpmcountshigh[,47:49],1,mean)
cpmcounts[,18] <- apply(cpmcountshigh[,50:52],1,mean)
cpmcounts[,21] <- apply(cpmcountshigh[,53:55],1,mean)
cpmcounts[,20] <- apply(cpmcountshigh[,56:58],1,mean)
cpmcounts[,23] <- apply(cpmcountshigh[,59:61],1,mean)
cpmcounts[,24] <- apply(cpmcountshigh[,62:64],1,mean)
cpmcounts[,22] <- apply(cpmcountshigh[,65:67],1,mean)
cpmcounts[,25] <- apply(cpmcountshigh[,68:70],1,mean)
cpmcounts[,26] <- apply(cpmcountshigh[,71:73],1,mean)
cpmcounts[,27] <- cpmcountshigh[,74]
cpmcounts[,15] <- cpmcountshigh[,75]
#cpmcounts <- as.data.frame(cpmcounts[,-1:-3]) #####delete gua,ba,bing
cpmcounts <- as.data.frame(cpmcounts[,-9])    #######combined cpmcounts for all tissue samples
colnames(cpmcounts) <- name
log2cpmcounts <- log2(cpmcounts+1)
#######get the lncRNA and gene cpm and log2 counts for heatmap analysis#############
lnccpmcounts <- cpmcounts[lncnamehigh,]                
lnclog2counts <- log2(lnccpmcounts+1)     ###lncRNA log2counts
print(head(lnclog2counts))
genecpmcounts <- cpmcounts[genenamehigh,]  
genelog2counts <- log2(genecpmcounts+1)    ####gene log2 counts


###identigy kmeans value#######
kmidentify <- function(x){
wss <- (nrow(x)-1)*sum(apply(x,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(x,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
	}
############perform Kmeans, Hcluster heatmap###############
heatmapfunc <- function(counts,filename){
	####kmeans clustering####
	set.seed(1)
	heat <- pheatmap(counts,kmeans_k=6,scale="row",color=colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 6)
	#heat <- pheatmap(lnclog2counts,scale="row",color=colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 4)
	heat1 <- ((as.matrix(counts[names((sort(heat$kmeans$cluster))),])))
	write.table(cbind(sort(heat$kmeans$cluster),counts[names((sort(heat$kmeans$cluster))),]),file=paste(filename,"kcluster.txt",sep=""),quote=F,sep="\t" )
	print(table(heat$kmeans$cluster))         ####check the gene number of each cluster
	pdf(paste(filename,"kheatmap.pdf",sep=""))
	col1 <- colorRampPalette(c("green","black","red"))(100)  ##set the color for heatmap
	pheatmap(t(heat1),scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 8)
	dev.off()
	####Hierachical clustering#############
	hclus <- agnes(t(counts),method = "gaverage")
	pdf(paste(filename,"hheatmap.pdf",sep=""))
	pheatmap(hclus$height,scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 8)
	dev.off()
}
#heatmapfunc(log2cpmcounts,"lnc_gene")
heatmapfunc(genelog2counts,"gene")
heatmapfunc(lnclog2counts,"lnc")
heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% uclass$gene_id),],"lncu")
heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% iclass$gene_id),],"lnci")
heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% oclass$gene_id),],"lnco")
heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% jclass$gene_id),],"lncj")
heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% xclass$gene_id),],"lncx")
############estimate the common distribution of the data###########
#cds <- estimateCommonDisp(cds,verbose=TRUE)
	#names(cds)
	#sqrt(cds$common.dispersion)   ##the common distribution value
#cds <- estimateTagwiseDisp(cds)
	#plotBCV(cds)
	#abline(h=sqrt(cds$common.dispersion), col="firebrick", lwd=3)

#et.RCSC <-exactTest(cds,pair=c(group[1],group[length(group)]))
#et.RCSC <-exactTest(cds,pair=c(1,2))
#toptag <-as.matrix(cpmcounts)[(rownames(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1))),]
#write.table(file=paste(name,"_diff.txt",sep=""),cbind(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1),toptag),quote=F,sep="\t")
#return(toptag)
#}
#DE(cpmcounts,allgroup,4,tmp)
#pdf(paste(name,"pdf",sep="."))
#plotMDS(cpmcounts)#, labels=group)  ##### check the replicate distribution
#dev.off
