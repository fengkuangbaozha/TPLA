#args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
#allcount.avs.high <-read.delim(file=args[1],row.names=1,header=F)
#allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
#allgroup <- allgroup[,2]
#filename <- args[3]

allcount <-read.delim("~/sunyd/identify/tomato_rnaseq/SRRm82/cuffnorm.light/isoforms.fpkm_table",row.names=1,header=T)    ###read in the cuffnorm gene count info
combgroup <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRm82/cuffnorm.light/samples.table.group.com",header=F)     ###combine the same group or tissue together info
allcount1 = allcount
colnames(allcount1) <- combgroup[,1]   ##put the group name for the colnames of the data

############combine the replicate counts together and use mean value##################
allcount.avs = data.frame(t(apply(allcount1,1,tapply,names(allcount1),mean)))    ########calculate the group mean fpkm value
write.table(allcount.avs,"normlnc_expression_mean.txt",quote=F,sep="\t",col.names=TRUE,row.names=TRUE)
explnccpm <- t(apply(allcount.avs[which(rownames(allcount.avs) %in% lncnamepure),],1,summary))   ######find fpkm max, min, medium, mean
explnccpm <- explnccpm[,c(1,3,4,6)]
write.table(explnccpm,"normlnc_expression_summary.txt",quote=F,sep="\t",col.names=TRUE,row.names=TRUE)

#########read in the gene names and lncrna names for analysis
genename <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRm82/cuff85.combined.gtf-correspond-mrnaname",header=F,stringsAsFactors = FALSE) 
genenamepure <- unique(as.matrix(genename[,2]))
lncname <- read.delim(file = "~/sunyd/identify/tomato_rnaseq/SRRm82/cuff85-CPC_left.lncname.iuox.del",header=F)
lncnamedel <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRm82/cuff85-CPC_left.lncname.ej.del",header=F,stringsAsFactors = FALSE)
lncnamepure <- unique(as.matrix(lncname[,2]))
lncnamepure.del <- unique(as.matrix(lncnamedel[,2]))
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

##########extract the fpkm value high gene name##############
log2countsorigin <- log2(allcount.avs+1)
keep <- apply(allcount.avs,1,function(x){any(x > 1)})    #####only extract FPKM > 1 in at least one sample
allcount.avs.high <- allcount.avs[keep,]
log2allcount.avs.high <- log2(allcount.avs.high+1)
genenamepure.high <- as.matrix(rownames(log2allcount.avs.high)[which(rownames(log2allcount.avs.high) %in% genenamepure)])   #####extract the high expressed gene name
lncnamepure.high <- as.matrix(rownames(log2allcount.avs.high)[which(rownames(log2allcount.avs.high) %in% lncnamepure)])
lncnamepure.high.del <- as.matrix(rownames(log2allcount.avs.high)[which(rownames(log2allcount.avs.high) %in% lncnamepure.del)])
gene.high <- log2allcount.avs.high[which(rownames(log2allcount.avs.high) %in% genenamepure),]  ######extract the high expressed genename expression, log 2 value
lnc.high <- log2allcount.avs.high[which(rownames(log2allcount.avs.high) %in% lncnamepure),]
lnc.high.del <- log2allcount.avs.high[which(rownames(log2allcount.avs.high) %in% lncnamepure.del),]
############calculate maximum cpm value for gene expression plot analysis###########
expgenehigh <- apply(gene.high,1,max)
explnchigh <- apply(lnc.high,1,max)
explnchigh.del <- apply(lnc.high.del,1,max)
valu <- c(data.matrix(expgenehigh),data.matrix(explnchigh),data.matrix(explnchigh.del))
tim <- c(length(expgenehigh),length(explnchigh),length(explnchigh.del))   
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
pdf("normlnc_expression.pdf")
boxplot(values ~ vars, data = df, col = color3,main="Transcript expression distribution (Maximum)")
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
keep.all <- apply(allcount,1,function(x){any(x > 1)})
allcounts.high <- allcount[keep.all,]
tcounts <- t(allcounts.high)
pca.total <- prcomp(log2(tcounts+1), retx=TRUE)
pdf(paste("normhtseq","pca",".pdf",sep=""))
	c <- round(100*summary(pca.total)$importance[2,1],digits=2)
	d <- round(100*summary(pca.total)$importance[2,2],digits=2)
	plot(pca.total$x[,1:2], pch=c(1:25), col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),main="PCA Plot of Samples")
	legend("topleft",cex=0.5,border=F, c(row.names(pca.total$x)),pch=c(1:25),col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)),bty="n")
dev.off()

#######get the lncRNA and gene cpm and log2 counts for heatmap analysis#############
#genelog2counts.order <- genelog2counts[,c(2,3,4,5,1,10,11,6,7,8,9,12:19)]
#lnclog2counts.order <- lnclog2counts[,c(2,3,4,5,1,10,11,6,7,8,9,12:19)]
#lnclog2counts.del.order <- lnclog2counts.del[,c(2,3,4,5,1,10,11,8,9,6,7,12:19)]
#colnames(genelog2counts.order)[16] <- "leaf"
#colnames(lnclog2counts.order)[16] <- "leaf"
#colnames(lnclog2counts.del.order)[16] <- "leaf"

############calculate tissue specific score###########3
ts <- ROKU(log2allcount.avs.high, upper.limit = 0.25, sort = FALSE)
tsorigin <- as.matrix(ts$H)
write.table(tsorigin,"normlnc_tissuespecific.txt",quote=F,sep="\t",col.names=FALSE,row.names=TRUE)
tsgenehigh <- as.matrix(tsorigin[which(rownames(tsorigin) %in% genenamepure.high),])
tslnchigh <- as.matrix(tsorigin[which(rownames(tsorigin) %in% lncnamepure.high),])       
tslnchigh.del <- as.matrix(tsorigin[which(rownames(tsorigin) %in% lncnamepure.high.del),])       

valu <- c(tsgenehigh,tslnchigh,tslnchigh.del)
tim <- c(nrow(tsgenehigh),nrow(tslnchigh),nrow(tslnchigh.del))
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
pdf("normlnc_tissuespecific.pdf")
boxplot(values ~ vars, data = df, col = color3,main="Transcript tissue specific score distribution")
dev.off()
#itslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% iclass),]) 
#utslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% uclass),]) 
#otslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% oclass),]) 
#xtslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% xclass),]) 
#jtslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% jclass),]) 
#etslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% eclass),]) 
#valu <- c(tsgenehigh,tslnchigh,itslnchigh,utslnchigh,otslnchigh,xtslnchigh,jtslnchigh,etslnchigh)
#tim <- c(nrow(tsgenehigh),nrow(tslnchigh),nrow(itslnchigh),nrow(utslnchigh),nrow(otslnchigh),nrow(xtslnchigh),nrow(jtslnchigh),nrow(etslnchigh))

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
	pdf(paste(filename,"kheatmap.pdf",sep=""))
	col1 <- colorRampPalette(c("green","black","red"))(100)  ##set the color for heatmap
	pheatmap(t(heat1),scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 8)
	dev.off()
	####Hierachical clustering#############
#	hclus <- agnes(t(counts),method = "gaverage")
#	pdf(paste(filename,"hheatmap.pdf",sep=""))
#	pheatmap(hclus$height,scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 8)
#	dev.off()
}
#heatmapfunc(log2allcount.avs.high,"lnc_gene")
heatmapfunc(gene.high,6,"normgene")
heatmapfunc(lnc.high,6,"normlnc")
#heatmapfunc(lr,6,"normlnclow")
#heatmapfunc(hr,6,"normlnchigh")
#heatmapfunc(lnclog2counts[which(rownames(lnclog2counts) %in% iclass),],6,"normlnci")
#heatmapfunc(lnclog2counts[which(rownames(lnclog2counts) %in% uclass),],6,"normlncu")
#heatmapfunc(lnclog2counts[which(rownames(lnclog2counts) %in% oclass),],6,"normlnco")
#heatmapfunc(lnclog2counts[which(rownames(lnclog2counts) %in% xclass),],6,"normlncx")
#heatmapfunc(lnclog2counts[which(rownames(lnclog2counts) %in% iclass | rownames(lnclog2counts) %in% uclass | rownames(lnclog2counts) %in% xclass),],6,"normlnciux")   #####keep only i,u,o,x classcode  (good)
#heatmapfunc(lnclog2counts[which(rownames(lnclog2counts) %in% iclass | rownames(lnclog2counts) %in% uclass | rownames(lnclog2counts) %in% xclass | rownames(lnclog2counts) %in% oclass),],6,"normlnciuox")   #####keep only i,u,o,x classcode  (good)
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% jclass),],6,"normlncdelj")
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% jclass | rownames(lnclog2counts) %in% hkrna | rownames(lnclog2counts) %in% te),],6,"normlncdeljhkte")
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% jclass | rownames(lnclog2counts) %in% eclass | rownames(lnclog2counts) %in% cclass | rownames(lnclog2counts) %in% hkrna | rownames(lnclog2counts) %in% te),],6,"normlnciuxdelhkte")

##################divide different tissues and compare their within differences###########
partheat <- function(x,col,filename,num){
	x <- x[,col]
	keep.x <- apply(x,1,function(x){any(x > 1)})
	x <- x[keep.x,]
	print(nrow(x))
	heatmapfunc(x,num,filename)
	}
onlygua <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/normlnciuoxkcluster.gua",header=F,stringsAsFactors = FALSE)
partheat(lnclog2counts[as.matrix(onlygua[,1]),],c(3,4,5),"normlnconlygua",3)
onlyfruit <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/normlnciuoxkcluster.fruit",header=F,stringsAsFactors = FALSE)
partheat(lnclog2counts[as.matrix(onlyfruit[,1]),],c(6,7),"normlnconlyfruit",2)
onlytritrome <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/normlnciuoxkcluster.tritrome",header=F,stringsAsFactors = FALSE)
partheat(lnclog2counts[as.matrix(onlytritrome[,1]),],c(8,9),"normlnconlytritrome",2)

            ##########flower xiaolanzhang analysis################
flowercount <-read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuffnorm75/isoforms.fpkm_table-flower",row.names=1,header=T)    ###read in the cuffnorm gene count info
keep.flower <- apply(flowercount,1,function(x){any(x > 1)})
flowercounthigh <- flowercount[keep.flower,]
fallcount.avs.high <- matrix(NA, nrow=nrow(flowercounthigh),ncol=11)   ###combine the replicate together
rownames(fallcount.avs.high) <- rownames(flowercounthigh)
colnames(fallcount.avs.high) <- c("Normal_0_DAP","Super_0_DAP","Super_1_DAP","Normal_1_DAP","Super_2_DAP","Normal_2_DAP","Super_4_DAP","Normal_4_DAP","Super_5_DAP","Super_6_DAP","Super_8_DAP")
fallcount.avs.high[,1] <- apply(flowercounthigh[,1:3],1,mean)
fallcount.avs.high[,2] <- apply(flowercounthigh[,4:6],1,mean)
fallcount.avs.high[,3] <- apply(flowercounthigh[,7:9],1,mean)
fallcount.avs.high[,4] <- apply(flowercounthigh[,10:12],1,mean)
fallcount.avs.high[,5] <- apply(flowercounthigh[,13:15],1,mean)
fallcount.avs.high[,6] <- apply(flowercounthigh[,16:18],1,mean)
fallcount.avs.high[,7] <- apply(flowercounthigh[,19:21],1,mean)
fallcount.avs.high[,8] <- apply(flowercounthigh[,22:24],1,mean)
fallcount.avs.high[,9] <- apply(flowercounthigh[,25:27],1,mean)
fallcount.avs.high[,10] <- apply(flowercounthigh[,28:30],1,mean)
fallcount.avs.high[,11] <- apply(flowercounthigh[,31:33],1,mean)
flncallcount.avs.high <- fallcount.avs.high[which(rownames(fallcount.avs.high) %in% lncnamepure),]
flnclog2counts <- log2(flncallcount.avs.high+1)     ###lncRNA log2counts

onlyflower <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/normlnciuoxkcluster.flower",header=F,stringsAsFactors = FALSE)
heatmapfunc(flnclog2counts[which(rownames(flnclog2counts) %in% onlyflower[,1]),],6,"normlnconlyflower")                                                                 
#partheat(flnclog2counts[which(rownames(flnclog2counts) %in% onlyflower[,1]),],c(1:11),"normlnconlyflower",6)
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% eclass | rownames(lnclog2counts) %in% jclass | rownames(lnclog2counts) %in% cclass | rownames(lnclog2counts) %in% as.matrix(hkrna) | rownames(lnclog2counts) %in% as.matrix(te)),],6,"normlncdelejc") ###keep only i,o,u,x and delete hkrna
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% as.matrix(delrna[,2])),],6,"normlncdel")
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% as.matrix(hkrna) | rownames(lnclog2counts) %in% as.matrix(mirna)),],"normlncdelhk")
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% as.matrix(hkrna)),],6,"normlncdelhk")
#heatmapfunc(lnclog2counts[-which(rownames(lnclog2counts) %in% as.matrix(TE)),],6,"normlncdelTE")

############estimate the common distribution of the data###########
#cds <- estimateCommonDisp(cds,verbose=TRUE)
	#names(cds)
	#sqrt(cds$common.dispersion)   ##the common distribution value
#cds <- estimateTagwiseDisp(cds)
	#plotBCV(cds)
	#abline(h=sqrt(cds$common.dispersion), col="firebrick", lwd=3)

#et.RCSC <-exactTest(cds,pair=c(group[1],group[length(group)]))
#et.RCSC <-exactTest(cds,pair=c(1,2))
#toptag <-as.matrix(allcount.avs.high)[(rownames(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1))),]
#write.table(file=paste(name,"_diff.txt",sep=""),cbind(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1),toptag),quote=F,sep="\t")
#return(toptag)
#}
#DE(allcount.avs.high,allgroup,4,tmp)
#pdf(paste(name,"pdf",sep="."))
#plotMDS(allcount.avs.high)#, labels=group)  ##### check the replicate distribution
#dev.off
########delete geneid and transcriptid, only keep pure name############
pureext <- function(x){
        class1 <- strsplit(as.character(x[,1])," ")        
        class2 <- strsplit(as.character(x[,2])," ")
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
#genenamepure <- unique(pureext(genename))  ###only keep pure gene or transcript name
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
 colnames(z) <- c("gene_id","transcript_id","class_code","oId","explnchighgth","start","end","strand")
 return(z)
}
#lncgenename <- nameext(lncgtf)  ######calculate the lncRNA info, extract name info
#lncnamepure <- unique(pureext(lncgenename))  ###only keep pure gene or transcript name

#############classify gene class depending on genomic location#############
codeext <- function(x,y){
    class <- x[which(x[,3]== y),]
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
