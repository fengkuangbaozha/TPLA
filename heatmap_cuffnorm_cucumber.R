#args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
#cpmcounts <-read.delim(file=args[1],row.names=1,header=F)
#allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
#allgroup <- allgroup[,2]
#filename <- args[3]

allcount <-read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuffnorm75del/isoforms.fpkm_tableorder",row.names=1,header=T)    ###read in the cuffnorm gene count info
combgroup <- read.delim("~/sunyd/identify/cucumber_rnaseq/sample_sheet2com.txt",header=F)     ###combine the same group or tissue together info
name <- combgroup[,3]   ##only extract the group info
genename <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75.combined.gtf-correspond-mrnaname-ano",header=F,stringsAsFactors = FALSE) 
genenamepure <- unique(as.matrix(genename[,2]))
lncname <- read.delim(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.lncnamedel-iuox",header=F)
lncnamedel <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.lncnamedel-cej",header=F,stringsAsFactors = FALSE)
lncnamepure <- unique(as.matrix(lncname[,2]))
lncnamepure.del <- unique(as.matrix(lncnamedel[,2]))
#delrna <- read.delim(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75_hkte.name",header=F)
#hkrna <- as.matrix(delrna[which(delrna[,4]=="HousekeepRNA"),2])
#mirna <- as.matrix(delrna[which(delrna[,4]=="SmallRNA"),2])
#te <- read.delim(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75_TE.name",header=F)
#te <- as.matrix(te[,2])
#allcount <-read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuffnorm75/genes.count_tablegroup",row.names=1,header=T)    ###read in the cuffnorm gene count info
#allgroup <- read.delim("~/sunyd/identify/cucumber_rnaseq/sample_sheet2.txt",header=F)     ####the group info, second col is group 
#allgroup <- allgroup[,2]        ##only extract the group info
#name <- c(allgroup[1],paste(allgroup[2:4],1:3,sep=""),paste(allgroup[5:7],1:3,sep=""),paste(allgroup[8:20],1:13,sep=""),paste(allgroup[21:22],1:2,sep=""),paste(allgroup[23:24],1:2,sep=""),paste(allgroup[25:34],1:10,sep=""),paste(allgroup[35],1,sep=""),paste(allgroup[36:37],1:2,sep=""),paste(allgroup[38],1,sep=""),paste(allgroup[39],1,sep=""),paste(allgroup[40:42],1:3,sep=""),paste(allgroup[43:75],1:33,sep=""))
#name <- c(allgroup[1],paste(allgroup[2:4],1:3,sep=""),paste(allgroup[5:7],1:3,sep=""),paste(allgroup[8:20],1:13,sep=""),paste(allgroup[21:28],1:8,sep=""),paste(allgroup[29],1,sep=""),paste(allgroup[30:31],1:2,sep=""),paste(allgroup[32],1,sep=""),paste(allgroup[33],1,sep=""),paste(allgroup[34:36],1:3,sep=""),paste(allgroup[37:69],1:33,sep=""))
#attr <-read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuffnorm75/genes.attr_table",header=T)
#genename <- attr[grep("Csa",attr$gene_short_name),c(1,5)]
#lncgtf <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.gtf",header=F)   ##lncRNA name and classification
iclass <- unique(as.matrix(lncname[which(lncname[,3]=="i"),2]))
uclass <- unique(as.matrix(lncname[which(lncname[,3]=="u"),2]))
xclass <- unique(as.matrix(lncname[which(lncname[,3]=="x"),2]))
oclass <- unique(as.matrix(lncname[which(lncname[,3]=="o"),2]))
jclass <- unique(as.matrix(lncnamedel[which(lncnamedel[,3]=="j"),2]))
eclass <- unique(as.matrix(lncnamedel[which(lncnamedel[,3]=="="),2]))
cclass <- unique(as.matrix(lncnamedel[which(lncnamedel[,3]=="c"),2]))
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

log2countsorigin <- log2(allcount+1)
keep <- apply(allcount,1,function(x){any(x > 1)})    #####only extract FPKM > 1 in at least one sample
cpmcountshigh <- allcount[keep,]
log2countshigh <- log2(cpmcountshigh+1)
genenamepure.high <- as.matrix(rownames(log2countshigh)[which(rownames(log2countshigh) %in% genenamepure)])
lncnamepure.high <- as.matrix(rownames(log2countshigh)[which(rownames(log2countshigh) %in% lncnamepure)])
lncnamepure.high.del <- as.matrix(rownames(log2countshigh)[which(rownames(log2countshigh) %in% lncnamepure.del)])

#cds <- DGEList(allcount,group = allgroup)
#cds <- calcnormFactors(cds)
##############get the counts and log2counts for all lncRNA and genes ###########
#cpmorigin <- cpm(cds$counts)
#log2countsorigin <- log2(cpmorigin+1)
##############get the counts and log2counts for high expressed lncRNA and genes ########### 
#keep <- rowSums(cpm(cds)>1)>=1   ####value > 1 and at least $num numbers, high expression value ones
#cds <-cds[keep,]      
#print(dim(cds))
#cpmcountshigh <- cpm(cds$counts) ###get the count number
#log2countshigh <- log2(cpmcountshigh+1)
#which(rownames(log2countshigh) %in% unique(lncnamepure[,1])) <- rownames(cpmcountshigh)[grep("XLOC",rownames(cpmcountshigh))]   #####only extract lncRNA name rows and delete gene name rows
#which(rownames(log2countshigh) %in% unique(genenamepure[,1])) <- rownames(cpmcountshigh)[grep("Csa",rownames(cpmcountshigh))]   #####only extract lncRNA name rows and delete gene name rows
############calculate maximum cpm value for gene expression plot analysis###########
expgeneorigin <- apply(log2countsorigin[genenamepure,],1,max)
explncorigin <- apply(log2countsorigin[lncnamepure,],1,max)
expgenehigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% genenamepure),],1,max)
explnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% lncnamepure),],1,max)
explnchigh.del <- apply(log2countshigh[which(rownames(log2countshigh) %in% lncnamepure.del),],1,max)
valu <- c(data.matrix(expgenehigh),data.matrix(explnchigh),data.matrix(explnchigh.del))
tim <- c(length(expgenehigh),length(explnchigh),length(explnchigh.del))   
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
pdf("normlnc_expression.pdf")
boxplot(values ~ vars, data = df, col = color3,main="Transcript expression distribution (Maximum)")
dev.off()
#iexplnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% iclass),],1,max)
#uexplnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% uclass),],1,max)
#oexplnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% oclass),],1,max)
#xexplnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% xclass),],1,max)
#jexplnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% jclass),],1,max)
#eexplnchigh <- apply(log2countshigh[which(rownames(log2countshigh) %in% eclass),],1,max)

plot(density(explnchigh),type = "l", xlab = "log2FPKM",ylab = "Density",col="red",main="Expression level distribution (Maximum)")
lines(density(expgenehigh),col = "black")
lines(density(iexplnchigh),col = "blue")
lines(density(uexplnchigh),col = "green")
lines(density(oexplnchigh),col = "purple")
lines(density(xexplnchigh),col = "pink")
lines(density(jexplnchigh),col = "yellow")
lines(density(eexplnchigh),col = "grey")
legend("topright",legend=legendname,col=color8,bg="white",lwd=2)

#itslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% iclass),]) 
#utslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% uclass),]) 
#otslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% oclass),]) 
#xtslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% xclass),]) 
#jtslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% jclass),]) 
#etslnchigh <- as.matrix(tslnchigh[which(rownames(tslnchigh) %in% eclass),]) 
#valu <- c(tsgenehigh,tslnchigh,itslnchigh,utslnchigh,otslnchigh,xtslnchigh,jtslnchigh,etslnchigh)
#tim <- c(nrow(tsgenehigh),nrow(tslnchigh),nrow(itslnchigh),nrow(utslnchigh),nrow(otslnchigh),nrow(xtslnchigh),nrow(jtslnchigh),nrow(etslnchigh))


####################the PCA analysis#################
tcounts <- t(cpmcountshigh)
pca.total <- prcomp(log2(tcounts+1), retx=TRUE)
pdf(paste("normhtseq","pca",".pdf",sep=""))
	c <- round(100*summary(pca.total)$importance[2,1],digits=2)
	d <- round(100*summary(pca.total)$importance[2,2],digits=2)
	plot(pca.total$x[,1:2], pch=c(1:25), col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)), xlab=paste("PC1(",c,"% Proportion of Variance)"),ylab=paste("PC2(",d,"%) Proportion of Variance"),main="PCA Plot of Samples")
	legend("topleft",cex=0.5,border=F, c(row.names(pca.total$x)),pch=c(1:25),col=c(rep("black",9),rep("red",13),rep("blue",8),rep("green",7),rep("purple",33),rep("pink",3),rep("darkgrey",2)),bty="n")
dev.off()

############combine the replicate counts together and use mean value##################
cpmcounts <- matrix(NA, nrow=nrow(cpmcountshigh),ncol=19)   ###combine the replicate together
rownames(cpmcounts) <- rownames(cpmcountshigh)
cpmcounts[,1] <- apply(cpmcountshigh[,1:3],1,mean)
cpmcounts[,2] <- apply(cpmcountshigh[,4:6],1,mean)
cpmcounts[,3] <- apply(cpmcountshigh[,7:8],1,mean)
cpmcounts[,4] <- apply(cpmcountshigh[,9:10],1,mean)
cpmcounts[,5] <- apply(cpmcountshigh[,11:12],1,mean)
cpmcounts[,6] <- apply(cpmcountshigh[,13:14],1,mean) 
cpmcounts[,7] <- apply(cpmcountshigh[,15:16],1,mean) 
cpmcounts[,8] <- apply(cpmcountshigh[,17:18],1,mean) 
cpmcounts[,9] <- apply(cpmcountshigh[,19:20],1,mean)
cpmcounts[,10] <- apply(cpmcountshigh[,21:26],1,mean)
cpmcounts[,11] <- apply(cpmcountshigh[,27:32],1,mean)
cpmcounts[,12] <- cpmcountshigh[,33]
cpmcounts[,13] <- cpmcountshigh[,34]
cpmcounts[,14] <- apply(cpmcountshigh[,35:37],1,mean)
cpmcounts[,15] <- apply(cpmcountshigh[,38:39],1,mean)
cpmcounts[,16] <- cpmcountshigh[,40]
cpmcounts[,17] <- cpmcountshigh[,41]
cpmcounts[,18] <- cpmcountshigh[,42]
cpmcounts[,19] <- apply(cpmcountshigh[,43:44],1,mean)
#cpmcounts <- as.data.frame(cpmcounts[,-1:-3]) #####delete gua,ba,bing
cpmcounts <- as.data.frame(cpmcounts)    #######combined cpmcounts for all tissue samples
colnames(cpmcounts) <- name
write.table(cpmcounts,"normlnc_expression_mean.txt",quote=F,sep="\t",col.names=TRUE,row.names=TRUE)
explnccpm <- t(apply(cpmcounts[which(rownames(cpmcounts) %in% lncnamepure),],1,summary))   ######find fpkm max, min, medium, mean
explnccpm <- explnccpm[,c(1,3,4,6)]
write.table(explnccpm,"normlnc_expression_summary.txt",quote=F,sep="\t",col.names=TRUE,row.names=TRUE)
log2cpmcounts <- log2(cpmcounts+1)
#######get the lncRNA and gene cpm and log2 counts for heatmap analysis#############
lnccpmcounts <- cpmcounts[which(rownames(cpmcounts) %in% lncnamepure),]                
lnclog2counts <- log2(lnccpmcounts+1)     ###lncRNA log2counts
lnccpmcounts.del <- cpmcounts[which(rownames(cpmcounts) %in% lncnamepure.del),]                
lnclog2counts.del <- log2(lnccpmcounts.del+1)     ###lncRNA log2counts
genecpmcounts <- cpmcounts[which(rownames(cpmcounts) %in% genenamepure),]  
genelog2counts <- log2(genecpmcounts+1)    ####gene log2 counts

genelog2counts.order <- genelog2counts[,c(2,3,4,5,1,10,11,6,7,8,9,12:19)]
lnclog2counts.order <- lnclog2counts[,c(2,3,4,5,1,10,11,6,7,8,9,12:19)]
lnclog2counts.del.order <- lnclog2counts.del[,c(2,3,4,5,1,10,11,8,9,6,7,12:19)]
colnames(genelog2counts.order)[16] <- "leaf"
colnames(lnclog2counts.order)[16] <- "leaf"
colnames(lnclog2counts.del.order)[16] <- "leaf"

############calculate tissue specific score###########3
ts <- ROKU(log2cpmcounts, upper.limit = 0.25, sort = FALSE)
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


###identigy kmeans value#######
kmidentify <- function(x){
wss <- (nrow(x)-1)*sum(apply(x,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(x,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
	}
kmidentify(lnccpmcounts)
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
#heatmapfunc(log2cpmcounts,"lnc_gene")
heatmapfunc(genelog2counts.order,6,"normgene")
heatmapfunc(rbind(lnclog2counts.order,lnclog2counts.del.order),6,"normlnc")
heatmapfunc(lnclog2counts.order,6,"normlnciuox")
heatmapfunc(lnclog2counts.del.order,6,"normlncdel")
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
fcpmcounts <- matrix(NA, nrow=nrow(flowercounthigh),ncol=11)   ###combine the replicate together
rownames(fcpmcounts) <- rownames(flowercounthigh)
colnames(fcpmcounts) <- c("Normal_0_DAP","Super_0_DAP","Super_1_DAP","Normal_1_DAP","Super_2_DAP","Normal_2_DAP","Super_4_DAP","Normal_4_DAP","Super_5_DAP","Super_6_DAP","Super_8_DAP")
fcpmcounts[,1] <- apply(flowercounthigh[,1:3],1,mean)
fcpmcounts[,2] <- apply(flowercounthigh[,4:6],1,mean)
fcpmcounts[,3] <- apply(flowercounthigh[,7:9],1,mean)
fcpmcounts[,4] <- apply(flowercounthigh[,10:12],1,mean)
fcpmcounts[,5] <- apply(flowercounthigh[,13:15],1,mean)
fcpmcounts[,6] <- apply(flowercounthigh[,16:18],1,mean)
fcpmcounts[,7] <- apply(flowercounthigh[,19:21],1,mean)
fcpmcounts[,8] <- apply(flowercounthigh[,22:24],1,mean)
fcpmcounts[,9] <- apply(flowercounthigh[,25:27],1,mean)
fcpmcounts[,10] <- apply(flowercounthigh[,28:30],1,mean)
fcpmcounts[,11] <- apply(flowercounthigh[,31:33],1,mean)
flnccpmcounts <- fcpmcounts[which(rownames(fcpmcounts) %in% lncnamepure),]
flnclog2counts <- log2(flnccpmcounts+1)     ###lncRNA log2counts

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
#toptag <-as.matrix(cpmcounts)[(rownames(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1))),]
#write.table(file=paste(name,"_diff.txt",sep=""),cbind(subset(topTags(et.RCSC,n=10000)$table,FDR<0.05&abs(logFC) >=  1),toptag),quote=F,sep="\t")
#return(toptag)
#}
#DE(cpmcounts,allgroup,4,tmp)
#pdf(paste(name,"pdf",sep="."))
#plotMDS(cpmcounts)#, labels=group)  ##### check the replicate distribution
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
