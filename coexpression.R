#args <- commandArgs(trailingOnly = TRUE)
#########read in the gene names and lncrna names for analysis
#gene.high <- read.delim(file=args[1],header=T,stringsAsFactors = FALSE) 
#genenamepure.high<- as.matrix(rownames(gene.high))
#lnc.high <- read.delim(file=args[2],header=T,stringsAsFactors = FALSE) 
#lncnamepure.high<- as.matrix(rownames(lnc.high))
#lnc.high.del <- read.delim(file=args[3],header=T,stringsAsFactors = FALSE) 
#lncnamepure.high.del <- as.matrix(rownames(lnc.high.del))
#genename <- read.delim(file=args[4],header=F,stringsAsFactors = FALSE) 
#kmean.onlyleaf <- read.delim(file=args[5],header=T,stringsAsFactors = FALSE)
#filename <- args[6]
#number <- args[7]
gene.high <- read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/coexpression/normlnc_expression_gene.high.go.txt",header=T,stringsAsFactors = FALSE) 
genenamepure.high<- as.matrix(rownames(gene.high))
lnc.high <- read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/coexpression/normlnc_expression_lnc.high.txt",header=T,stringsAsFactors = FALSE) 
lncnamepure.high<- as.matrix(rownames(lnc.high))
#lnc.high.del <- read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/coexpression/normlnc_expression_lnc.high.del.txt",header=T,stringsAsFactors = FALSE) 
#lncnamepure.high.del <- as.matrix(rownames(lnc.high.del))
filename="~/sunyd/identify/oryza_rnaseq/SRRpi/coexpression/coexpression/network"
number=6
log2allcount.avs.high <- rbind(gene.high,lnc.high)
legendname <- c("Genes","LncRNAs","Intron","Inter","Exon","Anti","iso","match")
color8 <- c("black","red","blue","green","purple","pink","yellow","grey")
legendname3 <- c("Genes","HC-lncRNAs","Other_lncRNAs")
color3 <- c("blue","red","green")

##################start coexpression analysis############
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(64)
library(flashClust)

tom <- function(datExpr,num){
gene.names=rownames(t(datExpr))
softPower = num;
		#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", power = softPower);
		#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
colnames(TOM) =rownames(TOM) = gene.names
dissTOM=1-TOM
return(dissTOM)
}
##################Coexpression Analysis#################
coexpression.dynamic <- function(datExpr,dissTOM,filename,num){
########## Choosing the soft-thresholding power: analysis of network topology###############
gene.names=rownames(t(datExpr))
softPower = num;
#######Module detection#########
geneTree = flashClust(as.dist(dissTOM),method="average")   #hierarchical clustering of the genes based on the TOM dissimilarity measure
#plot(geneTree, xlab="", sub="",cex=0.3)	#plot the resulting clustering tree (dendrogram)
minModuleSize = 30	#set the minimum module size relatively high
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize) 	# Module identification using dynamic tree cut
print(table(dynamicMods))	#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
print(table(dynamicColors))
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")      ########plot the original color and dynamic tree
module_colors= unique(dynamicColors)
for (color in module_colors){
    module=gene.names[which(dynamicColors==color)]
    write.table(module, paste(filename,color,".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
###########Look at expression patterns of these genes, as they're clustered in dynamicColors########
module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
names(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))       ####print the color order, from the bottom to the top
table(dynamicColors)
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
pdf(paste(filename,"heatmap.pdf",sep=""))
heatmap(t(m),zlim=c(0,1),col=colorRampPalette(c("green","black","red"))(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
dev.off()
return(dynamicColors)

}
coexpression.merge <- function(datExpr,dissTOM,dynamicColors,filename,num){
##########Quantify module similarity by eigengene correlation.merge similar modules#########
gene.names=rownames(t(datExpr))
geneTree = flashClust(as.dist(dissTOM),method="average")   #hierarchical clustering of the genes based on the TOM dissimilarity measure
MEList = moduleEigengenes(datExpr, colors = dynamicColors)       #######calculate eigengenes
MEs = MEList$eigengenes                   
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method = "average")     #######hcluster modules
#plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")    

MEDissThres = 0.25
#abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)    ###merge modules under the threshold line
mergedColors = merge$colors
mergedMEs = merge$newMEs      

#sizeGrWindow(12, 9)
pdf(paste(filename,"cluster.pdf",sep=""))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors           # Rename to moduleColors
colorOrder = c("grey", standardColors(50))              # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs;
#save(MEs, moduleLabels, moduleColors, geneTree, file = "module.color.txt")     # Save module colors and labels for use in subsequent parts

############Extract modules##############
module_colors= unique(mergedColors)
for (color in module_colors){
    module=gene.names[which(mergedColors==color)]
    write.table(module, paste(filename,color, ".merge.txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

###########Look at expression patterns of these genes, as they're clustered in merged colors########
module.order <- unlist(tapply(1:ncol(datExpr),as.factor(mergedColors),I))
names(tapply(1:ncol(datExpr),as.factor(mergedColors),I))       ####print the color order, from the bottom to the top
table(mergedColors)
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
pdf(paste(filename,"heatmap.merge.pdf",sep=""))
heatmap(t(m),zlim=c(0,1),col=colorRampPalette(c("green","black","red"))(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=mergedColors[module.order])
dev.off()
return(mergedColors)
}

#kmean.leaf.log2counts <- rbind(lnc.high[as.matrix(rownames(kmean.onlyleaf)),],gene.high)
kmean.leaf.log2counts <- log2allcount.avs.high
system.time(leaf.tom <- tom(t(kmean.leaf.log2counts),16))
leaf.coexpress.dynamic <- coexpression.dynamic(t(kmean.leaf.log2counts),leaf.tom,filename,16)
leaf.coexpress.merge <- coexpression.merge(t(kmean.leaf.log2counts),leaf.tom,leaf.coexpress.dynamic,filename,16)
####show the corresponding color module for the lncRNA######
colorname <- function(datExpr,dynamicColors,onlyleaf){
	module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))      #######names is color and content is gene num
	extract <- module.order[which(module.order %in% as.matrix(1:nrow(onlyleaf)))]     ####beacause lncRNA is from the front, so only extract the front position and show color
#	print(names(extract))
	gene.names=rownames(t(datExpr))
	lncname <- as.data.frame(names(extract))      ######show the lncRNA color 
	rownames(lncname) <- gene.names[extract]      ######show the color corresponding gene names
	return(lncname)
}

#############show the coresponding color gene names##########
color.location <- function(datExpr,dynamicColors,color,filename){
	color.location <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
	coexpress.color <-  rownames(t(datExpr))[color.location$color]
	write.table(coexpress.color, paste(filename,color,".txt",sep="", sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE))
	} 

kmean.leaf.color <- colorname(t(kmean.leaf.log2counts),leaf.coexpress.dynamic,kmean.onlyleaf)     ####show the corresponding color module for the lncRNA######
write.table(kmean.leaf.color, paste(filename,".coloronly.txt",sep=""), append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(as.matrix(table(leaf.coexpress.dynamic)), paste(filename,".coloronly.all.txt",sep=""), append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
kmean.leaf.color.merge <- colorname(t(kmean.leaf.log2counts),leaf.coexpress.merge,kmean.onlyleaf)     ####show the corresponding color module for the lncRNA######
write.table(kmean.leaf.color.merge, paste(filename,".coloronly.merge.txt",sep=""), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(as.matrix(table(leaf.coexpress.merge)), paste(filename,".coloronly.merge.all.txt",sep=""), append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
