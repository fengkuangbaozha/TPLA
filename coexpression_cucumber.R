allcount <-read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuffnorm75del/isoforms.fpkm_tableorder",row.names=1,header=T)    ###read in the cuffnorm gene count info
combgroup <- read.delim("~/sunyd/identify/cucumber_rnaseq/sample_sheet2com.txt",header=F)     ###combine the same group or tissue together info
name <- combgroup[,3]   ##only extract the group info
genename <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75.combined.gtf-correspond-genename",header=F,stringsAsFactors = FALSE) 
genenamepure <- unique(as.matrix(genename[,2]))
lncname <- read.delim(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.lncnamedel-iuox",header=F)
lncnamedel <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.lncnamedel-cej",header=F,stringsAsFactors = FALSE)
lncnamepure <- unique(as.matrix(lncname[,2]))
lncnamepure.del <- unique(as.matrix(lncnamedel[,2]))

iclass <- unique(as.matrix(lncname[which(lncname[,3]=="i"),2]))
uclass <- unique(as.matrix(lncname[which(lncname[,3]=="u"),2]))
xclass <- unique(as.matrix(lncname[which(lncname[,3]=="x"),2]))
oclass <- unique(as.matrix(lncname[which(lncname[,3]=="o"),2]))
jclass <- unique(as.matrix(lncnamedel[which(lncnamedel[,3]=="j"),2]))
eclass <- unique(as.matrix(lncnamedel[which(lncnamedel[,3]=="="),2]))
cclass <- unique(as.matrix(lncnamedel[which(lncnamedel[,3]=="c"),2]))
legendname3 <- c("Genes","HC-lncRNAs","Other_lncRNAs")
color3 <- c("blue","red","green")

log2countsorigin <- log2(allcount+1)
keep <- apply(allcount,1,function(x){any(x > 1)})    #####only extract FPKM > 1 in at least one sample
cpmcountshigh <- allcount[keep,]
log2countshigh <- log2(cpmcountshigh+1)
genenamepurehigh <- as.matrix(rownames(log2countshigh)[which(rownames(log2countshigh) %in% genenamepure)])
lncnamepurehigh <- as.matrix(rownames(log2countshigh)[which(rownames(log2countshigh) %in% lncnamepure)])
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
log2cpmcounts <- log2(cpmcounts+1)
#######get the lncRNA and gene cpm and log2 counts for heatmap analysis#############
lnccpmcounts <- cpmcounts[which(rownames(cpmcounts) %in% lncnamepure),]                
lnclog2counts <- log2(lnccpmcounts+1)     ###lncRNA log2counts
lnccpmcounts.del <- cpmcounts[which(rownames(cpmcounts) %in% lncnamepure.del),]                
lnclog2counts.del <- log2(lnccpmcounts.del+1)     ###lncRNA log2counts
genecpmcounts <- cpmcounts[which(rownames(cpmcounts) %in% genenamepure),]  
genelog2counts <- log2(genecpmcounts+1)    ####gene log2 counts

############get the correlation number between lncrna and gene###########
lnc_gene <- rbind(lnclog2counts,genelog2counts)      ######lnclog2counts.del
adj.test = adjacency(t(lnc_gene),selectCols = c(1:nrow(lnclog2counts)),type = "signed", power = 4)
adj.test_lnc <- t(adj.test[as.matrix(rownames(genelog2counts)),])     #########only extract lncRNA_gene correlation
cor.genename <- genename[which(as.matrix(genename[,2]) %in% as.matrix(colnames(adj.test_lnc))),c(2,4)]     #########correspond TCONS and Csa
cor.genename <- cor.genename[order(cor.genename[,1]),]               ###########order the name
colnames(adj.test_lnc) <- cor.genename[,2]           ######give the cucumber gff name to the correlation matrix
adj.test_lnc.top <- t(apply(adj.test_lnc,1,function(x){y <- colnames(adj.test_lnc)[order(x,decreasing = TRUE)]; y <- y[1:40]; return(y)}))
write.table(adj.test_lnc.top, "wgcna.correlation.txt", sep="\t", row.names=TRUE, col.names=FALSE,quote=FALSE)

##################start coexpression analysis############
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
library(flashClust)

##################Coexpression Analysis#################
coexpression <- function(datExpr,filename){
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

############Generating adjacency and TOM similarity matrices based on the selected softpower#########
softPower = 6;
		#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", power = softPower);
		#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
colnames(TOM) =rownames(TOM) = gene.names
dissTOM=1-TOM

#######Module detection#########
geneTree = flashClust(as.dist(dissTOM),method="average")   #hierarchical clustering of the genes based on the TOM dissimilarity measure
#plot(geneTree, xlab="", sub="",cex=0.3)	#plot the resulting clustering tree (dendrogram)
minModuleSize = 30	#set the minimum module size relatively high
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize) 	# Module identification using dynamic tree cut
print(table(dynamicMods))	#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
print(table(dynamicColors))
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")      ########plot the original color and dynamic tree

##########Quantify module similarity by eigengene correlation.merge similar modules#########
MEList = moduleEigengenes(datExpr, colors = dynamicColors)       #######calculate eigengenes
MEs = MEList$eigengenes                   
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method = "average")     #######hcluster modules
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")    

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
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
    write.table(module, paste(filename,color, ".txt",sep="."), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

###########Look at expression patterns of these genes, as they're clustered in merged colors########
module.order <- unlist(tapply(1:ncol(datExpr),as.factor(mergedColors),I))
names(tapply(1:ncol(datExpr),as.factor(mergedColors),I))       ####print the color order, from the bottom to the top
table(mergedColors)
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
pdf(paste(filename,"heatmap.merge.pdf",sep=""))
heatmap(t(m),zlim=c(0,1),col=colorRampPalette(c("green","black","red"))(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=mergedColors[module.order])
dev.off()

###########Look at expression patterns of these genes, as they're clustered in dynamicColors########
module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
names(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))       ####print the color order, from the bottom to the top
table(dynamicColors)
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
pdf(paste(filename,"heatmap.pdf",sep=""))
heatmap(t(m),zlim=c(0,1),col=colorRampPalette(c("green","black","red"))(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
return(dynamicColors)
}

kmean.onlygua <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnciuoxkcluster.gua",header=F,stringsAsFactors = FALSE)
kmean.gua.log2counts <- rbind(lnclog2counts[as.matrix(kmean.onlygua[,1]),],genelog2counts)
kmean.gua <- coexpression(t(kmean.gua.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.kmean.gua.")
kmean.onlyphloem <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnciuoxkcluster.phloem",header=F,stringsAsFactors = FALSE)
kmean.phloem.log2counts <- rbind(lnclog2counts[as.matrix(kmean.onlyphloem[,1]),],genelog2counts)
kmean.phloem <- coexpression(t(kmean.phloem.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.kmean.phloem.")
kmean.onlyflower <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnciuoxkcluster.flower",header=F,stringsAsFactors = FALSE)
kmean.flower.log2counts <- rbind(lnclog2counts[as.matrix(kmean.onlyflower[,1]),],genelog2counts)
kmean.flower <- coexpression(t(kmean.flower.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.kmean.flower.")
kmean.onlyfruit <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnciuoxkcluster.fruit",header=F,stringsAsFactors = FALSE)
kmean.fruit.log2counts <- rbind(lnclog2counts[as.matrix(kmean.onlyfruit[,1]),],genelog2counts)
kmean.fruit <- coexpression(t(kmean.fruit.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.kmean.fruit.")
kmean.onlytritrome <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnciuoxkcluster.tritrome",header=F,stringsAsFactors = FALSE)
kmean.tritrome.log2counts <- rbind(lnclog2counts[as.matrix(kmean.onlytritrome[,1]),],genelog2counts)
kmean.tritrome <- coexpression(t(kmean.tritrome.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.kmean.tritrome.")
#hclus.onlygua <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/module_royalblue.sort.txt",header=F,stringsAsFactors = FALSE)
#both.onlygua <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/spec.onlygua.name",header=F,stringsAsFactors = FALSE)
#hclus.gua.log2counts <- rbind(lnclog2counts[as.matrix(hclus.onlygua),],genelog2counts)
#both.gua.log2counts <- rbind(lnclog2counts[as.matrix(both.onlygua),],genelog2counts)
#hclus.gua <- coexpression(t(hclus.gua.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.hclus.gua.")
#both.gua <- coexpression(t(both.gua.log2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.both.gua.")

iuox <- coexpression(t(lnclog2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.lnciuox.")
merge = mergeCloseModules(t(lnclog2counts), iuox, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors
write.table(gene.names[which(mergedColors=="black")], "heatmap/hclus.vascular.cluster.txt", sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)  ##write tissue specific color
write.table(gene.names[which(mergedColors=="brown")], "heatmap/hclus.phloem.cluster.txt", sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
write.table(gene.names[which(mergedColors=="darkorange")], "heatmap/hclus.flower.cluster.txt", sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
write.table(gene.names[which(mergedColors=="lightgreen")], "heatmap/hclus.gua.cluster.txt", sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

gene <- coexpression(t(genelog2counts),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.gene.")
lnc_gene <- coexpression(t(lnc_gene),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.lnc_gene.")
lnc <- coexpression(t(rbind(lnclog2counts,lnclog2counts.del)),"/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.lnc.")

####show the corresponding color module for the lncRNA######
colorname <- function(datExpr,dynamicColors,onlygua){
	module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))      #######names is color and content is gene num
	extract <- module.order[which(module.order %in% as.matrix(1:nrow(onlygua)))]     ####beacause lncRNA is from the front, so only extract the front position and show color
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

kmean.gua.color <- colorname(t(kmean.gua.log2counts),kmean.gua,kmean.onlygua)     ####show the corresponding color module for the lncRNA######
write.table(kmean.gua.color, "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyguakcluster.coloronly.txt", append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(as.matrix(table(kmean.gua)), "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyguakcluster.coloronly.all.txt", append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
kmean.onlygua.divide <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyguakcluster.txt",header=T,stringsAsFactors = FALSE,row.names=1)
kmean.gua.merge <- merge(kmean.onlygua.divide,kmean.gua.color,by="row.names",all.x=TRUE,sort=FALSE)       ###merge the color file and the cluster file 
write.table(kmean.gua.merge, "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyguakcluster.color.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
gua.color.location <- tapply(1:ncol(t(kmean.gua.log2counts)),as.factor(kmean.gua),I)
gua.coexpress.turquoise <- rownames(kmean.gua.log2counts)[gua.color.location$turquoise]        ############this is for kmeans heatmap specific genes
write.table(gua.coexpress.turquoise, paste("wgcna.kmean.gua.spec.module_","turquoise", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
gua.coexpress.greenyellow <- rownames(kmean.gua.log2counts)[gua.color.location$greenyellow]        ############this is for kmeans heatmap specific genes
write.table(gua.coexpress.greenyellow, paste("wgcna.kmean.gua.spec.module_","greenyellow.pedicle", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
gua.coexpress.royalblue <- rownames(kmean.gua.log2counts)[gua.color.location$royalblue]        ############this is for kmeans heatmap specific genes
write.table(gua.coexpress.royalblue, paste("wgcna.kmean.gua.spec.module_","royalblue.fruit", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

kmean.flower.color <- colorname(t(kmean.flower.log2counts),kmean.flower,kmean.onlyflower)     ####show the corresponding color module for the lncRNA######
write.table(kmean.flower.color, "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyflowerkcluster.coloronly.txt", append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(as.matrix(table(kmean.flower)), "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyflowerkcluster.coloronly.all.txt", append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
kmean.onlyflower.divide <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyflowerkcluster.txt",header=T,stringsAsFactors = FALSE,row.names=1)
kmean.flower.merge <- merge(kmean.onlyflower.divide,kmean.flower.color,by="row.names",all.x=TRUE,sort=FALSE)       ###merge the color file and the cluster file 
write.table(kmean.flower.merge, "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyflowerkcluster.color.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
flower.color.location <- tapply(1:ncol(t(kmean.flower.log2counts)),as.factor(kmean.flower),I)
flower.coexpress.blue <- rownames(kmean.flower.log2counts)[flower.color.location$blue]        ############this is for kmeans heatmap specific genes
write.table(flower.coexpress.blue, paste("wgcna.kmean.flower.spec.module_","blue", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
flower.coexpress.green <- rownames(kmean.flower.log2counts)[flower.color.location$green]        ############this is for kmeans heatmap specific genes
write.table(flower.coexpress.green, paste("wgcna.kmean.flower.spec.module_","green", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
flower.coexpress.greenyellow <- rownames(kmean.flower.log2counts)[flower.color.location$greenyellow]        ############this is for kmeans heatmap specific genes
write.table(flower.coexpress.greenyellow, paste("wgcna.kmean.flower.spec.module_","greenyellow", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

kmean.phloem.color <- colorname(t(kmean.phloem.log2counts),kmean.phloem,kmean.onlyphloem)     ####show the corresponding color module for the lncRNA######
write.table(kmean.phloem.color, "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyphloemkcluster.coloronly.txt", append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(as.matrix(table(kmean.phloem)), "~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/normlnconlyphloemkcluster.coloronly.all.txt", append=T, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
phloem.color.location <- tapply(1:ncol(t(kmean.phloem.log2counts)),as.factor(kmean.phloem),I)
phloem.coexpress.brown <- rownames(kmean.phloem.log2counts)[phloem.color.location$brown]        ############this is for kmeans heatmap specific genes
write.table(phloem.coexpress.brown, paste("wgcna.kmean.phloem.spec.module_","brown", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
phloem.coexpress.yellow <- rownames(kmean.phloem.log2counts)[phloem.color.location$yellow]        ############this is for kmeans heatmap specific genes
write.table(phloem.coexpress.yellow, paste("wgcna.kmean.phloem.spec.module_","yellow", ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)

#discard the unassigned genes, and focus on the rest to build module########
#restGenes= (dynamicColors != "grey")
#diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)
#colnames(diss1) =rownames(diss1) = gene.names[restGenes]
#hier1=flashClust(as.dist(diss1), method="average" )
#plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
	#set the diagonal of the dissimilarity to NA 
#diag(dissTOM) = NA;
	#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
#pdf("~/sunyd/identify/cucumber_rnaseq/SRRadd/tomplot.pdf")
#sizeGrWindow(7,7)
#TOMplot(dissTOM, geneTree, as.character(dynamicColors))
#dev.off()
