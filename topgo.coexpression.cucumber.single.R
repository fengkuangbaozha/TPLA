                ###########make GO annotation for lncrna coexpressed gene###########
library(topGO)
library(ALL)
library(parallel)
anoGO.new <- readMappings("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuc.new.map")
anoGO.old <- readMappings("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuc.map")

lncrnaInterestingGenes <- as.matrix(read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/wgcna.correlation.spec.txt",header=F,row.names=1,stringsAsFactors=FALSE))[,1:20]
lncrnaInterestingGenes.list <- split(lncrnaInterestingGenes,rownames(lncrnaInterestingGenes))  ######convert to a list to perform lapply
#selected <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/heatmap/explnc.all",header=F,row.names=1,stringsAsFactors=FALSE)
#s <- c(1:19,335:352,582:600,757:769,943:968,2868:2891,3732:3767,4296:4316,4896:4925,5281:5301)
#selected2 <- as.matrix(rownames(selected[s,]))
#selected.InterestingGenes <- lncrnaInterestingGenes[selected2,]
#selected.InterestingGenes.list <- split(selected.InterestingGenes,rownames(selected.InterestingGenes))  ######convert to a list to perform lapply
ano <- function(myInterestingGenes,anoGO,filename){
mat <- data.frame(BP_term=NA,BP_pvalue=NA,MF_term=NA,MF_pvalue=NA,CC_term=NA,CC_pvalue=NA)                                                                                   
        for (i in 1:nrow(myInterestingGenes)){
        geneNames <- names(anoGO)
        InterestingGenes <- myInterestingGenes[i,]
        names(InterestingGenes) <-InterestingGenes
        geneList <- factor(as.integer(geneNames %in% InterestingGenes))        ###extract the GO of my interesting genes, labeled as 1.
        names(geneList) <- geneNames
#####use biological process to build topGO data##############
        GOBP <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)   #####use biological process to build topGO data
        GOMF <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)
        GOCC <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)
######perform GO annotation test############
        BPfis <- runTest(GOBP, algorithm = "weight", statistic = "fisher")               #####perform GO annotation test
        MFfis <- runTest(GOMF, algorithm = "weight", statistic = "fisher")
        CCfis <- runTest(GOCC, algorithm = "weight", statistic = "fisher")
#######analyze the GO annotation result#####################
        BPtable <- GenTable(GOBP, pvalue = BPfis, topNodes = 40)
        MFtable <- GenTable(GOMF, pvalue =MFfis, topNodes = 40)
        CCtable <- GenTable(GOCC, pvalue = CCfis, topNodes = 40)
        mat1 <- matrix(c(BPtable[1,2],BPtable[1,6],MFtable[1,2],MFtable[1,6],CCtable[1,2],CCtable[1,6]),nrow=1)
        mat[i,] <- mat1
        }
        rownames(mat) <- rownames(myInterestingGenes)
write.table(mat, filename,append=T, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
        return(mat)
}
#myInterestingGenes.topgo <- ano(lncrnaInterestingGenes,anoGO.new,"wgcna.correlation.topGO.txt")
#selected.InterestingGenes.topgo <- ano(selected.InterestingGenes,anoGO.new,"wgcna.correlation.topGO.selected.txt")

ano.parallel <- function(x){
library(topGO)
library(ALL)
	anoGO <- readMappings("/psc/bioinformatics/sunyd/identify/cucumber_rnaseq/SRRadd/cuc.new.map")
        geneNames <- names(anoGO)
        InterestingGenes <- x
        names(InterestingGenes) <-InterestingGenes
        geneList <- factor(as.integer(geneNames %in% InterestingGenes))        ###extract the GO of my interesting genes, labeled as 1.
        names(geneList) <- geneNames
#####use biological process to build topGO data##############
        GOBP <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)   #####use biological process to build topGO data
        GOMF <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)
        GOCC <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)
######perform GO annotation test############
        BPfis <- runTest(GOBP, algorithm = "weight", statistic = "fisher")               #####perform GO annotation test
        MFfis <- runTest(GOMF, algorithm = "weight", statistic = "fisher")
        CCfis <- runTest(GOCC, algorithm = "weight", statistic = "fisher")
#######analyze the GO annotation result#####################
        BPtable <- GenTable(GOBP, pvalue = BPfis, topNodes = 40)
        MFtable <- GenTable(GOMF, pvalue =MFfis, topNodes = 40)
        CCtable <- GenTable(GOCC, pvalue = CCfis, topNodes = 40)
        mat1 <- matrix(c(BPtable[1,2],BPtable[1,6],MFtable[1,2],MFtable[1,6],CCtable[1,2],CCtable[1,6]),nrow=1)
	mat1
}
cl <- makeCluster(getOption("cl.cores", 32))
system.time(lncrna.go <- parLapply(cl, lncrnaInterestingGenes.list, ano.parallel))
lncrna.go <- data.frame(t(sapply(lncrna.go,c)))
write.table(lncrna.go, "wgcna.correlation.topGO.txt",append=T, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

