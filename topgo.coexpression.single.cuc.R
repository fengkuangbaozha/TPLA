args <- commandArgs(trailingOnly = TRUE)

                ###########make GO annotation for lncrna coexpressed gene###########
library(topGO)
library(ALL)
library(parallel)
anoGO <- readMappings(file=args[1])
filename <- args[3]
lncrnaInterestingGenes <- as.matrix(read.delim(file=args[2],header=F,row.names=1,stringsAsFactors=FALSE))[,1:20]
lncrnaInterestingGenes.list <- split(lncrnaInterestingGenes,rownames(lncrnaInterestingGenes))  ######convert to a list to perform lapply
ano.parallel <- function(x){
library(topGO)
library(ALL)
	anoGO <- readMappings("/psc/bioinformatics/sunyd/genome/Cucumber/cuc.new.map")
        geneNames <- names(anoGO)
        InterestingGenes <- x
        names(InterestingGenes) <-InterestingGenes
        print(class(InterestingGenes))
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
write.table(lncrna.go, paste(filename,".correlation.topGO.txt",sep=""),append=F, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

