
#R --vanilla < ~/tmp/R/scripts/topgo4arab_Genelist.R 
#R --vanilla < ~/tmp/R/scripts/topgo4m6ack.R   
args <- commandArgs(trailingOnly = TRUE)

library(topGO)
library(ALL)
#library(RColorBrewer)
anoGO <- readMappings(file=args[1])
str(head(anoGO))
brown.InterestingGenes <- read.delim(file=args[2],header=F,stringsAsFactors = FALSE)  ####read my interesting genes, both hclus and kmeans
filename <- args[3]
topnode <- args[4]
num1 <- args[5]
num2 <- args[6]
num3 <- args[7]
#anoGO <- readMappings("/psc/bioinformatics/sunyd/genome/tomato/ITAG2.4.go.csv")
#setwd("/home/liyq/m6A/")
#input= "flowerck_vs_leafck.diff.1"
#input = commandArgs()[3]
#base = input
#base = sub('.txt$', '', input)
#de <- read.table(file=input,header=F,row.names=1)

##################performing topGO analysis################
ano <- function(myInterestingGenes,anoGO,filename,num1,num2,num3){
geneNames <- names(anoGO)
InterestingGenes <- myInterestingGenes[,1]
names(InterestingGenes) <- InterestingGenes
geneList <- factor(as.integer(geneNames %in% InterestingGenes))        ###extract the GO of my interesting genes, labeled as 1.
names(geneList) <- geneNames
#str(geneList)

#####use biological process to build topGO data##############
GOBP <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)   #####use biological process to build topGO data
GOMF <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)
GOCC <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = anoGO)

######perform GO annotation test############
BPfis <- runTest(GOBP, algorithm = "weight", statistic = "fisher")               #####perform GO annotation test
MFfis <- runTest(GOMF, algorithm = "weight", statistic = "fisher")
CCfis <- runTest(GOCC, algorithm = "weight", statistic = "fisher")

#######analyze the GO annotation result#####################
BPtable <- GenTable(GOBP, pvalue = BPfis, topNodes = topnode)
if(sum(BPtable$pvalue =="< 1e-30")>0){BPtable[BPtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
BPtable <- BPtable[as.numeric(BPtable$pvalue) < num1 ,]
BPtable$Ontology="BP"
a <- BPtable[order(as.numeric(BPtable$pvalue),decreasing=T),]
MFtable <- GenTable(GOMF, pvalue =MFfis, topNodes = topnode)
if(sum(MFtable$pvalue =="< 1e-30")>0){MFtable[MFtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
MFtable <- MFtable[as.numeric(MFtable$pvalue) < num2,]
MFtable$Ontology="MF"
b <- MFtable[order(as.numeric(MFtable$pvalue),decreasing=T),]
CCtable <- GenTable(GOCC, pvalue = CCfis, topNodes = topnode)
if(sum(CCtable$pvalue =="< 1e-30")>0){CCtable[CCtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
CCtable <- CCtable[as.numeric(CCtable$pvalue) < num3,]
CCtable$Ontology="CC"
c <- CCtable[order(as.numeric(CCtable$pvalue),decreasing=T),]

o <- rbind(a,b,c) 
up.table <- merge(merge(BPtable,MFtable,all=T),CCtable,all=T)
up.table <- up.table[order(up.table$Ontology),]

write.table(file=paste(filename,".upgo",sep=""),o,quote=F,sep="\t",row.names=F)     #####the significant GO annontation 
#system("perl ~/sunyd/identify/script/print_topgo.pl ~/sunyd/identify/script/GO.terms_and_ids *.upgo")
scoresInTerm(GOBP,a[,1],use.names=T) ->bp1
scoresInTerm(GOMF,b[,1],use.names=T) ->bp2
scoresInTerm(GOCC,c[,1],use.names=T) ->bp3
#rbindlist(bp1,bp2,bp3) ->bp
lapply(bp1,function(x) return(names(x[x==2]))) -> mylist1
lapply(bp2,function(x) return(names(x[x==2]))) -> mylist2
lapply(bp3,function(x) return(names(x[x==2]))) -> mylist3

write(sapply(names(mylist1),function(x) paste(x,paste(mylist1[[x]],collapse="//"))),file=paste(filename,"upgo.genelist",sep=""),append=T)
write(sapply(names(mylist2),function(x) paste(x,paste(mylist2[[x]],collapse="//"))),file=paste(filename,"upgo.genelist",sep=""),append=T)
write(sapply(names(mylist3),function(x) paste(x,paste(mylist3[[x]],collapse="//"))),file=paste(filename,"upgo.genelist",sep=""),append=T)

png(paste(filename,"plotGoBarplot.png",sep=""))
library(RColorBrewer)
#layout(matrix(1,1), widths=lcm(7), heights=lcm(15))
# BASE GO TERMS PROPORTIONS
par(mar=c(1,18,3,1))
#t<-read.delim(file=paste(filename,".upgo",sep=""),header=T,sep="\t")
t <- o
go<-t$Term
pvalue<- as.numeric(t$pvalue)
names(pvalue) <-go
func<-t$Ontology
as.vector(t$Ontology) -> t$Ontology

brewer.pal(9,"Set1") ->col
t$col[t$Ontology=="BP"] <- col[3]
t$col[t$Ontology=="MF"] <-col[2] 
t$col[t$Ontology=="CC"] <- col[1]
#enrich<-c(go,-log10(pvalue))
#barplot(-(log10(pvalue)),cex.names=.6,las=1,horiz=T,col=t$col,axes=T,border=F,xlab =expression(paste("-Log10(",italic(P),"value)")),xlim=c(0,max(-(log10(pvalue)))+1))
barplot(-(log10(pvalue)),width=1,space=1,cex.names=1,cex.lab=1,las=1,horiz=T,col=t$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(pvalue)))+1))
legend("right",legend=c("GOBP","GOMF","GOCC"),fill=col[c(3,2,1)],bty="n",cex=.8,border=F)
axis(3,cex.axis=1,mgp=c(0.5,0.3,-0.6),tck=-0.02)
#axis(3,cex.axis=.6,tck=-0.02)
mtext(expression(paste("-Log10(",italic(P),"value)")),side=3,line=1,cex=1.5)
dev.off()
return(pvalue)
}


ano(brown.InterestingGenes,anoGO,filename,num1,num2,num3)
#legend("right",legend=c("GOBP","GOMF","GOCC"),fill=col[c(3,2,1)],bty="n",cex=.6,border=F)
#mtext("Enriched Go terms of differential expressed genes",side=3,line=2,cex=1.1,adj=0)

#png(paste(filename,"plotGoBarplot.png",sep=""))
#dev.off()
plotgo <- function(t,num,filename){
par(mar=c(1,22,3,1))
go<-t$Term
pvalue<- as.numeric(t$pvalue)
names(pvalue) <-go
func<-t$Ontology
as.vector(t$Ontology) -> t$Ontology

t$col[t$Ontology=="BP"] <- "red"
t$col[t$Ontology=="MF"] <-"blue" 
t$col[t$Ontology=="CC"] <- "green"
#enrich<-c(go,-log10(pvalue))
#barplot(-(log10(pvalue)),cex.names=.6,las=1,horiz=T,col=t$col,axes=T,border=F,xlab =expression(paste("-Log10(",italic(P),"value)")),xlim=c(0,max(-(log10(pvalue)))+1))
barplot(-(log10(pvalue)),width=1,space=1,cex.names=num,las=1,horiz=T,col=t$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(pvalue)))+1))
legend("right",legend=c("GOBP","GOMF","GOCC"),fill=c("red","blue","green"),bty="n",cex=.6,border=F)
axis(3,cex.axis=.6,mgp=c(0.5,0.3,-0.6),tck=-0.02)
}
