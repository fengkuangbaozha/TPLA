
#R --vanilla < ~/tmp/R/scripts/topgo4arab_Genelist.R 
#R --vanilla < ~/tmp/R/scripts/topgo4m6ack.R   

library(topGO)
library(ALL)
library(RColorBrewer)
maizegeneID2GO <- readMappings(file = system.file("examples/cuc.map", package = "topGO"))
str(head(maizegeneID2GO))
#setwd("/home/liyq/m6A/")
#input= "flowerck_vs_leafck.diff.1"
input = commandArgs()[3]
#base = input
base = sub('.txt$', '', input)
de <- read.table(file=input,header=F,row.names=1)
geneNames <- names(maizegeneID2GO)
myInterestingGenes <- row.names(de)
names(myInterestingGenes) <-myInterestingGenes
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
#str(geneList)


GOBP <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = maizegeneID2GO)

GOMF <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = maizegeneID2GO)
GOCC <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = maizegeneID2GO)

BPfis <- runTest(GOBP, algorithm = "weight", statistic = "fisher")
MFfis <- runTest(GOMF, algorithm = "weight", statistic = "fisher")
CCfis <- runTest(GOCC, algorithm = "weight", statistic = "fisher")

BPtable <- GenTable(GOBP, pvalue = BPfis, topNodes = 40)
if(sum(BPtable$pvalue =="< 1e-30")>0){BPtable[BPtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
BPtable <- BPtable[as.numeric(BPtable$pvalue) < 0.05 ,]
BPtable$Ontology="BP"
MFtable <- GenTable(GOMF, pvalue =MFfis, topNodes = 40)
if(sum(MFtable$pvalue =="< 1e-30")>0){MFtable[MFtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
MFtable$Ontology="MF"
MFtable <- MFtable[as.numeric(MFtable$pvalue) < 0.01,]
CCtable <- GenTable(GOCC, pvalue = CCfis, topNodes = 40)
if(sum(CCtable$pvalue =="< 1e-30")>0){CCtable[CCtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
CCtable <- CCtable[as.numeric(CCtable$pvalue) < 0.05,]
CCtable$Ontology="CC"

a <- CCtable[order(as.numeric(CCtable$pvalue),decreasing=T),]
b <- MFtable[order(as.numeric(MFtable$pvalue),decreasing=T),]
c <- BPtable[order(as.numeric(BPtable$pvalue),decreasing=T),]
o <- rbind(a,b,c) 

merge(merge(BPtable,MFtable,all=T),CCtable,all=T) -> up.table
up.table[order(up.table[,7]),]->up.table
######output the gene list of sig terms
scoresInTerm(GOBP,c[,1],use.names=T) ->bp1
scoresInTerm(GOMF,b[,1],use.names=T) ->bp2
scoresInTerm(GOCC,a[,1],use.names=T) ->bp3
#rbindlist(bp1,bp2,bp3) ->bp
lapply(bp1,function(x) return(names(x[x==2]))) -> mylist1
lapply(bp2,function(x) return(names(x[x==2]))) -> mylist2
lapply(bp3,function(x) return(names(x[x==2]))) -> mylist3

write(sapply(names(mylist1),function(x) paste(x,paste(mylist1[[x]],collapse="//"))),file=paste(base,".upgo.genelist",sep=""),append=T)
write(sapply(names(mylist2),function(x) paste(x,paste(mylist2[[x]],collapse="//"))),file=paste(base,".upgo.genelist",sep=""),append=T)
write(sapply(names(mylist3),function(x) paste(x,paste(mylist3[[x]],collapse="//"))),file=paste(base,".upgo.genelist",sep=""),append=T)


#write.table(file=paste(base,".downgo.txt",sep=""),down.table,quote=F,sep="\t",row.names=F)
write.table(file=paste(base,".upgo.txt",sep=""),up.table,quote=F,sep="\t",row.names=F)

#write.table(file=paste(base,".downgo",sep=""),t,quote=F,sep="\t",row.names=F)
write.table(file=paste(base,".upgo",sep=""),o,quote=F,sep="\t",row.names=F)

system("perl ~/pro/print_topgo.pl /home/liyq/Public/earDGE/GO.terms_and_ids *go") 

#pdf(paste(base,"plotGoBarplot.pdf",sep=""),height=14,width=13)
pdf(paste(base,"plotGoBarplot.pdf",sep=""))
 library(RColorBrewer)
#layout(matrix(1,1), widths=lcm(7), heights=lcm(15))
# BASE GO TERMS PROPORTIONS
par(mar=c(1,22,3,1))
t<-read.table(file=paste(base,".upgo.1",sep=""),header=T,sep="\t")


go<-t$Term
pvalue<- as.numeric(t$pvalue)
names(pvalue) <-go
func<-t$Ontology
as.vector(t$Ontology) -> t$Ontology


brewer.pal(9,"Set1") ->col
t$col[t$Ontology=="BP"] <- col[3]
t$col[t$Ontology=="MF"] <-col[2] 
t$col[t$Ontology=="CC"] <- col[1]



#t$col[t$Ontology=="BP"] <- "cadetblue3"
#t$col[t$Ontology=="MF"] <-"skyblue3" 
#t$col[t$Ontology=="CC"] <- "cyan4"
#enrich<-c(go,-log10(pvalue))
#barplot(-(log10(pvalue)),cex.names=.6,las=1,horiz=T,col=t$col,axes=T,border=F,xlab =expression(paste("-Log10(",italic(P),"value)")),xlim=c(0,max(-(log10(pvalue)))+1))
barplot(-(log10(pvalue)),cex.names=.6,las=1,horiz=T,col=t$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(pvalue)))+1))


legend("right",legend=c("GOBP","GOMF","GOCC"),fill=col[c(3,2,1)],bty="n",cex=.6,border=F)
axis(3,cex.axis=.6,mgp=c(0.5,0.3,-0.6),

tck=-0.02)

#mtext("Enriched Go terms of differential expressed genes",side=3,line=2,cex=1.1,adj=0)
mtext(expression(paste("-Log10(",italic(P),"value)")),side=3,line=1,cex=1)


dev.off()
