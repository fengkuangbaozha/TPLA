
#R --vanilla < ~/tmp/R/scripts/topgo4arab_Genelist.R 
#R --vanilla < ~/tmp/R/scripts/topgo4m6ack.R   

library(topGO)
library(ALL)
#library(RColorBrewer)
anoGO.new <- readMappings("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuc.new.map")
anoGO.old <- readMappings("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuc.map")
str(head(anoGO.new))
str(head(anoGO.old))
#setwd("/home/liyq/m6A/")
#input= "flowerck_vs_leafck.diff.1"
#input = commandArgs()[3]
#base = input
#base = sub('.txt$', '', input)
#de <- read.table(file=input,header=F,row.names=1)

##################performing topGO analysis################
ano <- function(myInterestingGenes,anoGO,filename){
geneNames <- names(anoGO)
InterestingGenes <- myInterestingGenes[,1]
names(InterestingGenes) <-InterestingGenes
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
BPtable <- GenTable(GOBP, pvalue = BPfis, topNodes = 40)
if(sum(BPtable$pvalue =="< 1e-30")>0){BPtable[BPtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
BPtable <- BPtable[as.numeric(BPtable$pvalue) < 0.01 ,]
BPtable$Ontology="BP"
a <- BPtable[order(as.numeric(BPtable$pvalue),decreasing=T),]
MFtable <- GenTable(GOMF, pvalue =MFfis, topNodes = 40)
if(sum(MFtable$pvalue =="< 1e-30")>0){MFtable[MFtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
MFtable <- MFtable[as.numeric(MFtable$pvalue) < 0.01,]
MFtable$Ontology="MF"
b <- MFtable[order(as.numeric(MFtable$pvalue),decreasing=T),]
CCtable <- GenTable(GOCC, pvalue = CCfis, topNodes = 40)
if(sum(CCtable$pvalue =="< 1e-30")>0){CCtable[CCtable$pvalue =="< 1e-30",]$pvalue="1e-30"}
CCtable <- CCtable[as.numeric(CCtable$pvalue) < 0.05,]
CCtable$Ontology="CC"
c <- CCtable[order(as.numeric(CCtable$pvalue),decreasing=T),]

o <- rbind(a,b,c) 
merge(merge(BPtable,MFtable,all=T),CCtable,all=T) -> up.table
up.table[order(up.table[,7]),]->up.table

######output the gene list of sig terms
#scoresInTerm(GOBP,c[,1],use.names=T) ->bp1
#scoresInTerm(GOMF,b[,1],use.names=T) ->bp2
#scoresInTerm(GOCC,a[,1],use.names=T) ->bp3
##rbindlist(bp1,bp2,bp3) ->bp
#lapply(bp1,function(x) return(names(x[x==2]))) -> mylist1
#lapply(bp2,function(x) return(names(x[x==2]))) -> mylist2
#lapply(bp3,function(x) return(names(x[x==2]))) -> mylist3

#write(sapply(names(mylist1),function(x) paste(x,paste(mylist1[[x]],collapse="//"))),file=paste(filename,".upgo.genelist",sep=""),append=T)
#write(sapply(names(mylist2),function(x) paste(x,paste(mylist2[[x]],collapse="//"))),file=paste(filename,".upgo.genelist",sep=""),append=T)
#write(sapply(names(mylist3),function(x) paste(x,paste(mylist3[[x]],collapse="//"))),file=paste(filename,".upgo.genelist",sep=""),append=T)

write.table(file=paste(filename,".upgo",sep=""),o,quote=F,sep="\t",row.names=F)

#write.table(file=paste(filename,".downgo.txt",sep=""),down.table,quote=F,sep="\t",row.names=F)
#write.table(file=paste(filename,".upgo.txt",sep=""),up.table,quote=F,sep="\t",row.names=F)

#write.table(file=paste(filename,".downgo",sep=""),t,quote=F,sep="\t",row.names=F)

#system("perl ~/sunyd/identify/script/print_topgo.pl /home/liyq/Public/earDGE/GO.terms_and_ids *go") 

#pdf(paste(filename,"plotGoBarplot.pdf",sep=""),height=14,width=13)
png(paste(filename,"plotGoBarplot.png",sep=""))
 #library(RColorBrewer)
#layout(matrix(1,1), widths=lcm(7), heights=lcm(15))
# BASE GO TERMS PROPORTIONS
par(mar=c(1,22,3,1))
t<-read.delim(file=paste(filename,".upgo",sep=""),header=T,sep="\t")

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
barplot(-(log10(pvalue)),width=1,space=1,cex.names=.6,las=1,horiz=T,col=t$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(pvalue)))+1))
legend("right",legend=c("GOBP","GOMF","GOCC"),fill=c("red","blue","green"),bty="n",cex=.6,border=F)
axis(3,cex.axis=.6,mgp=c(0.5,0.3,-0.6),tck=-0.02)
mtext(expression(paste("-Log10(",italic(P),"value)")),side=3,line=1,cex=1)
dev.off()
return(pvalue)
}

brown.InterestingGenes <- read.delim("wgcna.kmean.phloem.spec.module_brown.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, both hclus and kmeans
yellow.InterestingGenes <- read.delim("wgcna.kmean.phloem.spec.module_yellow.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans
turquoise.InterestingGenes <- read.delim("wgcna.kmean.gua.spec.module_turquoise.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans not good
greenyellow.InterestingGenes <- read.delim("wgcna.kmean.gua.spec.module_greenyellow.pedicle.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans not good
royalblue.InterestingGenes <- read.delim("wgcna.kmean.gua.spec.module_royalblue.fruit.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans not good
blue.InterestingGenes <- read.delim("wgcna.kmean.flower.spec.module_blue.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans not good
green.InterestingGenes <- read.delim("wgcna.kmean.flower.spec.module_green.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans not good
greenyellow.InterestingGenes <- read.delim("wgcna.kmean.flower.spec.module_greenyellow.genename",header=F,stringsAsFactors = FALSE)  ####read my interesting genes, kmeans not good

ano(brown.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/phloem.coexpress.new.brown")
ano(yellow.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/phloem.coexpress.new.yellow")
ano(turquoise.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/gua.coexpress.new.turquoise")
ano(greenyellow.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/gua.coexpress.new.greenyellow")
ano(royalblue.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/gua.coexpress.new.royalblue")
ano(green.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/flower.coexpress.new.green")
ano(blue.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/flower.coexpress.new.blue")
ano(greenyellow.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/flower.coexpress.new.greenyellow")

ano(brown.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/phloem.coexpress.old.brown")
ano(yellow.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/phloem.coexpress.old.yellow")
ano(turquoise.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/gua.coexpress.old.turquoise")
ano(greenyellow.InterestingGenes,anoGO.old,"~/sunyd/identify/cucumber_rnaseq/SRRadd/gua.coexpress.old.greenyellow")
ano(royalblue.InterestingGenes,anoGO.old,"~/sunyd/identify/cucumber_rnaseq/SRRadd/gua.coexpress.old.royalblue")
ano(green.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/flower.coexpress.old.green")
ano(blue.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/flower.coexpress.old.blue")
ano(greenyellow.InterestingGenes,anoGO.new,"~/sunyd/identify/cucumber_rnaseq/SRRadd/flower.coexpress.old.greenyellow")
#brewer.pal(9,"Set1") ->col
#t$col[t$Ontology=="BP"] <- col[3]
#t$col[t$Ontology=="MF"] <-col[2] 
#t$col[t$Ontology=="CC"] <- col[1]
#legend("right",legend=c("GOBP","GOMF","GOCC"),fill=col[c(3,2,1)],bty="n",cex=.6,border=F)
#mtext("Enriched Go terms of differential expressed genes",side=3,line=2,cex=1.1,adj=0)

png(paste(filename,"plotGoBarplot.png",sep=""))
dev.off()
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
