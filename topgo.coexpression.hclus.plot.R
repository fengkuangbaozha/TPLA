args <- commandArgs(trailingOnly = TRUE)
o<-read.delim(file=args[1],header=T)
filename <- args[2]
#o<-read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/coexpression/root.de.lncrna20.upgo.1",header=T)
library(RColorBrewer)
o$Ontology = as.character(o$Ontology)
o$Ontology[o$Ontology=="CC"] <- "OC"
t <- o[order(o[,7],o[,6],decreasing=T),]

go<-t$Term
pvalue<- t$pvalue 
names(pvalue) <-go
func<-t$Ontology
as.vector(t$Ontology) -> t$Ontology

brewer.pal(9,"Set1") ->col
t$col[t$Ontology=="BP"] <- col[1]
t$col[t$Ontology=="MF"] <-col[2] 
t$col[t$Ontology=="OC"] <- col[3]
pdf(paste(filename,"plotGoBarplot.pdf",sep=""))
par(mar=c(1,22,3,1))
barplot(-(log10(pvalue)),width=1,space=1,cex.names=1.2,las=1,horiz=T,col=t$col,axes=F,bg="white",border=F,xlim=c(0,max(-(log10(pvalue)))+1),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
legend("bottomright",legend=c("GOBP","GOMF","GOCC"),fill=col[c(1,2,3)],bty="n",pt.cex=1,cex=1.2,border=F)
axis(3,cex.axis=1.2,mgp=c(0.5,0.3,-0.6),tck=-0.02)
#axis(3,cex.axis=1.5,tck=-0.02)
mtext(expression(paste("-Log10(",italic(P),"value)")),side=3,line=1,cex=1)
dev.off()
