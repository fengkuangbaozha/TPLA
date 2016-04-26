args <- commandArgs(trailingOnly = TRUE)
up.down <- read.delim(file=args[1],row.names=1,header=F)
filename <- args[2]
#up.down <- read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/cuffdiff.028766_root_tmp/root.up.down.bar",row.names=1,header=F)
#filename <- "~/sunyd/identify/oryza_rnaseq/SRRpi/cuffdiff.028766_root_tmp/root.up.down.bar"

up.down <- up.down[c(1,11,9,10,12,2,3,4,7,8,5,6),]
t.up.down <- t(up.down)
num <- ceiling(max(c(t.up.down[1,],t.up.down[2,])/100)) * 100  ####caculate the max num of the matrix and set the border value
#library(RColorBrewer)
#col <- brewer.pal(9,"Paired")
#col <-"gray"
pdf(paste(filename,".pdf",sep=""))
COLA="#ab2023"
COLB="#1c76b1"
######plot the up part plot######
bar <- barplot(t.up.down[1,],ylim=c(-num,num),space=0.9,width=0.5,main="Number of DE-lncRNAs under root Pi starvation",axisnames=FALSE, col=COLA,border="grey",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
#mtext("Up regulation",side=2, line=2, adj=0.8, col=COLA, cex=1.2)
axis(1,at=bar,labels=F)
text(bar,t.up.down[1,]+10,t.up.down[1,],cex=1.5)
text(x=bar+.25,y =bar-425,labels=colnames(t.up.down),cex=1.5,srt=45,pos=2,xpd=TRUE)
box(which="plot",lty="solid",col="black")
#######plot the down part plot#####
barplot(-t.up.down[2,], add=T,col=COLB,border="grey",space=0.9,width=0.5,axisnames=F,cex.axis=1.5,cex.lab=1.5)
#mtext("Down regulation",side=2, line=2, adj=0.2, col=COLB, cex=1.2)
text(bar,-t.up.down[2,]-10,t.up.down[2,],cex=1.5)
#text(bar,-(num),rownames(up.down),cex=0.9)
abline(h=0,col="black")
legend("topleft",border=F, c("Up regulation", "Down regulation"),pch=15, col=c(COLA,COLB),bty="n",pt.cex=1,cex=1.5)
dev.off()
