args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
info <- read.delim(file=input,header=T,stringsAsFactors = FALSE,row.names=1)
#info <- read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/coexpression/co2gene.TCONS_00131442",header=T,stringsAsFactors = FALSE,row.names=1)
#filename <- args[2]
#info_root <- info[,c(13,20,18,19,21,14,15,17,16,4,11,9,10,12,5,6,8,7,1,3,2)]
#info_shoot <- info[,c(34,41,39,40,42,35,36,38,37,25,32,30,31,33,26,27,29,28,22,24,23)]
info_root <- info[,c(13,4,20,11,18,9,19,10,21,12,14,5,15,6,1,17,8,3,16,7,2)]
info_shoot <- info[,c(34,25,41,32,39,30,40,31,42,33,35,26,36,27,22,38,29,24,37,28,23)]
COLA="#ab2023"
COLB="#1c76b1"
num=ncol(info)/2
pdf(paste(input,"_coexpression.pdf",sep=""))
plot(x=c(1:num),y=info_root[1,],type="l",ylim=c(0,max(info_root[1:2,])+2),xaxt='n',col=COLA,xlab="",ylab="log2FPKM",main="RNA-seq",lty=1,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
lines(x=c(1:num),y=info_root[2,],type="l",col=COLB,lty=1,lwd=2)
##lines(x=c(1:num),y=info_shoot[1,],type="l",col=COLA,lty=5,lwd=2)
##lines(x=c(1:num),y=info_shoot[2,],type="l",col=COLB,lty=5,lwd=2)
axis(side=1,at=seq(1,21,by=1),labels = FALSE,tck=-0.01)
##axis(side=1, at=seq(1,21,by=1), labels=colnames(info_root), tck=-0.00,cex.axis=1)
#text(x=seq(1.5,21.5,by=1),y =-0.55,labels=colnames(info_root),cex=1.2,srt=45,pos=2,xpd=TRUE)
legend("topright",legend=rownames(info),col=c(COLA,COLB),bg="white",lwd=2,pt.cex=1,cex=1.2)
dev.off()
