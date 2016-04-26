args <- commandArgs(trailingOnly = TRUE)
actin <- read.delim(file=args[1],header=F,stringsAsFactors = FALSE,row.names=1)
gene <- read.delim(file=args[2],header=F,stringsAsFactors = FALSE,row.names=1)
genename = args[3]
lnc <- read.delim(file=args[4],header=F,stringsAsFactors = FALSE,row.names=1)
lncname = args[5]
filename <- args[6]
dim(actin)
dim(gene)
#print(gene)
dim(lnc)
#info <- read.delim(file="~/sunyd/identify/oryza_rnaseq/SRRpi/qrt-pcr/result0808",header=F,stringsAsFactors = FALSE,row.names=1)
#gene = info[,c(1:3)]
#lnc = info[,c((colnum-1):6)]
#actin = info[,tail(colnames(info),3)]
sig <- function(x){
	x[which(x <= 0.001),] = "**"
	x[which(x <= 0.05),] = "*"
	x[which(x > 0.05),] = NA
	return(x)
}
order = c(13,1,14,2,15,3,16,4,17,5,18,6,19,7,10,20,8,11,21,9,12)
library("Hmisc")
COLA="#ab2023"
COLB="#1c76b1"
control = c(13:21)
control21d = c(19:21)
gene_delta = (gene - actin)
gene_delta_control = rbind(gene_delta[control,],gene_delta[control21d,],gene_delta[control,])
gene_delta_delta = 2^(-(gene_delta - gene_delta_control))
gene_delta_delta = cbind(gene_delta_delta,mean=apply(gene_delta_delta,1,mean),sd=apply(gene_delta_delta,1,sd))
gene_delta_delta_control = rbind(gene_delta_delta[control,],gene_delta_delta[control21d,],gene_delta_delta[control,])
gene_delta_delta_test = cbind(gene_delta_delta[,1:3],gene_delta_delta_control[,1:3])
gene_test = as.matrix(apply(gene_delta_delta_test[1:9,],1,function(x){t.test(x[4:6],x[1:3])$p.value}))
gene_test = sig(gene_test)
lnc_delta = (lnc - actin)
lnc_delta_control = rbind(lnc_delta[control,],lnc_delta[control21d,],lnc_delta[control,])
lnc_delta_delta = 2^(-(lnc_delta - lnc_delta_control))
lnc_delta_delta = cbind(lnc_delta_delta,mean=apply(lnc_delta_delta,1,mean),sd=apply(lnc_delta_delta,1,sd))
lnc_delta_delta_control = rbind(lnc_delta_delta[control,],lnc_delta_delta[control21d,],lnc_delta_delta[control,])
lnc_delta_delta_test = cbind(lnc_delta_delta[,1:3],lnc_delta_delta_control[,1:3])
lnc_test = as.matrix(apply(lnc_delta_delta_test[1:9,],1,function(x){t.test(x[4:6],x[1:3])$p.value}))
lnc_test = sig(lnc_test)

pdf(paste(filename,"_qrt.pdf",sep=""))
num=nrow(lnc_delta_delta)
colnum=ncol(lnc_delta_delta)
signum=c(2,4,6,8,10,12,14,17,20)
plot(x=c(1:num),y=lnc_delta_delta[order,(colnum-1)],type="l",ylim=c(0,max(c(lnc_delta_delta[,(colnum-1)],gene_delta_delta[,(colnum-1)])) + 2),xaxt='n',col=COLA,xlab="",ylab="Fold change",main="qRT-PCR",lty=1,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
#plot(x=c(1:num),y=lnc_delta_delta[order,(colnum-1)],type="l",ylim=c(0,10),xaxt='n',col=COLA,xlab="",ylab="Fold change",main="qRT-PCR expression in various Pi samples",lty=1,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
lines(x=c(1:num),y=gene_delta_delta[order,(colnum-1)],type="l",col=COLB,lty=1,lwd=2)
text(x=signum-0.4,y =lnc_delta_delta[order,(colnum-1)][signum],labels=lnc_test,cex=2,col=COLA,xpd=TRUE)
text(x=signum+0.4,y =gene_delta_delta[order,(colnum-1)][signum],labels=gene_test,cex=2,col=COLB,xpd=TRUE)
axis(side=1,at=seq(1,num,by=1),labels = FALSE,tck=-0.01)
text(x=seq(1.5,num+0.5,by=1),y =-.3,labels=rownames(lnc_delta_delta)[order],cex=1.2,srt=45,pos=2,xpd=TRUE)
errbar(x=c(1:num),y=lnc_delta_delta[order,(colnum-1)],(lnc_delta_delta[order,(colnum-1)]+(1.96 * lnc_delta_delta[order,colnum] / 2)),lnc_delta_delta[order,(colnum-1)],add=T,pch=20)  ##(lnc_delta_delta[order,(colnum-1)]-(1.96 * lnc_delta_delta[order,colnum] / 2))
errbar(x=c(1:num),y=gene_delta_delta[order,(colnum-1)],(gene_delta_delta[order,(colnum-1)]+(1.96 * gene_delta_delta[order,colnum] / 2)),gene_delta_delta[order,(colnum-1)],add=T,pch=20)  ##(gene_delta_delta[order,(colnum-1)]-(1.96 * gene_delta_delta[order,colnum] / 2)),
legend("topright",legend=c(lncname,genename),col=c(COLA,COLB),bg="white",lwd=2,pt.cex=1,cex=1.2)
dev.off()









#m = t(cbind(gene_delta_delta[control,(colnum-1)],gene_delta_delta[c(1:9),(colnum-1)],lnc_delta_delta[control,(colnum-1)],lnc_delta_delta[c(1:9),(colnum-1)]))
#sd = t(cbind(gene_delta_delta[control,colnum],gene_delta_delta[c(1:9),colnum],lnc_delta_delta[control,colnum],lnc_delta_delta[c(1:9),colnum]))
#se = 1.96 * sd / 2
#axis(side=1, at=seq(1,21,by=1), labels=colnames(info_root), tck=-0.00,cex.axis=1)
#bar=barplot(m,beside=T,ylim=c(0,colnum),names.arg=LETTERS[1:9],col=c(COLB,COLA),space=rep(c(3,0,0,0),9))
#box(which="plot",lty="solid",col="black")
