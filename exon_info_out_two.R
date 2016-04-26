args <- commandArgs(trailingOnly = TRUE)
#combinedgtf <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75.combined.gtf",header=F)
gtf <- read.delim(file=args[1],header=F,stringsAsFactors = FALSE) 
#gtf <- gtf[-which(gtf[,4]>10000),]
colnames(gtf) <- c("gene_id","transcript_id","class_code","exonlength","start","end","strand","chr_loc")
seqcomp <- read.delim(file=args[2],header=F,stringsAsFactors = FALSE,row.names=1)
genename <- read.delim(file=args[3],header=F,stringsAsFactors = FALSE) 
genename <- as.matrix(genename[,2])
lncname <- read.delim(file=args[4],header=F,stringsAsFactors = FALSE)
lncname <- as.matrix(lncname[,2])
filename <- args[5]
num1 <- as.numeric(args[6])
num2 <- as.numeric(args[7])
COLA="#ab2023"
COLB="#1c76b1"
legendname3 <- c("lncRNA","PCgene")
color3 <- c(COLA,COLB)

##########count the exon number################
exoncount <- function(x){    #######exon count number for all transcripts(isoforms)######
	exonname <- unique(x[,1:2])   ##extract the name of geneid and transcriptid
	exon <- as.matrix(x[,2])   ###for genename input
	count <- cbind(exonname,as.matrix(table(exon)))
	count <- as.data.frame(count,stringsAsFactors=F)
	return(count)
	}
#combinedexon <- exoncount(combinedgenename)   #####exon counts for the all genes
exon <- data.matrix(table(as.matrix(gtf[,2])))
geneexon <- exon[which(rownames(exon) %in% genename),,drop=FALSE] ###combinedexon[which(as.character(combinedexon$transcript_id) %in% genename$transcript_id),]
lncexon <- exon[which(rownames(exon) %in% lncname),,drop=FALSE] ###combinedexon[which(as.character(combinedexon$transcript_id) %in% genename$transcript_id),]
lncexon2 <- lncexon[-which(lncexon == 1),]
lncexon2all <- gtf[which(gtf[,2] %in% as.matrix(names(lncexon2))),,drop=FALSE]
for (i in 1:nrow(lncexon2all)){
    if (identical(lncexon2all[i,6],lncexon2all[i+1,5])){
              print(lncexon2all[i,2])
# }else{
 #             print(" ")
		}
	}
write.table(lncexon2all,paste(filename,"lnc_intron.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
##########calculate the cumulative number class of each trans exon num for barplot##########
barinput <- function(x,num){              ######the num differs depends on the data you operate
	tmp <- data.matrix(table(x[,1]))      ###how many trans contain 2 exon and so on
	tmp <- (c(tmp[1:num,],sum(tmp[(num+1):nrow(tmp),])) / sum(tmp)) * 100     ####display 1 to 15 and over 15
	return(tmp)
	}
lncbar <- barinput(lncexon,num1)
lncbar <- c(lncbar,rep(0,6-num1))
genebar <- barinput(geneexon,6)
bar <- rbind(lncbar,genebar)
colnames(bar) <- NULL
png(paste(filename,"lnc_exonnum.png",sep=""))
mp = barplot(bar,beside=T,main="Exons per transcript",ylab="Proportion", col=color3,ylim=c(0,100),cex.axis=1.5,cex.lab=1.5,cex.main=2,cex.sub=1.5)
box(which="plot",lty="solid",col="black")
colnames(bar) <- c(1:6,c("> 6"))
axis(1,at=seq(2,20,by=3),labels=colnames(bar),cex.axis=1.5)
legend("topright",legend=legendname3,col=color3,bg="white",lwd=2,pt.cex=1,cex=1.5)
dev.off()

#############the exon length distribution##########
exonlen <- as.matrix(gtf[,4])      ###########combinedgenename[which(as.character(combinedgenename$transcript_id) %in% genename$transcript_id),]
rownames(exonlen) <- gtf[,2]
geneexon.len <- exonlen[which(rownames(exonlen) %in% genename),,drop=FALSE]      ###
lncexon.len <- exonlen[which(rownames(exonlen) %in% lncname),,drop=FALSE]      ###
#iexonlen <- as.matrix(lncname[which(lncname[,2] %in% iclass),4])         ###combinedgenename[which(as.character(combinedgenename$transcript_id) %in% lncname$transcript_id),]

png(paste(filename,"lnc_exonlen.png",sep=""))
plot(density(geneexon.len),type="l", xlab = "Exon length",ylab = "Density",col=COLB,main="Exon length distribution",xlim=c(0,5000),lty=1,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
lines(density(lncexon.len),type="l",col=COLA,lty=1,lwd=2)
legend("topright",legend=legendname3,col=color3,bg="white",lwd=2,pt.cex=1,cex=1.5)
dev.off()

plot2densityline <- function(input,inputtwo,xname,mainname,filename,legendname){
pdf(filename)
plot(density(input),type="l", xlab = xname,ylab = "Density",col=COLB,main=mainname,lty=1,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
lines(density(inputtwo),type="l",col=COLA,lty=1,lwd=2)
legend("topright",legend=legendname,col=color3,bg="white",lwd=2,pt.cex=1,cex=1.5)
dev.off()
}
#plot2densityline(geneexonlen,lncexonlen,"exon length","Exon length distribution","lnc_exonlen.png",c("Long noncoding RNAs","mRNA"))

#############the sequence composition information###############
#genename <- as.matrix(unique(gtf[,2]))                  ###############codeext(genename)
#lncname <- as.matrix(unique(lncname[,2]))                    ############codeext(lncname)

##########sequence length distribution plot##########
genelen <- seqcomp[genename,1,drop=FALSE] 
lnclen <- seqcomp[lncname,1,drop=FALSE]
valu <- c(data.matrix(lnclen),data.matrix(genelen))
tim <- c(nrow(lnclen),nrow(genelen))
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
png(paste(filename,"lnc_length.png",sep=""))
boxplot(values ~ vars, data = df, col = color3, main="Transcript length distribution",cex.axis=1.5,cex.lab=1.5,cex.main=2)
dev.off()
png(paste(filename,"lnc_length_only.png",sep=""))
hist(data.matrix(lnclen),main="Histogram of lncRNAs length",xlab="Length",breaks=20,freq=FALSE)
dev.off()

print(summary(genelen))
print(summary(lnclen))
#valu <- c(data.matrix(genelen),data.matrix(lnclen),data.matrix(lnclen[iclass,,drop=FALSE]),data.matrix(lnclen[uclass,,drop=FALSE]),data.matrix(lnclen[oclass,,drop=FALSE]),data.matrix(lnclen[xclass,,drop=FALSE]),data.matrix(lnclen[jclass,,drop=FALSE]),data.matrix(lnclen[eclass,,drop=FALSE]))
#tim <- c(nrow(genelen),nrow(lnclen),nrow(lnclen[iclass,,drop=FALSE]),nrow(lnclen[uclass,,drop=FALSE]),nrow(lnclen[oclass,,drop=FALSE]),nrow(lnclen[xclass,,drop=FALSE]),nrow(lnclen[jclass,,drop=FALSE]),nrow(lnclen[eclass,,drop=FALSE]))

#plot2densityline(genelen,lnclen,"Length","Transcript length distribution","lnc_length.png",c("mRNA","Long noncoding RNAs"))
#lnclen[which(lnclen>=400 & lnclen<=1000),]

#########AU,GC content plot#############
#######calculate the cumsum for cumulative curve########
cumplot <- function(x){
	breaks <- seq(0,1,by=0.01)
	xcut <- cut(x,breaks,right=FALSE)
	xfreq <- table(xcut)
	xcum <- c(0,cumsum(xfreq)) / length(x)
	mat <- cbind(breaks,xcum)
	return(mat)
	}

genegccum <- cumplot(seqcomp[genename,2])
lncgccum <- cumplot(seqcomp[lncname,2])

png(paste(filename,"lnc_gc.png",sep=""))
plot(genegccum[,1],genegccum[,2],type="l",col=COLB,xlab="GC content",ylab="Cumulative frequency",main="The GC content of transcripts",lty=1,lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=2)
lines(lncgccum[,1],lncgccum[,2],type="l",col=COLA,lty=1,lwd=2)
legend("bottomright",legend=legendname3,col=color3,bg="white",lwd=2,pt.cex=1,cex=1.5)
dev.off()



contentplot <- function(col,xname,mainname,filename){
	lncgc <- seqcomp[lncname,col]
	genegc <- seqcomp[genename,col]
	lncgccum <- cumplot(lncgc)
	genegccum <- cumplot(genegc)
png(filename)
	plot(genegccum[,1],genegccum[,2],type="l",col=COLB,xlab=xname,ylab="Cumulative frequency",main=mainname)
	lines(lncgccum[,1],lncgccum[,2],type="l",col=COLA)
	dev.off()
	}
