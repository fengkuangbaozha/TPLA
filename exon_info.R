#combinedgtf <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75.combined.gtf",header=F)
gtf <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRall/cuff226.expression.gene.seqinfo.gtf-exon",header=F,stringsAsFactors = FALSE) 
gtf <- gtf[-which(gtf[,4]>10000),]
colnames(gtf) <- c("gene_id","transcript_id","class_code","exonlength","start","end","strand","chr_loc")
seqcomp <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRall/cuff226.sequence_info.txt",header=F,stringsAsFactors = FALSE,row.names=1)
genename <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRall/cuff226.genename.ano",header=F,stringsAsFactors = FALSE) 
genename <- as.matrix(genename[,2])
lncname <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRall/cuff226-CPC_left.lncname.iuox.del",header=F,stringsAsFactors = FALSE)
lncnamedel <- read.delim("~/sunyd/identify/tomato_rnaseq/SRRall/cuff226-CPC_left.lncname.ej.del",header=F,stringsAsFactors = FALSE)
lncname <- as.matrix(lncname[,2])
lncnamedel <- as.matrix(lncnamedel[,2])
num1 <- 5
num2 <- 4
legendname <- c("Genes","LncRNAs","Intron","Inter","Anti")
legendname3 <- c("Genes","HC-lncRNAs","Other_lncRNAs")
color8 <- c("black","red","blue","green","purple","yellow")
color3 <- c("blue","red","green")

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
lncexondel <- exon[which(rownames(exon) %in% lncnamedel),,drop=FALSE] ###combinedexon[which(as.character(combinedexon$transcript_id) %in% genename$transcript_id),]
lncexon2 <- lncexon[-which(lncexon == 1),]
lncexon2all <- gtf[which(gtf[,2] %in% as.matrix(names(lncexon2))),,drop=FALSE]
for (i in 1:nrow(lncexon2all)){
    if (identical(lncexon2all[i,6],lncexon2all[i+1,5])){
              print(lncexon2all[i,2])
# }else{
 #             print(" ")
}
}
write.table(lncexon2all,"lnc_intron.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
##########calculate the cumulative number class of each trans exon num for barplot##########
barinput <- function(x,num){              ######the num differs depends on the data you operate
	tmp <- data.matrix(table(x[,1]))      ###how many trans contain 2 exon and so on
	tmp <- (c(tmp[1:num,],sum(tmp[(num+1):nrow(tmp),])) / sum(tmp)) * 100     ####display 1 to 15 and over 15
	return(tmp)
	}
lncbar <- barinput(lncexon,num1)
lncbar <- c(lncbar,rep(0,6-num1))
lncbardel <- barinput(lncexondel,num2)
lncbardel <- c(lncbardel,rep(0,6-num2))
genebar <- barinput(geneexon,6)
bar <- rbind(genebar,lncbar,lncbardel)
colnames(bar) <- c(1:6,c("> 6"))
pdf("lnc_exonnum.pdf")
barplot(bar,beside=T,main="Exons per transcript",ylab="Proportion",col=c("blue","red","green"),legend=c("Genes","HC-lncRNAs","Other_lncRNAs"))
dev.off()

#############the exon length distribution##########
exonlen <- as.matrix(gtf[,4])      ###########combinedgenename[which(as.character(combinedgenename$transcript_id) %in% genename$transcript_id),]
rownames(exonlen) <- gtf[,2]
geneexon.len <- exonlen[which(rownames(exonlen) %in% genename),,drop=FALSE]      ###
lncexon.len <- exonlen[which(rownames(exonlen) %in% lncname),,drop=FALSE]      ###
lncexon.dellen <- exonlen[which(rownames(exonlen) %in% lncnamedel),,drop=FALSE]      ###
#iexonlen <- as.matrix(lncname[which(lncname[,2] %in% iclass),4])         ###combinedgenename[which(as.character(combinedgenename$transcript_id) %in% lncname$transcript_id),]

pdf("lnc_exonlen.pdf")
plot(density(geneexon.len),type="l", xlab = "Exon length",ylab = "Density",col="blue",main="Exon length distribution")
lines(density(lncexon.len),col = "red")
lines(density(lncexon.dellen),col = "green")
#lines(density(uexonlen),col = "green")
#lines(density(oexonlen),col = "purple")
#lines(density(xexonlen),col = "pink")
#lines(density(jexonlen),col = "yellow")
#lines(density(eexonlen),col = "grey")
legend("topright",legend=legendname3,col=color3,bg="white",lwd=2)
dev.off()

plot2densityline <- function(input,inputtwo,xname,mainname,filename,legendname){
pdf(filename)
plot(density(input),type="l", xlab = xname,ylab = "Density",col="blue",main=mainname)
lines(density(inputtwo),col = "red")
legend("topright",legend=legendname,col=c("red","blue"),bg="white",lwd=2)
dev.off()
}
#plot2densityline(geneexonlen,lncexonlen,"exon length","Exon length distribution","lnc_exonlen.pdf",c("Long noncoding RNAs","Genes"))

#############the sequence composition information###############
#genename <- as.matrix(unique(gtf[,2]))                  ###############codeext(genename)
#lncname <- as.matrix(unique(lncname[,2]))                    ############codeext(lncname)

##########sequence length distribution plot##########
genelen <- seqcomp[genename,1,drop=FALSE] 
lnclen <- seqcomp[lncname,1,drop=FALSE]
lnclendel <- seqcomp[lncnamedel,1,drop=FALSE]
valu <- c(data.matrix(genelen),data.matrix(lnclen),data.matrix(lnclendel))
tim <- c(nrow(genelen),nrow(lnclen),nrow(lnclendel))
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
pdf("lnc_length.pdf")
boxplot(values ~ vars, data = df, col = color3, main="Transcript Length distribution")
dev.off()
pdf("lnc_length_only.pdf")
hist(data.matrix(lnclen),main="Histogram of lncRNAs length",xlab="Length",breaks=20,freq=FALSE)
dev.off()

print(summary(genelen))
print(summary(lnclen))
#valu <- c(data.matrix(genelen),data.matrix(lnclen),data.matrix(lnclen[iclass,,drop=FALSE]),data.matrix(lnclen[uclass,,drop=FALSE]),data.matrix(lnclen[oclass,,drop=FALSE]),data.matrix(lnclen[xclass,,drop=FALSE]),data.matrix(lnclen[jclass,,drop=FALSE]),data.matrix(lnclen[eclass,,drop=FALSE]))
#tim <- c(nrow(genelen),nrow(lnclen),nrow(lnclen[iclass,,drop=FALSE]),nrow(lnclen[uclass,,drop=FALSE]),nrow(lnclen[oclass,,drop=FALSE]),nrow(lnclen[xclass,,drop=FALSE]),nrow(lnclen[jclass,,drop=FALSE]),nrow(lnclen[eclass,,drop=FALSE]))

#plot2densityline(genelen,lnclen,"Length","Transcript length distribution","lnc_length.pdf",c("Genes","Long noncoding RNAs"))
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
lncgccumdel <- cumplot(seqcomp[lncnamedel,2])
#igccum <- cumplot(seqcomp[iclass,3])
#ugccum <- cumplot(seqcomp[uclass,3])
#ogccum <- cumplot(seqcomp[oclass,3])
#xgccum <- cumplot(seqcomp[xclass,3])
#jgccum <- cumplot(seqcomp[jclass,3])
#egccum <- cumplot(seqcomp[eclass,3])

pdf("lnc_au.pdf")
plot(genegccum[,1],genegccum[,2],type="l",col="blue",xlab="GC content",ylab="Cumulative frequency",main="The GC content of transcripts")
lines(lncgccum[,1],lncgccum[,2],type="l",col="red")
lines(lncgccumdel[,1],lncgccumdel[,2],type="l",col="green")
legend("bottomright",legend=legendname3,col=color3,bg="white",lwd=2)
dev.off()
#lines(igccum[,1],igccum[,2],type="l",col="blue")
#lines(ugccum[,1],ugccum[,2],type="l",col="green")
#lines(ogccum[,1],ogccum[,2],type="l",col="purple")
#lines(xgccum[,1],xgccum[,2],type="l",col="pink")
#lines(jgccum[,1],jgccum[,2],type="l",col="yellow")
#lines(egccum[,1],egccum[,2],type="l",col="grey")



contentplot <- function(col,xname,mainname,filename){
	lncgc <- seqcomp[lncname,col]
	genegc <- seqcomp[genename,col]
	lncgccum <- cumplot(lncgc)
	genegccum <- cumplot(genegc)
	pdf(filename)
	plot(genegccum[,1],genegccum[,2],type="l",col="blue",xlab=xname,ylab="Cumulative frequency",main=mainname)
	lines(lncgccum[,1],lncgccum[,2],type="l",col="red")
	dev.off()
	}
#contentplot(2,"GC content","The GC content of transcripts","lnc_gc.pdf")
#contentplot(3,"AU content","The AU content of transcripts","lnc_au.pdf")
#plot2densityline(data.matrix(seqcomp[genename,4]),data.matrix(seqcomp[lncname,4]),"GC/AU","GC/AU distribution","lnc_gc_au.pdf",c("Long noncoding RNAs","Genes"))

########calculate the class code count and extract different class categories###############
#classcode <- as.matrix(genename[,3])
#classcodecount <- as.matrix(table(classcode))

######show the transcript name and exon location info#############
#gtforder <- combinedgenename[order(combinedgenename[,1],combinedgenename[,6]),]
#exoncombine <- function(a){
#tmp <- data.frame(gene_id=NA,transcript_id=NA,class_code=NA,oId=NA,exonlength=NA,start=NA,end=NA,strand=NA)
#for (i in 1:nrow(a)){                #(nrow(gtforder)-1)){
#	x <- a[i,] 
#	y <- a[(i+1),]
#  if (x[1] == y[1] && y[6] >= x[6] && y[7] <= x[7]){
#	  tmp[i,] <- a[i,]
#	  a <- a[-(i+1),]
#  }else if (x[1] == y[1] && y[6] >= x[6] && y[7] >= x[7]){
#    tmp[i,] <- a[i,]
#	tmp[i,7] <- a[(i+1),7]
 #   a <- a[-(i+1),]
#  }else if (x[1] == y[1] && y[6] <= x[6] && y[7] <= x[7]){
#	tmp[i,]  <- a[i,]
#	tmp[i,6] <- a[(i+1),6]
#    a <- a[-(i+1),]
 # }else{
  #  tmp[i,] <- a[i,]
##	  }
#	}
#	return(tmp)
#}
#gtfcombine <- exoncombine(gtforder)


#gtforder[(i+1),1] == gtforder[i,1] && gtforder[(i+1),5] <= gtforder[i,5] && gtforder[(i+1),6] <= gtforder[i,6]
#lncrnafasta <- read.fasta(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.fa")
##loci <- read.delim("~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75_gene.loci",header=F)
##loci <- loci[,c(1,2,4,6,7)]
#library("seqinr")
#combinedfasta <- read.fasta(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75.combinedonelinename.fa") 
#########extract only geneid, transid, and classcode#######
