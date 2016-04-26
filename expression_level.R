args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
#allcount.avs.high <-read.delim(file=args[1],row.names=1,header=F)
#allgroup <- read.delim(args[2],header=F)     ####c("Ve",rep("EC",3),rep("Sp",3),rep("Ve",2))
#allgroup <- allgroup[,2]
filename <- args[3]

#########read in the gene names and lncrna names for analysis
gene.high <- read.delim(file=args[1],header=T,row.names=1) 
genenamepure.high<- as.matrix(rownames(gene.high))
lnc.high <- read.delim(file=args[2],header=T,stringsAsFactors = FALSE,row.names=1) 
lncnamepure.high<- as.matrix(rownames(lnc.high))
log2allcount.avs.high <- rbind(gene.high,lnc.high)
#legendname <- c("Genes","LncRNAs","Intron","Inter","Exon","Anti","iso","match")
#color8 <- c("black","red","blue","green","purple","pink","yellow","grey")
COLA="#ab2023"
COLB="#1c76b1"
legendname3 <- c("lncRNA","PCgene")
color3 <- c(COLA,COLB)

############calculate maximum cpm value for gene expression plot analysis###########
expgenehigh <- apply(gene.high,1,max)
explnchigh <- apply(lnc.high,1,max)
valu <- c(data.matrix(explnchigh),data.matrix(expgenehigh))
tim <- c(length(explnchigh),length(expgenehigh))   
df <- data.frame(values = valu,vars = rep(legendname3, times = tim))
pdf(paste(filename,"_expression.pdf",sep=""))
boxplot(values ~ vars, data = df, col = color3,main="Transcript expression distribution (Maximum)",ylab="log2FPKM Maximum",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
dev.off()

