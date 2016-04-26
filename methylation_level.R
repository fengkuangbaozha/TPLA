args <- commandArgs(trailingOnly = TRUE)
file <- read.delim(args[1],header=F,stringsAsFactors=F)
filename = args[2]
name = args[3]
#file <- read.delim("~/sunyd/identify/oryza_rnaseq/SRRpi/methylation/all.windows.intersect.cg",header=F,stringsAsFactors=F)
COLA="#ab2023"
COLB="#1c76b1"

promotor = file[grep("promotor",file[,1]),]
genebody = file[grep("genebody",file[,1]),]
terminal = file[grep("terminal",file[,1]),]
class_count <- function(x,class){
	x_count <-  aggregate(x['ratio'],by=x['class'],mean)
	x_count$type <-  rep(class,nrow(x_count))
	return(x_count)
}

type_count <- function(file,num){
c_count = aggregate(file['V3'],by=file['V1'],sum)
ct_count = aggregate(file['V4'],by=file['V1'],sum)
class = aggregate(file['V2'],by=file['V1'],mean)
methy_level = cbind(class[order(class[,1]),,drop=F],c(c_count[order(c_count[,1]),2]/ct_count[order(ct_count[,1]),2]))
colnames(methy_level) = c("genename","class","ratio")
lncrna = methy_level[grep("TCONS",methy_level[,1]),]
lncrna_class_count = class_count(lncrna,"lncRNA")
gene = rbind(methy_level[grep("LOC",methy_level[,1]),],methy_level[grep("Chr",methy_level[,1]),])
gene_class_count = class_count(gene,"gene")
all_count = rbind(gene_class_count,lncrna_class_count)
all_count$class = all_count$class + num
return(all_count)
#write.table(methy_level,args[2],quote=F,sep="\t",col.names=FALSE,row.names=FALSE)
}
promotor_count = type_count(promotor,0)
genebody_count = type_count(genebody,100)
terminal_count = type_count(terminal,200)
count_all = rbind(promotor_count,genebody_count,terminal_count)
png(paste(filename,".png",sep=""))
plot(x=NULL,y=NULL,xlim=c(1,300),xaxt="n",ylim=c(0.0,0.6),type="l",col="green",lty=1,lwd=2,xlab="",ylab="",main="",cex.axis=1.5)
axis(side=1, at=seq(1,300,by=100), labels = FALSE, tck=-0.01)
axis(side=1, at=seq(50,250,by=100), labels=c("upstream 1kb","genebody","downstream 1kb"), tck=-0.00,cex.axis=1.5)
lines(smooth.spline(count_all[which(count_all$type=="lncRNA"),1:2]),type="l",col=COLA,lty=1,lwd=2)
lines(smooth.spline(count_all[which(count_all$type=="gene"),1:2]),type="l",col=COLB,lty=1,lwd=2)
abline(v=100,col="green",lty=2)
abline(v=200,col="green",lty=2)
text(150, 0.6, labels=name,pos=1,cex=2)
#legend("topright",col=c(COLB,COLA),legend=c("PCgene","lncRNA"),lty=1,cex=1.2)
dev.off()
#le_name <- c("rf","ORFpls","Boruta","RRF","treebag","bagFDA","ada","nodeHarvest","svmRadial","knn","nb","J48","rpart","C5.0","JRip","nnet","earth","gam")
#legend("bottomright",legend=le_name,lty=c(rep(1,9),rep(3,9)),col=c("red","black","yellow","blue","green","orange","purple","darkgrey","pink"),bg="white",lwd=2)
#####plot gene and lncrna methylation level#######
#valu <- as.numeric(c(data.matrix(gene[,2]),data.matrix(lncrna[,2])))
#tim <- c(nrow(gene),nrow(lncrna))
#df <- data.frame(values = valu,vars = rep(c("Gene","LncRNA"), times = tim))
#png(paste(args[3],"png",sep="."))
#boxplot(values ~ vars, data = df, col = c("blue","red"), main="DNA methylation distribution")
#dev.off()

