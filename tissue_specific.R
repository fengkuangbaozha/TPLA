args <- commandArgs(trailingOnly = TRUE)
exp <- read.delim(args[1], head = TRUE, dec=".", row.names = 1, as.is = FALSE)
spec <- read.delim(args[2], head = FALSE, dec=".", row.names = 1, as.is = FALSE)
filename = args[3]
name = as.matrix(rownames(spec)[which(spec[,1] <= 3)])  ######select tissue specific score less than 3
nrow(name)
filter = exp[name,,drop=F]
tissue = as.matrix(apply(as.matrix(apply(filter,1, which.max)),1,function(x){colnames(filter)[x]}))
nrow(tissue)
write.table(tissue, filename,sep="\t", row.names=TRUE, col.names=FALSE,quote=FALSE)
non = spec[-which(spec[,1] <= 3),,drop=F]
if(nrow(non)== 0){
	non = NULL
	}else{
	non[,1] = "non-specific"
	nrow(non)
	write.table(non, filename,sep="\t",append=T, row.names=TRUE, col.names=FALSE,quote=FALSE)
}
