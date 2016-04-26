args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
filename <- args[2]

#########read in the gene names and lncrna names for analysis
lnc.high <- read.delim(file=args[1],header=T,row.names=1,stringsAsFactors = FALSE) 
lncnamepure.high<- as.matrix(rownames(lnc.high))
##########built a DGElist and store your data in it#################
library(TCC)

############calculate tissue specific score###########3
ts <- ROKU(lnc.high, upper.limit = 0.25, sort = FALSE)
out = ts$outlier
mod = as.matrix(ts$modH)
name = as.matrix(rownames(mod)[which(mod[,1] <= 2.5)])
filter = out[name,,drop=F]
coln = colnames(filter)
tissue.high = apply(filter,1,function(x){y <- coln[which(x == 1)];return(y)})
high.file = sapply(names(tissue.high),function(x) paste(x,paste(tissue.high[[x]],collapse=" ")))
write(high.file,filename)
tissue.low = apply(filter,1,function(x){y <- coln[which(x == -1)];return(y)})
low.file = sapply(names(tissue.low),function(x) paste(x,paste(tissue.low[[x]],collapse=" ")))
write(low.file,filename,append=T)
