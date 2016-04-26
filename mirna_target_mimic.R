#mimic.candi <- read.delim(file = "~/sunyd/identify/oryza_rnaseq/SRRstressheat/rice.mirna/psrobot.mirna.oneline1",header=F,stringsAsFactors=FALSE)
#filename <- "~/sunyd/identify/oryza_rnaseq/SRRstressheat/rice.mirna/psrobot.mirna1.tmp"
args <- commandArgs(trailingOnly = TRUE)
filename=args[2]
mimic.candi <- read.delim(file = args[1],header=F,stringsAsFactors=FALSE)
library(parallel)
cl <- makeCluster(getOption("cl.cores", 32))
target.location = apply(mimic.candi[,5,drop=FALSE],1,function(x) strsplit(x,"")[[1]])      ######extract the miRNA target section and split into single character
#target.keep = mimic.candi[which(parLapply(cl,target.location,function(x) any(x[9:12]=="-") && sum(x[9:15] == "-")==3)==TRUE),]     #####only keep 9-12 position is -
target.keep = mimic.candi[which(parLapply(cl,target.location,function(x) sum(x[9:12] == "-")==3 || sum(x[10:13] == "-")==3 || sum(x[11:14] == "-")==3 || sum(x[12:15] == "-")==3)==TRUE),]     #####only keep 9-12 position is -
align.location = apply(target.keep[,10,drop=FALSE],1,function(x) strsplit(x,"")[[1]])   #####extract the aligment quality column and split into single character
align.keep = target.keep[which(parLapply(cl,align.location,function(x) sum(x != "|")) <= 6),]           ######only keep mismatch and GU pair less 6
align.keep.location = apply(align.keep[,10,drop=FALSE],1,function(x) strsplit(x,"")[[1]])   #####extract the aligment quality column and split into single character
align.keep.last = align.keep[which(parLapply(cl,align.keep.location,function(x) all(x[2:8] == "|") || x[1] == "|") == TRUE),]           ######only keep 2:8 is perfect match
write.table(align.keep.last,paste(filename,".mimic.name",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
