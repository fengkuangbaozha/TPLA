args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     


all_length <-read.delim(file=args[1],header=F)    ###read in the cuffnorm gene count info
filename <- args[2]

pdf(paste(filename,".pdf",sep=""))
hist(all_length[,1])
dev.off()
