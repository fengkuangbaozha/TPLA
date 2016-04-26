args <- commandArgs(trailingOnly = TRUE)
gtf <- read.delim(file=args[1],header=F,stringsAsFactors = FALSE) 
summary(gtf)
