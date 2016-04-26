args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
top <- read.delim("/psc/bioinformatics/sunyd/lncrna/com1_top50", head = FALSE, quote="\"", as.is = FALSE)
com1 <- read.delim(args[1], head = TRUE, quote="\"", dec=".", row.names = 1, as.is = FALSE)      ##validation dataset

com <- com1[,as.matrix(top)]
name <- matrix(c("RNAname",colnames(com)),nrow=1)
write.table(name, file = args[2], append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(com, file = args[2], append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

