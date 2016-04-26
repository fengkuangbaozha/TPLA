###########read the files from the command line in certain order ################
args <- commandArgs(trailingOnly = TRUE)
tracking <- read.delim(args[1], head = FALSE, quote="\"", dec=".", row.names = 1)###########reading the tracking loci file
filename <- args[2]
#tracking1 <- tracking[,4:ncol(tracking)]

#print(c("The tracking files contains ",nrow(tracking)))
############delete those FPKM < 2 in all exps in tracking file####################
for (i in 4:ncol(tracking)){
	tracking[,i] <- apply(tracking[,i,drop=FALSE],1,function(x){return(unlist(strsplit(as.character(x),"[|]"))[4])})  ###extract the FPKM value in tracking
}
tracking[is.na(tracking)] <- 0                  ######convert the lost FPKM value into 0
tracking2 <- apply(tracking[,4:ncol(tracking)],1,function(x){all(x < 2)})      ####all FPKM < 2 are labeled as TRUE and are deleted
tracking3 <- tracking[!tracking2,]
print(c("The FPKM remain tracking files contains ",nrow(tracking3)))
keepnames <- rownames(tracking3)
write.table(as.matrix(keepnames), filename, append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
############extract the keeping FPKM names in gtf file ############
#combined.gtf1 <- apply(combined.gtf[,9,drop=FALSE],1,function(x){return(unlist(strsplit(as.character(x),"[ ]"))[4])}) ###########extract the TCONS number in combined.gtf    
#combined.gtf1 <- gsub(";","",combined.gtf1)                              #############replace the ;
#location <- match(combined.gtf1,keepnames)
#combined.gtf2 <- combined.gtf[!is.na(location),]
#print(c("The final remain combined.gtf files contains ",nrow(combined.gtf2)))
#write.table(combined.gtf2, filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
##"/psc/bioinformatics/sunyd/lncrna/combine.feature"
#keepnum <- tracking3[,1]
#combined.gtftmp <- apply(combined.gtf[,9,drop=FALSE],1,function(x){return(unlist(strsplit(as.character(x),"[ ]"))[2])})
#combined.gtftmp <- gsub(";","",combined.gtftmp)
#locationtmp <- match(combined.gtftmp,keepnum)
#combined.gtf2tmp <- combined.gtf[!is.na(locationtmp),]
#write.table(combined.gtf2tmp, "cuffRC.combinedhq_XLOC.gtf", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#index <- which(duplicated(combined.gtf1))   #######check the combined.gtf duplicated line for TCONS
#tmp <- combined.gtf1[-index]                #######delete the duplicated rows

