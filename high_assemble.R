###########read the files from the command line in certain order ################
args <- commandArgs(trailingOnly = TRUE)
loci <- read.delim(args[1], head = FALSE, quote="\"", dec=".", row.names = 1)###########reading the compare loci file
tracking <- read.delim(args[2], head = FALSE, quote="\"", dec=".", row.names = 1)###########reading the tracking loci file
combined.gtf <- read.delim(args[3], head = FALSE, quote="\"", dec=".")###########reading the tracking loci file
number <- args[5]        ###sample number -1
print(c("The loci files contains ",nrow(loci)))
print(c("The tracking files contains ",nrow(tracking)))
print(c("The combined.gtf files contains ",nrow(combined.gtf)))
#############extract the transcripts appeared at least two experiments in loci file and only keep those in tracking file###############
loci1 <- apply(loci,1,function(x){return(length(which(x == "-")) > number)})    ####extract the one appered transcripts, deleting those "-" appeared over sample-1 times
locikeep <- as.matrix(names(loci1[loci1 == FALSE]))                    ##########extract the name of over one experiment appeared transcripts
twoexp_remain <- match(tracking[,1],locikeep)                          ##########match the two_exps transcripts with the tracking file first column
tracking1 <- tracking[!is.na(twoexp_remain),]                            ###########keep the two_exps transcripts in tracking file
print(c("The two exps remain tracking files contains ",nrow(tracking1)))
#############delete the classcode for "e","p","r","s","." in tracking file since they are low quality for our study################
del <- c("e","p","r","s",".")
del_track <- match(tracking1[,3],del)
tracking2 <- tracking1[is.na(del_track),]
print(c("The class code remain tracking files contains ",nrow(tracking2)))
############delete those FPKM < 2 in all exps in tracking file####################
for (i in 4:ncol(tracking2)){
	tracking2[,i] <- apply(tracking2[,i,drop=FALSE],1,function(x){return(unlist(strsplit(as.character(x),"[|]"))[4])})  ###extract the FPKM value in tracking2
}
tracking2[is.na(tracking2)] <- 0                  ######convert the lost FPKM value into 0
#tmp1 <- apply(tmp[,4:6] < 2,1,any)
tracking3 <- apply(tracking2[,4:ncol(tracking2)],1,function(x){all(x < 2)})      ####all FPKM < 2 are labeled as TRUE and are deleted
tracking3 <- tracking2[!tracking3,]
print(c("The FPKM remain tracking files contains ",nrow(tracking3)))
keepnames <- rownames(tracking3)
############extract the keeping FPKM names in gtf file ############
combined.gtf1 <- apply(combined.gtf[,9,drop=FALSE],1,function(x){return(unlist(strsplit(as.character(x),"[ ]"))[4])}) ###########extract the TCONS number in combined.gtf    
combined.gtf1 <- gsub(";","",combined.gtf1)                              #############replace the ;
location <- match(combined.gtf1,keepnames)
combined.gtf2 <- combined.gtf[!is.na(location),]
print(c("The final remain combined.gtf files contains ",nrow(combined.gtf2)))
write.table(combined.gtf2, args[4], append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
##"/psc/bioinformatics/sunyd/lncrna/combine.feature"
#keepnum <- tracking3[,1]
#combined.gtftmp <- apply(combined.gtf[,9,drop=FALSE],1,function(x){return(unlist(strsplit(as.character(x),"[ ]"))[2])})
#combined.gtftmp <- gsub(";","",combined.gtftmp)
#locationtmp <- match(combined.gtftmp,keepnum)
#combined.gtf2tmp <- combined.gtf[!is.na(locationtmp),]
#write.table(combined.gtf2tmp, "cuffRC.combinedhq_XLOC.gtf", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
index <- which(duplicated(combined.gtf1))   #######check the combined.gtf duplicated line for TCONS
tmp <- combined.gtf1[-index]                #######delete the duplicated rows

