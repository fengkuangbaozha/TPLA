
###########read the files from the command line in certain order ################
args <- commandArgs(trailingOnly = TRUE)
cpc <- read.delim(args[1], head = TRUE, quote="\"", dec=".", row.names = 1)
blast <- read.delim(args[2], head = TRUE, quote="\"", dec=".", row.names = 1)
#ff <- read.delim(args[3], head = TRUE, quote="\"", dec=".", row.names = 1)
comp <- read.delim(args[3], head = TRUE, quote="\"", dec=".", row.names = 1)
fold <- read.delim(args[4], head = TRUE, quote="\"", dec=".", row.names = 1)
lfold <- read.delim(args[5], head = TRUE, quote="\"", dec=".", row.names = 1)
filename <- args[6]
print(nrow(cpc))
print(nrow(blast))
#print(nrow(ff))
print(nrow(comp))
print(ncol(fold))
print(ncol(lfold))
#print(head(fold[,6]))
#print(head(lfold[,23]))

cpccom <- cbind(blast[,1:4],cpc[,2])
colnames(cpccom) <- c("blast_hit_number","blast_hit_segment_number","mean_hit_score_of_three_frames","frame_variation","coding_potential_value")
#print(head(cpccom))
#compare <- cbind(ff[,1],ff[,2],comp[,1],comp$orf_length,comp$orf_coverage)
#print(head(compare))
len <- as.matrix(comp[,1])
#gc <- as.matrix(comp[,2])
print(nrow(len))
#############divide the parameters of fold and lfold  by length################
loop <- function(x,n,length){
    x <- as.matrix(x)
     for (i in n){
        x[,i] <- x[,i] / length
	}
    return(x)
 }

foldn <- c(1,2,3,5,6)
#lfoldn <- c(9,10,11,12,13,14,23)
foldlen <- loop(fold,foldn,len)
#lfoldlen <- loop(lfold,lfoldn,len)
colnames(foldlen)[foldn] <- paste(colnames(foldlen)[foldn],"/Length",sep="")
#colnames(lfoldlen)[lfoldn] <- paste(colnames(lfoldlen)[lfoldn],"/Length",sep="")
#print(head(foldlen))
#print(head(lfoldlen))

feature <- cbind(cpccom,comp,foldlen,lfold)

print(head(feature))
name <- matrix(c("RNAname",colnames(feature)),nrow=1)

write.table(name, args[6], append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(feature, args[6], append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
