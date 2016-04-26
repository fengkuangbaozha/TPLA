args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
library("seqinr")
fold <- read.fasta(args[1])    ##read one gene information as a whole
lfold <- read.fasta(args[2])

name1 <- matrix(c("RNAname","minimum_free_energy","ensemble_energy","centroid_strucformatture","distance_to_ensemble_energy","MEA_energy","MEA_value","mfe_frequency","ensemble_diversity"),nrow=1)                       #####assign the colnames for the rnafold index
write.table(name1, args[3], sep='\t', append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
name2 <- matrix(c("RNAname", "locally_stable_number"),nrow=1)
write.table(name2, args[4], sep='\t', append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

################extract numbers only #######################
extractnum <- function(x){
	astring <- c2s(x)                                        ## convert the list into a string
	numposition <- gregexpr("[-+]?([0-9]*[.][0-9]+|[0-9]+)", astring, perl = TRUE)    ###find the position of the number and the match length
	number <- regmatches(astring, numposition)             ##extract the number 
	number <- matrix(unlist(number),nrow = 1)              ##conver to a matrix
	return(number)
}





#########  Extracting the number from the RNAfold result #################
rnafold <- function(x,filename){
for(i in 1:length(x)){
	a <- x[[i]]                                      ## extract one gene information to calculate
	rnafold <- extractnum(a)                         ## extract numbers only in a gene
	rnafold8 <- gregexpr("^-[0-9]*", rnafold[8], perl = TRUE)
#	rnafold8test <- attr(rnafold8[[1]], which="useBytes")
	if(rnafold8[[1]] == 1){
		p <- paste(rnafold[7],rnafold[8], sep = "e")     ##combine the 7 and 8 line because they Evalue
		rnafold <- matrix(c(rnafold[1:6],p,rnafold[9]), nrow=1)	##create all the indexes in the rnafold method
		}else{
			rnafold
			}
	rownames(rnafold) <- attr(a,which="name")               ##assign the gene name for the matrix
	colnames(rnafold) <- c("minimum_free_energy","ensemble_energy","centroid_strucformatture","distance_to_ensemble_energy","MEA_energy","MEA_value","mfe_frequency","ensemble_diversity")                       #####assign the colnames for the rnafold index
	name <- matrix(c("RNAname", colnames(rnafold)), nrow=1)
#write.table(name, filename, sep='\t', append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(rnafold, filename, sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
}
}

#rnafold(foldmrna2,"~/rna/lncrna/feature/rnafoldmrna_2format.txt")
std <- function(x) sd(x)/sqrt(length(x))  ####the standard error
############ extract the locally stable file number ########################                       
rnalfold <- function(x,filename){
for(i in 1:length(x)){
	a <- x[[i]]                                      ## extract one gene information to calculate
	rnalfold <- extractnum(a)                        ## extract numbers only in a gene
	len <- length(rnalfold)
	if(len > 1){
	countnum <- as.matrix((length(rnalfold)-1) / 3)             ## calculate the total number of locally stable strucformatture
	po <- seq(1, (length(rnalfold)-1), by = 3)		 ## 
	mfe <- matrix(summary(as.numeric(rnalfold[po])),nrow=1)   ##make all locally mfe in one line, and calculate the summary information of it, display as a one row matrix
	start <- matrix(summary(as.numeric(rnalfold[po + 1])),nrow=1)
	z_score <- matrix(summary(as.numeric(rnalfold[po + 2])),nrow=1)
	mfe_std <- as.matrix(std(as.numeric(rnalfold[po])))
	start_std <- as.matrix(std(as.numeric(rnalfold[po + 1])))
	z_score_std <- as.matrix(std(as.numeric(rnalfold[po + 2])))
	sumname <- c("min","1st_quarter","median","mean","3rd_quarter","max")
	mfe_name <- paste("locally_mfe", sumname, sep = "_")  ##assign the colnames for the matrix.
	start_name <- paste("locally_start_position", sumname, sep = "_")
	z_score_name <- paste("locally_z_score", sumname, sep = "_")
	result <- cbind(countnum,mfe,mfe_std,start,start_std,z_score,z_score_std,as.matrix(rnalfold[len]))
	rownames(result) <- attr(a,which="name")
	colnames(result) <- c("locally_stable_number", mfe_name, "mfe_std", start_name, "start_position_std", z_score_name, "z_score_std", "locally_stable_energy")
	result
#	name <- matrix(c("RNAname", "locally_stable_number", mfe_name, "mfe_std", start_name, "start_position_std", z_score_name, "z_score_std", "locally_stable_energy"), nrow=1)
#	write.table(name, filename, sep='\t', append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
	write.table(result[,1,drop=FALSE], filename, sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
	}else{
		result <- matrix(c(0,rep(0,21),rnalfold), nrow=1)
		sumname <- c("min","1st_quarter","median","mean","3rd_quarter","max")
		mfe_name <- paste("locally_mfe", sumname, sep = "_")  ##assign the colnames for the matrix.
		start_name <- paste("locally_start_position", sumname, sep = "_")
		z_score_name <- paste("locally_z_score", sumname, sep = "_")
		rownames(result) <- attr(a,which="name")
		colnames(result) <- c("locally_stable_number", mfe_name, "mfe_std", start_name, "start_position_std", z_score_name, "z_score_std", "locally_stable_energy")
		result
	write.table(result[,1,drop=FALSE], filename, sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
}
}
}

#rnalfold(lfoldmrna2, "~/rna/lncrna/feature/rnalfoldmrna_2format.txt")


rnafold(fold,args[3])    ##read one gene information as a whole
#write.table(name1, args[3], sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
rnalfold(lfold,args[4])
#write.table(lfoldfor, args[4], sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)






   
