args <- commandArgs(trailingOnly = TRUE)
library("seqinr")
library(parallel)
combinedfasta <- read.fasta(file = args[1])        ###"~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75.combinedonelinename.fa") 
filename <- args[2]
#lncrnafasta <- read.fasta(file = "~/sunyd/identify/cucumber_rnaseq/SRRadd/cuff75-CPC_left.fa")

#####calculating the sequence compostion information ################
cl <- makeCluster(getOption("cl.cores", 32))
len <- parLapply(cl, combinedfasta, length)
len.frame <- as.matrix(sapply(len,c))

gc <- parLapply(cl, combinedfasta, GC)
gc.frame <- as.matrix(sapply(gc,c))
au <- 1 - gc.frame
gc2au <- gc.frame / au

codon.sum <- function(a){
###calculating the base frequency
basecount <- function(x,num,s){
	countnum <- as.matrix(count(x, num, start = s, freq = TRUE))    #start means the start position and begins with 0.
	#countper <- countnum / sum(countnum)
	result <- t(countnum)
	return(result)
	}

########## calculating codon usage bias eff, freq, rsuc value #########
codonbias <- function(seq,s){
	countindices <- uco(seq, frame = s, index = c("eff", "freq", "rscu"), as.data.frame = TRUE,NA.rscu = 0) ##frame means the start position and begins with 0.
	triname <- rownames(countindices)
	effname <- paste(triname, "eff", sep = "_")   ###change the name of index
	eff <- matrix(countindices[,3],nrow=1)
	colnames(eff) <- effname                     ###convert the names into the changed name
	return(eff)
	}

#########  calculating  the tri-bases index #########
tri_index <- function(x,s){
	bs1 <- basecount(x,3,s)   ###calculating the tri-bases frequency in different positions
	cs1 <- codonbias(x,s)     ###measuring the codon usage bias indices
	m <- cbind(bs1,cs1)	###combine them into one row matrix
	return(m)	
	}
	
#######  choose the ORF order for three nucleotides #########
tri_sum <- function(x){
	index1 <- tri_index(a,0)
	index2 <- tri_index(a,1)
	index3 <- tri_index(a,2)
	index4 <- tri_index(a,3)
	index5 <- tri_index(a,4)
	index6 <- tri_index(a,5)
	mat <- (index1 + index2 + index3 + index4 + index5 + index6)
	return(mat)
}
matrix2 <- tri_sum()
return(matrix2)
}
code <- parLapply(cl, combinedfasta, codon.sum)
code.frame <- t(as.matrix(sapply(code,c)))
colnames(code.frame) <- colnames(code[[1]])

matrixsum <- cbind(len.frame,gc.frame,au,gc2au,code.frame)

colnames(matrixsum) <- c("length","GC","AU","GC2AU",colnames(code.frame))
write.table(matrixsum, file = args[2], append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)   ####write the header in the first line


sequencecompostion <- function(x,filename){
 matall <- matrix(NA,ncol=132,nrow=length(x)) 
 for (i in 1:length(x)){
  a <- x[[i]]
  aname <- attr(a,"name")   ###the name of this gene
  len <- length(a)  ### calculating the sequence length
  print(len)
####calculating the  GC, AU and the divided content ###########
gccontent <- function(x){
	gc <- GC(x)       ### calculating the GC content
	au <- 1 - gc      ### calculating the AU content
	gc2au <- gc / au  ### calculating the GC/AU value
	result <- matrix(c(gc,au,gc2au),nrow=1)
	return(result)
}
gc <- gccontent(a)   # calculating the GC content of the sequence
print(gc)
matrix1 <- cbind(c(len),gc)  ##put the length and GC content indices into a matrix
colnames(matrix1) <- c("sequence_length","GC_cotent","AU_content","GC/AU")
print(matrix1)

###calculating the base frequency
basecount <- function(x,num,s){
	countnum <- as.matrix(count(x, num, start = s, freq = TRUE))    #start means the start position and begins with 0.
	#countper <- countnum / sum(countnum)
	result <- t(countnum)
	return(result)
	}

########## calculating codon usage bias eff, freq, rsuc value #########
codonbias <- function(seq,s){
	countindices <- uco(seq, frame = s, index = c("eff", "freq", "rscu"), as.data.frame = TRUE,NA.rscu = 0) ##frame means the start position and begins with 0.
	triname <- rownames(countindices)
	effname <- paste(triname, "eff", sep = "_")   ###change the name of index
	eff <- matrix(countindices[,3],nrow=1)
	colnames(eff) <- effname                     ###convert the names into the changed name
	return(eff)
	}

#########  calculating  the tri-bases index #########
tri_index <- function(x,s){
	bs1 <- basecount(x,3,s)   ###calculating the tri-bases frequency in different positions
	cs1 <- codonbias(x,s)     ###measuring the codon usage bias indices
	m <- cbind(bs1,cs1)	###combine them into one row matrix
	return(m)	
	}
	
#######  choose the ORF order for three nucleotides #########
tri_sum <- function(x){
	index1 <- tri_index(a,0)
	index2 <- tri_index(a,1)
	index3 <- tri_index(a,2)
	index4 <- tri_index(a,3)
	index5 <- tri_index(a,4)
	index6 <- tri_index(a,5)
	mat <- (index1 + index2 + index3 + index4 + index5 + index6) 
	return(mat)
}
matrix2 <- tri_sum()
#print(matrix2)

######combine several several matrix into a sum and write to a file ##################
matrixsum <- cbind(matrix1,matrix2)
rownames(matrixsum) <- aname    ###put the gene name as the rowname
#print(matrixsum)
matall[i,] <- matrixsum
######### set the header for the features in the file ################
colsname <- matrix(colnames(matrixsum),nrow=1)      #######assign the name of the matrix to another one 
header <- cbind(c("RNA_name"),colsname)            #######print the header of the file
#write.table(header, file = filename, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)   ####write the header in the first line

write.table(matrixsum, file = filename, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)   ### write the value in the following lines
}
return(header)
}
#wholemat <- sequencecompostion(combinedfasta,args[2])           ##"~/sunyd/identify/cucumber_rnaseq/SRRadd/lnclength.txt")

#sequencecompostion(lncrna,"~/sunyd/identify/cucumber_rnaseq/SRRadd/length.txt")
#if ((x - 1) %% 3 == 0 || length(x) == 0) {
#  index1 <- tri_index(a)	
#	index4 <- tri_index(atriloc4)
#	mat1 <- cbind(index1,index4)    ####combine the two conditions together
#	return(mat1)
#	}else if ((x + 1) %% 3 == 0){
#	index2 <- tri_index(atriloc2)
#	index5 <- tri_index(atriloc5)
#	mat2 <- cbind(index2,index5)
#	return(mat2)
