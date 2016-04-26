args <- commandArgs(trailingOnly = TRUE)                                                                                                                                     
library("seqinr")
training <- args[1]
input <- read.fasta(file = args[2])
foldfile <- args[3]
lfoldfile <- args[4]
output <- args[5]

fold <- read.fasta(foldfile)    ##read one gene information as a whole
lfold <- read.fasta(lfoldfile)

#####calculating the sequence compostion information ################

header <- matrix(c("RNA_name","sequence_length","orf_number","orf_min","orf_1st_Quarter","orf_Median","orf_Mean","orf_3rd_Quarter","orf_standard_error","orf_length","orf_coverage","orf_start_position","aaa_eff","aac_eff","aag_eff","aat_eff","aca_eff","acc_eff","acg_eff","act_eff","aga_eff","agc_eff","agg_eff","agt_eff","ata_eff","atc_eff","atg_eff","att_eff","caa_eff","cac_eff","cag_eff","cat_eff","cca_eff","ccc_eff","ccg_eff","cct_eff","cga_eff","cgc_eff","cgg_eff","cgt_eff","cta_eff","ctc_eff","ctg_eff","ctt_eff","gaa_eff","gac_eff","gag_eff","gat_eff","gca_eff","gcc_eff","gcg_eff","gct_eff","gga_eff","ggc_eff","ggg_eff","ggt_eff","gta_eff","gtc_eff","gtg_eff","gtt_eff","taa_eff","tac_eff","tag_eff","tat_eff","tca_eff","tcc_eff","tcg_eff","tct_eff","tga_eff","tgc_eff","tgg_eff","tgt_eff","tta_eff","ttc_eff","ttg_eff","ttt_eff"),nrow=1)
write.table(header, file = output, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   ####write the header in the first line

sequencecompostion <- function(x,filename){
for (i in 1:length(x)){
a <- x[[i]]
aname <- attr(a,"name")   ###the name of this gene
len <- length(a)  ### calculating the sequence length
#print(len)



#calculating ORF 
#find the position of those codons
findPotentialStartsAndStops <- function(sequence)
  {
     # Define a vector with the sequences of potential start and stop codons
     codons <- c("atg", "taa", "tag", "tga")
     # Find the number of occurrences of each type of potential start or stop codon
     for (i in 1:4)
     {
        codon <- codons[i]
        # Find all occurrences of codon "codon" in sequence "sequence"
        occurrences <- vmatchPattern(codon, sequence)
        # Find the start positions of all occurrences of "codon" in sequence "sequence"
        codonpositions <- unlist(startIndex(occurrences))
        # Find the total number of potential start and stop codons in sequence "sequence"
        numoccurrences <- length(codonpositions)
        if (i == 1)
        {
           # Make a copy of vector "codonpositions" called "positions"
           positions <- codonpositions
           # Make a vector "types" containing "numoccurrences" copies of "codon"
           types <- rep(codon, numoccurrences)
        }
        else
        {
           # Add the vector "codonpositions" to the end of vector "positions":
           positions   <- append(positions, codonpositions, after=length(positions))
           # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
           types       <- append(types, rep(codon, numoccurrences), after=length(types))
        }
     }
     # Sort the vectors "positions" and "types" in order of position along the input sequence:
     if (length(positions) == 0)
		{
			mylist <- list(0,0)
			return(mylist)
		}else{
	 indices <- order(positions)
     positions <- positions[indices]
     types <- types[indices]
     # Return a list variable including vectors "positions" and "types":
     mylist <- list(positions,types)
     return(mylist)
		}
  }

#find the start, end and length of the ORFs
findORFsinSeq <- function(sequence)
  {
     require(Biostrings)
     # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
     mylist <- findPotentialStartsAndStops(sequence)
     if (mylist[[1]] == 0){
	 mylist <- list(0,0,0)
     return(mylist)
	 }else{
     positions <- mylist[[1]]
     types <- mylist[[2]]
     # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
     orfstarts <- numeric()
     orfstops <- numeric()
     # Make a vector "orflengths" to store the lengths of the ORFs
     orflengths <- numeric()
     # Print out the positions of ORFs in the sequence:
     # Find the length of vector "positions"
     numpositions <- length(positions)
     # There must be at least one start codon and one stop codon to have an ORF.
     if (numpositions >= 2)
     {
        for (i in 1:(numpositions-1))
        {
           posi <- positions[i]
           typei <- types[i]
           found <- 0
           while (found == 0)
           {
              for (j in (i+1):numpositions)
              {
                 posj  <- positions[j]
                 typej <- types[j]
                 posdiff <- posj - posi
                 posdiffmod3 <- posdiff %% 3
                 # Add in the length of the stop codon
                 orflength <- posj - posi + 3
                 if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
                 {
                    # Check if we have already used the stop codon at posj+2 in an ORF
                    numorfs <- length(orfstops)
                    usedstop <- -1
                    if (numorfs > 0)
                    {
                      for (k in 1:numorfs)
                      {
                          orfstopk <- orfstops[k]
                          if (orfstopk == (posj + 2)) { usedstop <- 1 }
                      }
                    }
                    if (usedstop == -1)
                    {
                       orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                       orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                       orflengths <- append(orflengths, orflength, after=length(orflengths))
                    }
                    found <- 1
                    break
                 }
                 if (j == numpositions) { found <- 1 }
              }
           }
        }
     }
     # Sort the final ORFs by start position:
     indices <- order(orfstarts)
     orfstarts <- orfstarts[indices]
     orfstops <- orfstops[indices]
     # Find the lengths of the ORFs that we have
     orflengths <- numeric()
     numorfs <- length(orfstarts)
     for (i in 1:numorfs)
     {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        orflength <- orfstop - orfstart + 1
        orflengths <- append(orflengths,orflength,after=length(orflengths))
     }
     mylist <- list(orfstarts, orfstops, orflengths)
     return(mylist)
		}
  }

astr <- c2s(a)
#print(findPotentialStartsAndStops(astr))
orf <- findORFsinSeq(astr)
#print(orf)
orf3 <- orf[[3]]     #######put all the different ORF length in a vector
orf1 <- orf[[1]]   #######put all the ORF start positions in a vector

#####calculating ORF number, length, start position
orfnum <- length(orf3)    #####The number of possible ORFs
#print(orfnum)
orfsum <- as.matrix(summary(orf3))  ####the summary 
#print(orfsum)
std <- function(x) sd(x)/sqrt(length(x))  ####the standard error
orfstd <- std(orf3)
#print(orfstd)
orflen <- max(orf3)           #####the longest ORF length
orflencov <- orflen / len
#print(orflencov)
orfstart <- orf1[which(orf3 == max(orf3))]   ####the start position of longest ORF
#print(orfstart)   
orfstartper <- orfstart / len  ####the start position percentage of longest ORF
#print(orfstartper)
orfsumfour <- t(orfsum[1:5])    ###collect the first four elements in orfsum.
if(length(orfstart) == 0 || orfstart == 0){
	matrix2 = matrix(rep(0, 10),nrow=1)
	}else{
	matrix2 <- cbind(c(orfnum),orfsumfour,c(orfstd),c(orflen),c(orflencov),c(orfstartper))   #####put the orf indices into one matrix
}
colnames(matrix2) <-  c("orf_number","orf_min","orf_1st_Quarter","orf_Median","orf_Mean","orf_3rd_Quarter","orf_standard_error","orf_length","orf_coverage","orf_start_position")
#print(matrix2)


###calculating the base frequency
basecount <- function(x,num,s){
	countnum <- as.matrix(count(x, num, start = s, freq = TRUE))    #start means the start position and begins with 0.
	#countper <- countnum / sum(countnum)
	result <- t(countnum)
	return(result)
	}


########## calculating codon usage bias eff, freq, rsuc value #########
codonbias <- function(seq,s){
	countindices <- uco(seq, frame = s, index = c("eff"), as.data.frame = TRUE,NA.rscu = 0) ##frame means the start position and begins with 0.
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
	index1 <- codonbias(a,0)
	index2 <- codonbias(a,1)
	index3 <- codonbias(a,2)
	index4 <- codonbias(rev(a),0)
	index5 <- codonbias(rev(a),1)
	index6 <- codonbias(rev(a),2)
	mat <- (index1 + index2 + index3 + index4 + index5 + index6) 
	return(mat)

}

matrix5 <- tri_sum()
#print(matrix5)

######combine several several matrix into a sum and write to a file ##################
matrixsum <- cbind(as.matrix(len),matrix2,matrix5)
rownames(matrixsum) <- aname    ###put the gene name as the rowname
print(aname)

write.table(matrixsum, file = filename, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)   ### write the value in the following lines
}
}
sequencecompostion(input,output)

##########format the RNA structure value#############
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
		p <- paste(rnafold[7],rnafold[8], sep = "e")     ##combine the 7 and 8 line because they Evalue: mfe_frequency
		rnafold <- matrix(c(rnafold[1:6],p,rnafold[9]), nrow=1)	##create all the indexes in the rnafold method
		}else{
			rnafold
			}
	rownames(rnafold) <- attr(a,which="name")               ##assign the gene name for the matrix
	colnames(rnafold) <- c("minimum_free_energy","ensemble_energy","centroid_strucformatture","distance_to_ensemble_energy","MEA_energy","MEA_value","mfe_frequency","ensemble_diversity")                       #####assign the colnames for the rnafold index
	name <- matrix(c("RNAname", colnames(rnafold)), nrow=1)
write.table(rnafold[,c(1,2,3,4,5,6,8),drop=F], filename, sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
}
}

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
	write.table(result, filename, sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
	}else{
		result <- matrix(c(0,rep(0,21),rnalfold), nrow=1)
		sumname <- c("min","1st_quarter","median","mean","3rd_quarter","max")
		mfe_name <- paste("locally_mfe", sumname, sep = "_")  ##assign the colnames for the matrix.
		start_name <- paste("locally_start_position", sumname, sep = "_")
		z_score_name <- paste("locally_z_score", sumname, sep = "_")
		rownames(result) <- attr(a,which="name")
		colnames(result) <- c("locally_stable_number", mfe_name, "mfe_std", start_name, "start_position_std", z_score_name, "z_score_std", "locally_stable_energy")
		result
	write.table(result, filename, sep='\t', append = TRUE, quote = FALSE, row.names = TRUE, col.names = FALSE)
}
}
}

name1 <- matrix(c("RNAname","minimum_free_energy","ensemble_energy","centroid_strucformatture","distance_to_ensemble_energy","MEA_energy","MEA_value","ensemble_diversity"),nrow=1)                       #####assign the colnames for the rnafold index
write.table(name1, paste(foldfile,".format",sep=""), sep='\t', append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
name2 <- matrix(c("RNAname","locally_stable_number","locally_mfe_min","locally_mfe_1st_quarter","locally_mfe_median","locally_mfe_mean","locally_mfe_3rd_quarter","locally_mfe_max","mfe_std","locally_start_position_min","locally_start_position_1st_quarter","locally_start_position_median","locally_start_position_mean","locally_start_position_3rd_quarter","locally_start_position_max","start_position_std","locally_z_score_min","locally_z_score_1st_quarter","locally_z_score_median","locally_z_score_mean","locally_z_score_3rd_quarter","locally_z_score_max","z_score_std","locally_stable_energy"),nrow=1)
write.table(name2, paste(lfoldfile,".format",sep=""), sep='\t', append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)

rnafold(fold,paste(foldfile,".format",sep=""))    ##read one gene information as a whole
rnalfold(lfold,paste(lfoldfile,".format",sep=""))

######combine the structure and sequence info together##########
###########read the files from the command line in certain order ################
comp <- read.delim(file=output, head = TRUE, quote="\"", dec=".", row.names = 1, as.is = FALSE)      ##validation dataset
fold <- read.delim(paste(foldfile,".format",sep=""), head = TRUE, quote="\"", dec=".", row.names = 1)
lfold <- read.delim(paste(lfoldfile,".format",sep=""), head = TRUE, quote="\"", dec=".", row.names = 1)
len <- as.matrix(comp[,1])
#############divide the parameters of fold and lfold  by length################
loop <- function(x,n,length){
    x <- as.matrix(x)
     for (i in n){
        x[,i] <- x[,i] / length
	}
    return(x)
 }

foldn <- c(1,2,3,5,6)
lfoldn <- c(9,10,11,12,13,14,23)
foldlen <- loop(fold,foldn,len)
lfoldlen <- loop(lfold,lfoldn,len)
colnames(foldlen)[foldn] <- paste(colnames(foldlen)[foldn],".Length",sep="")
colnames(lfoldlen)[lfoldn] <- paste(colnames(lfoldlen)[lfoldn],".Length",sep="")

feature <- cbind(comp[order(rownames(comp)),,drop=F],foldlen[order(rownames(foldlen)),,drop=F],lfoldlen[order(rownames(lfoldlen)),,drop=F])
print(head(feature))

#######select top 50 performance features################
com2 <- feature
com2[is.na(com2)] <- 0

##########start to train and predict########
library("caret")
###############performance measure of the predictions###############
performance_measure <- function(model,testset,filename){
	prediction1 <- predict(model, newdata=testset,type = "prob")
	prediction2 <- predict(model, newdata=testset)
    result <- cbind(prediction1[,1,drop=FALSE],as.matrix(prediction2))
    #colnames(result) <- parameter
    colnames(result) <- c("lncRNA_potential","Type")
    result <- result[(result[,2] == "lncRNA"),]
	write.table(result, file = filename, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
	return(result)
}

###########evaluating com1 performance on self test data###############"E:/machine_learning/caret_all"
library("randomForest")
load(training)
performance_measure(rf_model,com2,paste(output,"_score.txt",sep=""))     ##  --means implicit feature selection
#system(paste("rm", output))
   

