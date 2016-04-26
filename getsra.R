library(SRAdb)
args <- commandArgs(trailingOnly = TRUE)
file <- read.delim(args[1],header=F)
#file <- read.delim("~/sunyd/identify/tomato_rnaseq/sra.file",header=F)
setwd('/psc/bioinformatics/sunyd/identify/download')
#########building the sra sql library#########
sqlfile <- 'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
timeStart <- proc.time()
proc.time() - timeStart
sra_con <- dbConnect(SQLite(),sqlfile)

##########get your interested SRP info ######################
#conversion <- sraConvert( file[,1], sra_con = sra_con )     ######get the SRP corresponding info
#rs = listSRAfile ( conversion[,4], sra_con, fileType = 'sra', srcType='fasp' )           #######display the download info
#rs = getSRAinfo ( conversion[,4], sra_con, sraType = 'sra' )           #######display more download info, including size

########download sra files#########
#ascpCMD <- 'ascp -QT -l 300m -i /psc/home/sunyidan/.aspera/connect/etc/aspera_id_dsa.putty'
#getSRAfile( conversion[,4], sra_con, fileType = 'sra', srcType='fasp', ascpCMD = ascpCMD )     ####use aspera to download
getSRAfile( file[,1], sra_con, fileType ='sra', srcType = "ftp" )   ####use ftp to download
#write.table(out,args[1],row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

