args <- commandArgs(trailingOnly = TRUE)
com1 <- read.delim("/psc/bioinformatics/sunyd/lncrna/com1_top50.feature", head = TRUE, quote="\"", dec=".", row.names = 1, as.is = FALSE)      ####model built dataset
#com2 <- read.delim(args[2], head = TRUE, quote="\"", dec=".", row.names = 1, as.is = FALSE)      ##validation dataset
com2 <- read.delim(args[1], head = TRUE, dec=".", row.names = 1, as.is = FALSE)      ##noncode validation dataset
com3 <- args[2]        ##output file name

#com1 <- read.delim("/psc/bioinformatics/sunyd/lncrna/combinetest", head = TRUE, quote="\"", dec=".", row.names = 1,  as.is = FALSE)
#com1 <- read.delim("E:/machine_learning/combinetest", head = TRUE, quote="\"", dec=".", row.names = 1,  as.is = FALSE)
com1[is.na(com1)] <- 0
com2[is.na(com2)] <- 0
library("caret")
library("ROCR")

##################preprocessing############################
#hist(train$sequence_length,main="",xlab="sequence length")
#normalization <- preProcess(com1[,-1],method = c("center", "scale"))
#com1 <- predict(normalization, com1[,-1])
#trainstd <- predict(preobj,train[,-1])
#testsdd <- predict(preobj,test[,-1])



####################data splicing to make train and test data######################
set.seed(40000)
intrain <- createDataPartition(com1$class,p=0.75,list=FALSE)
#train <- com1[intrain,]
#test <- com1[-intrain,]
#fold <- createFolds(y=com1$class,k=10,list=TRUE,returnTrain=TRUE)
#sapply(fold,length)




#################train a model with cross validation##################
model_training <- function(x,mod){
	set.seed(40000)
	controltrain <- trainControl(method = "cv",number=10,p=0.8,classProbs = TRUE)
	model <- train(class ~.,data=x,method=mod,preProcess=c("center", "scale"),trControl=controltrain)
	return(model)
	}
	
###############performance measure of the predictions###############
performance_measure <- function(model,testset,filename){
	prediction1 <- predict(model, newdata=testset,type = "prob")
	prediction2 <- predict(model, newdata=testset)
	result <- cbind(prediction1[,1,drop=FALSE],as.matrix(prediction2))
	#colnames(result) <- parameter
	colnames(result) <- c("lncRNA_potential","Type")
	result <- result[(result[,2] == "lncRNA"),]
	write.table(result, file = filename, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
	return(result)
}


###########evaluating com1 performance on self test data###############"E:/machine_learning/caret_all"
train <- com1[intrain,]
test <- com1[-intrain,]

rf_model <- model_training(train,"rf")           ###train the model
performance_measure(rf_model,com2,com3)     ##  --means implicit feature selection
#eval_non <- evaluationall(rf_model,vali_noncode,"Noncode_v4_RF",com3)     ##  --means implicit feature selection
#eval_zea <- evaluationall(rf_model,vali_zea,"Zea_mays_RF",com3 )     ##  --means implicit feature selection






