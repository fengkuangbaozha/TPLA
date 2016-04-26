mydata <- d
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata,
									                                          centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
	      ylab="Within groups sum of squares")


read.table(file="RS.diff2.seq.1",header=T,row.names=1) ->de
library(RColorBrewer)
library(pheatmap)
set.seed(1)
pheatmap(de,kmeans_k=6,scale="row")->dea
((as.matrix(de[names((sort(dea$kmeans$cluster))),]))) -> dea1
write.table(file="cluster_arde.diff1.txt",cbind(sort(dea$kmeans$cluster),de[names((sort(dea$kmeans$cluster))),]),quote=F,sep="\t" )
table(dea$kmeans$cluster)
tiff("Cluster_RS.tiff",height=480,width=1000)
colorRampPalette(c("green","black","red"))(100) ->col1
pheatmap(t(dea1),scale="column",color=col1,cluster_rows = F, cluster_cols = F,show_colnames=F,show_rownames=T,fontsize = 15)
dev.off()
