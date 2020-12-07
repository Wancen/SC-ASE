library(mclust)
library(dynamicTreeCut)
genecluster<-function(ratio,nct,G=c(4,8,12,16,20),method="mclust",...){
  #PCA first
  pca<-prcomp(ratio,rank. = 2*nct) #use 2*nct
  # vari<-summary(pca)$importance[2,1:2*nct]
  # sum(vari)
  ratio_pca<-as.matrix(pca$x)
  
  if(method=="mclust"){
    d_clust<-Mclust(ratio_pca,G=G,modelNames = "EII",initialization = list(hcPairs = hc(ratio_pca,modelName = "EII", use = "VARS")),...)
    # summary(d_clust)
    # plot(d_clust)
    m.best <- dim(d_clust$z)[2]
    cat("model-based optimal number of clusters:", m.best, "\n")
    plot(ratio_pca[,1],ratio_pca[,2],col=d_clust$classification)
    # mod<-MclustDR(d_clust)
    # plot(mod)
    return(d_clust$classification)
  }
  if(method=="hierarchical"){
    my.dist <- dist(ratio_pca,method = "manhattan")
    my.tree <- hclust(my.dist, method="ward.D2")
    my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
    plot(ratio_pca[,1],ratio_pca[,2],col=my.clusters)
    return(my.clusters)
  }
}
