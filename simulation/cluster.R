genecluster<-function(ratio,nct,G){
  #PCA first
  pca<-prcomp(ratio,rank. = 2*nct) #use 2*nct
  # vari<-summary(pca)$importance[2,1:2*nct]
  # sum(vari)
  ratio_pca<-as.matrix(pca$x)
  d_clust<-Mclust(ratio_pca,modelNames =c("EII"), G=G,initialization = list(hcPairs = hc(ratio_pca, use = "SVD")),prior = priorControl())
  m.best <- dim(d_clust$z)[2]
  cat("model-based optimal number of clusters:", m.best, "\n")
  return(d_clust$classification)
}
