genecluster<-function(ratio,nct,G){
  #PCA first
  pca<-prcomp(ratio_pseudo,rank. = 2*nct) #use 2*nct
  # vari<-summary(pca)$importance[2,1:2*nct]
  # sum(vari)
  ratio_pca<-as.matrix(pca$x)
  # d_clust<-mclustBIC(ratio_pca,modelNames ="EEI",initialization = list(hcPairs = hc(ratio_pca,modelName = "EII", use = "VARS")))
  # d_clust<-mclustICL(ratio_pca,modelNames = c("EEI"),initialization = list(hcPairs = hc(ratio_pca,modelName = "EII", use = "VARS")))
  
  d_clust<-Mclust(ratio_pca,G=G,modelNames = "EEI",initialization = list(hcPairs = hc(ratio_pca,modelName = "EII", use = "VARS")))
  # summary(d_clust)
  # plot(d_clust)
  m.best <- dim(d_clust$z)[2]
  cat("model-based optimal number of clusters:", m.best, "\n")
  # plot(ratio_pca[,1],ratio_pca[,2],col=d_clust$classification)
  # mod<-MclustDR(d_clust)
  # plot(mod)
  return(d_clust$classification)
}
# dist<-dist(ratio_pca,method = "manhattan")
# nb<-cmdscale(dist)
# plot(nb[,1],nb[,2])
# km.res<-kmeans(ratio_pca,3,nstart = 25)
# km.res$cluster
