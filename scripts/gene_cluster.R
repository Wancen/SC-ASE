setwd("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1")
load(file="gene_cluster.rda")
library(pheatmap)
library(factoextra)
library(dplyr)
library(tibble)
library(mclust)
nct=10
### gene clustering ############################
pca<-prcomp(ratio_pseudo,rank. = 20)
vari<-summary(pca)$importance[2,1:20]
sum(vari)
ratio_pca<-as.matrix(pca$x)

## decide optimal number of clusters ###########

system.time(d_clust<-Mclust(ratio_pca,modelNames =c("EII"), G=c(8,10,12,14,16,17,18,20),initialization = list(hcPairs = hc(ratio_pca, use = "SVD"))))
d_clust$BIC
d_clust$parameters
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
summary(d_clust)
gene_cluster_mclust<-data.frame(gene=rownames(ratio_pca),cluster=d_clust$classification)
table(gene_cluster_mclust$cluster)
plot(ratio_pca[,1],ratio_pca[,2],col=d_clust$classification)


