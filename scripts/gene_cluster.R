setwd("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1")
load(file="gene_cluster.rda")
library(pheatmap)
library(factoextra)
library(dplyr)
library(tibble)
library(mclust)
nct=10
### gene clustering ############################
pca<-prcomp(ratio21,rank. = 20)
vari<-summary(pca)$importance[2,1:20]
sum(vari)
ratio_pca<-as.matrix(pca$x)

## decide optimal number of clusters ###########
mclust.options(emModelNames=c("VVE","VVI","VEI","VEE"))
mclust.options(hcModelName = "VVV")
system.time(d_clust<-Mclust(ratio_pca, G=c(8,10,12,14),initialization = list(hcPairs = hc(ratio_pca, use = "SVD"))))
d_clust$BIC
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
summary(d_clust)
gene_cluster_mclust<-data.frame(gene=rownames(ratio_pca),cluster=d_clust$classification)
table(gene_cluster_mclust$cluster)
gene_feat_mclust<-which(gene_cluster_mclust$cluster==7)
pheatmap(ratio21[gene_feat_mclust,], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = anno_df,show_colnames = F,show_rownames = F)
# check each gene pattern histogram ####################
hist(as.numeric(unlist(ratio21[gene_feat_km,])),main="Distribution of allelic ratio",xlab = "allelic ratio")
## kmeans #################
## use m.best as input for estimated k##############
km.res <- kmeans(ratio_pca, m.best,nstart = 10)
gene_cluster_km1<-data.frame(gene=rownames(ratio_pca),cluster=km.res$cluster)
table(gene_cluster_km1$cluster)
gene_feat_km<-which(gene_cluster_km1$cluster==1)

pheatmap(ratio21[gene_feat_km,], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = anno_df,show_colnames = F,show_rownames = F)

# save(ratio21, anno_df,file = "gene_cluster.rda")
