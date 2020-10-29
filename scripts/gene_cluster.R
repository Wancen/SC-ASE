setwd("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1")
load(file="gene_cluster.rda")
library(pheatmap)
library(factoextra)
library(dplyr)
library(tibble)
library(mclust)
nct=10
### gene clustering ############################
pca<-prcomp(ratio_psedo,rank. = 20)
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
pheatmap(ratio_psedo[gene_feat_mclust,], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = anno_df,show_colnames = F,show_rownames = F)

## Hierarchical clustering ##########################
start_time <- Sys.time()
gene_dist<-dist(ratio_psedo,method = "manhattan")
# check linage method to use #######################
# gene_hclust1 <- hclust(gene_dist, method = "complete")
gene_hclust <- hclust(gene_dist, method = "ward.D2")
jpeg(file="Dendrogram.jpeg",width = 5, height = 5,units = "in",res=450)
plot(gene_hclust,labels=FALSE)
abline(h = 180, col = "red", lwd = 1)
dev.off()

end_time <- Sys.time()
end_time - start_time

# c2=cophenetic(gene_hclust)
# c3=cophenetic(gene_hclust2)
# cor(gene_dist,c2)
# cor(gene_dist,c3)
abline(h = 180, col = "brown", lwd = 1) # add horizontal line to illustrate cutting dendrogram

# cut genes into 2 group
gene_cluster<-cutree(gene_hclust, h=180) %>%
  # turn the named vector into a tibble
  enframe(name="gene",value = "cluster")
table(gene_cluster$cluster)
gene_feat<-which(gene_cluster$cluster==7)

# jpeg(file="heatmap_183gene.jpeg",width = 7, height = 5,units = "in",res=450)
pheatmap(ratio_psedo[gene_feat,], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = anno_df,show_colnames = F,show_rownames = F)


