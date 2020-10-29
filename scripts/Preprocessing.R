setwd("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/project1")
load("Deng's.RData")
library(pheatmap)
library(factoextra)
library(dplyr)
library(tibble)
library(smurf)
library(M3C)
library(pbapply)
library(mclust)
library(parallel)
nct=10

## Preprocess ############################################
#C57BL/6J(mother) x CAST/EiJ(father)
B6<-read.table("B6.txt",skip = 1)
CAST<-read.table("CAST.txt",skip = 1)
genename<-read.table("genename.txt",skip = 1)
sampleinfo<-scan("sample.txt", character(), quote = "")
celltype<-sub("_.*", "", sampleinfo)

## Check categories of cells, remove unuseful cells and outlier cells indicated in Deng's Supplyment###################
ratio1<-(CAST)/(CAST+B6)

celltype[grep("zy",celltype)]="zy";
sampleinfo[grep("8cell_2pooled", sampleinfo)]="8cell_nd";
sampleinfo[grep("8cell_split", sampleinfo)]="8cell_nd";
sampleinfo[grep("16cell_2pooled", sampleinfo)]="16cell_nd";
sampleinfo[grep("16cell_split", sampleinfo)]="16cell_nd";
sampleinfo[grep("8cell_1", sampleinfo)]="8cell_nd";
sampleinfo[grep("8cell_5", sampleinfo)]="8cell_nd";
sampleinfo[grep("lateblast_2", sampleinfo)]="lateblast_nd";
midblast<-colMeans(ratio1[,grep("midblast", sampleinfo)],na.rm = T)
hist(midblast)
sampleinfo[order(midblast,decreasing = T)[1]]="midblast_nd";
eicell<-colMeans(ratio1[,grep("8cell", sampleinfo)],na.rm = T)
hist(eicell)
sampleinfo[order(eicell)[1]]="8cell_nd";
sixtcell<-colMeans(ratio1[,grep("16cell", sampleinfo)],na.rm = T)
hist(sixtcell)
sampleinfo[order(sixtcell,decreasing = T)[1]]="16cell_nd";
indices_not_reqd <- which(celltype=="BXC"   | celltype=="C57twocell" | celltype=="fibroblast" | sampleinfo =="8cell_nd"| sampleinfo =="lateblast_nd" |sampleinfo =="midblast_nd" | sampleinfo == "16cell_nd")
celltype<-celltype[-indices_not_reqd]
sampleinfo<-sampleinfo[-indices_not_reqd]
table(celltype)

#check whether all gene name are same across all samples and return True if it is
gene<-apply(genename, 1,unique)
length(gene)==nrow(genename)

#return the position where have duplicate gene name
genedui<-which(duplicated(gene) | duplicated(gene, fromLast = TRUE))

## Change the order based on development order ############################
cell_meta_unique <- c("zy","early2cell","mid2cell","late2cell","4cell","8cell","16cell","earlyblast","midblast","lateblast") ;
order_of_development <- order(match(celltype,cell_meta_unique))
sampleinfo<-sampleinfo[order_of_development]
celltype<- celltype[order_of_development]

## remove unrelated genes and cells and order count matrix ####
b6_uq<-B6[-genedui,-indices_not_reqd]
rownames(b6_uq)<-gene[-genedui]
b6_uq<-b6_uq[,order_of_development]
colnames(b6_uq)<-sampleinfo
cast_uq<-CAST[-genedui,-indices_not_reqd]
rownames(cast_uq)<-gene[-genedui]
cast_uq<-cast_uq[,order_of_development]
colnames(cast_uq)<-sampleinfo
total<-b6_uq+cast_uq


## Define embryo level for futher research ##############
cell_embryo<-sub("-.*", "", sampleinfo)
cell_embryo[grep("zy",cell_embryo)]="zy";
cell_embryo_unique<-unique(cell_embryo)

## create annotation df for heatmap ##########################
anno_df = data.frame(CellType = factor(celltype,levels = cell_meta_unique),row.names = sampleinfo)
anno_df2 = data.frame(CellType = factor(cell_embryo,levels = cell_embryo_unique),row.names = sampleinfo)

## define cell types except 3 mainly expressed cell types for filtering ###############
celltype_short<-celltype
indices_specific<-which(celltype_short%in%c("early2cell","late2cell","zy"))
celltype_short[-indices_specific]<-"others"

## filter out low expressed genes ###################
# take 1min27s
# genesum<-pbsapply(1:nrow(total),function(i){length(which(total[i,]>=5))}) 

# take 20.68s
no_cores <- detectCores()
clust <- makeCluster(no_cores) #This line will take time
clusterExport(clust, "total")
#The parallel version of lapply() is parLapply() and needs an additional cluster argument.
system.time(genesum<-parSapply(clust, 1:nrow(total),function(i){length(which(total[i,]>=5))}))
stopCluster(clust)

check<-which(genesum>ncol(total)*0.25)
total_fil<-total[check,]
# pheatmap(total[check,], cluster_rows = FALSE, cluster_cols = FALSE,
         # annotation_col = anno_df,show_colnames = F,show_rownames = F)


## filter out genes mainly expressed in three cell types #################
celltype_prop<- pbsapply(1:nrow(total_fil),function(i){tapply(unlist(total_fil[i,]),as.factor(celltype_short), FUN=sum)})
other_prop<-celltype_prop["others",]/colSums(celltype_prop)
h2<-hist(other_prop,main="7514 genes",
        xlab="prop of other cell types")
text(h2$mids,h2$counts,labels=h2$counts, adj=c(0.5, -0.5))
abline(v=0.8,col="red")
check2<-which(other_prop>0.8)
total_fil2<-total_fil[check2,]

totalexp<-rowSums(total_fil2)
ratio<-(cast_uq[check,])[check2,]/total_fil2
ratio_psedo<-((cast_uq[check,])[check2,]+2)/(total_fil2+4)

### Dimension decution ############################
# pca<-prcomp(ratio_psedo,rank. = 10)
# vari<-summary(pca)$importance[2,1:20]
# sum(vari)
# ratio_pca<-as.matrix(pca$x)


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
jpeg(file="heatmap_183gene.jpeg",width = 7, height = 5,units = "in",res=450)
pheatmap(ratio_psedo[gene_feat,], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = anno_df,show_colnames = F,show_rownames = F)
dev.off()
##look at embro level ###
# ratio3_sub<-ratio3[gene_feat,]
# total3_sub<-total_fil3[gene_feat,]


save(ratio_psedo,total_fil2, anno_df,file = "gene_cluster.rda")
save(totalexp,gene_feat,ratio,total_fil2,cell_meta_unique,celltype, anno_df,file = "gene_modeling.rda")
