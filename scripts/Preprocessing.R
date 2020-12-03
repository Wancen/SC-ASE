setwd("/proj/milovelab/mu/")
load("./data/Deng's.RData")
library(tidyverse)
library(scater)
library(pheatmap)
library(factoextra)
library(pbapply)
library(mclust)
library(parallel)
nct=10

## Preprocess ############################################
#C57BL/6J(mother) x CAST/EiJ(father)
B6<-read_delim("./data/B6.txt",delim =" ", skip = 1,col_names = F)
CAST<-read_delim("./data/CAST.txt",skip = 1,delim=" ",col_names = F)
genename<-read_delim("./data/genename.txt",skip = 1,delim=" ",col_names = F)
sampleinfo<-scan("./data/sample.txt", character(), quote = "")
celltype<-sub("_.*", "", sampleinfo)

#check whether all gene name are same across all samples and return True if it is
gene<-apply(genename, 1,unique)
length(gene)==nrow(genename)

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
hist(midblast) #one extremely high midblast
sampleinfo[order(midblast,decreasing = T)[1]]="midblast_nd";
eicell<-colMeans(ratio1[,grep("8cell", sampleinfo)],na.rm = T)
hist(eicell) #one extremely low 8cell
sampleinfo[order(eicell)[1]]="8cell_nd";
sixtcell<-colMeans(ratio1[,grep("16cell", sampleinfo)],na.rm = T)
hist(sixtcell) #one extremely high 16cell
sampleinfo[order(sixtcell,decreasing = T)[1]]="16cell_nd";
indices_not_reqd <- which(celltype=="BXC"   | celltype=="C57twocell" | celltype=="fibroblast" | sampleinfo =="8cell_nd"| sampleinfo =="lateblast_nd" |sampleinfo =="midblast_nd" | sampleinfo == "16cell_nd")
celltype<-celltype[-indices_not_reqd]
sampleinfo<-sampleinfo[-indices_not_reqd]
table(celltype)

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
anno_df = data.frame(CellType = factor(celltype,levels = cell_meta_unique),row.names = sampleinfo) #cell type level
anno_df2 = data.frame(CellType = factor(cell_embryo,levels = cell_embryo_unique),row.names = sampleinfo) #embryo level

## define cell types except 3 mainly expressed cell types for filtering ###############
celltype_short<-celltype
indices_specific<-which(celltype_short%in%c("early2cell","late2cell","zy"))
celltype_short[-indices_specific]<-"others"

## filter out low expressed genes (express 5 reads in at least 25% cells)###################
keep_feature <- nexprs(as.matrix(total), byrow=TRUE,  detection_limit=5) >= ncol(total)*0.25
table(keep_feature)
total_fil<-total[keep_feature,]

## filter out genes mainly expressed in three cell types #################
celltype_prop<- pbsapply(1:nrow(total_fil),function(i){tapply(unlist(total_fil[i,]),as.factor(celltype_short), FUN=sum)})
other_prop<-celltype_prop["others",]/colSums(celltype_prop)
h2<-hist(other_prop,main="7514 genes",
        xlab="prop of other cell types")
text(h2$mids,h2$counts,labels=h2$counts, adj=c(0.5, -0.5))
abline(v=0.8,col="red")
check2<-which(other_prop>0.8)
total_QCed<-total_fil[check2,]

ratio<-(cast_uq[keep_feature,])[check2,]/total_QCed
ratio_psedo<-((cast_uq[keep_feature,])[check2,]+2)/(total_QCed+4)

save(ratio,ratio_psedo,total_QCed, cell_meta_unique,celltype,anno_df,file = "./data/postprocess.rda")
