load("/proj/milovelab/mu/data/postprocess.rda")
source("../simulation/cluster.R")
source("../simulation/fusedlasso.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")


BiocManager::install("apeglm/1.12.0")
library(tidyverse)
library(dynamicTreeCut)
library(pheatmap)
library(boot)
library(smurf)
library(doParallel)
library(pbapply)
library(clue)
library(fastmatch)
nct=10
ncores=10
total1 <- total[keep_feature, ]
total_CellQC<-total1[,cellQC]
total_fil<-total_CellQC[keep_feature2,]
total_QCed<-total_fil[check2,]
ratio<-((cast_uq[keep_feature,cellQC])[keep_feature2,])[check2,]/total_QCed
ratio_psedo<-(((cast_uq[keep_feature,cellQC])[keep_feature2,])[check2,]+2)/(total_QCed+4)

cts<-apply(total_QCed, 1,mean, na.rm=T)
x <- matrix(1, ncol=1, nrow=nrow(total_QCed))
theta.hat <- 100
param <- cbind(theta.hat, cts)
fit.mle <- apeglm(Y=ratio, x=x, log.lik=NULL, param=param, no.shrink=TRUE, log.link=FALSE, method="betabinC")
theta.hat <- bbEstDisp(success=ratio, size=total_QCed, x=x, beta=fit.mle$map, minDisp=.01, maxDisp=5000)