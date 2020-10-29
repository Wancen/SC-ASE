library("smurf")
library(dplyr)
library(emdbook)
library(mclust)
library(tibble)
library(pbapply)
library(factoextra)
library(aricode)
library(clue)
library(pheatmap)
n <- 20 # cells per cluster
frac <- .5 # fraction of high counts
ngene <- 200
nct<-8
# f <- factor(rep(1:4,each=n/4)) # two clusters
f <- factor(rep(1:nct,each=n))
mu1 <- 5
mu2 <- 100
nb.disp <- 1/100
k=5
set.seed(1)
T=20
sd<-0.01
################################################################################
## Define true cell type label according to below simulations####################
ctype_tru<-rbind(matrix(rep(c(1,1,2,2,3,3,4,4),each=20),ncol=nct),matrix(rep(c(1,1,2,2,3,3,3,3),each=20),ncol=nct),matrix(rep(c(1,1,1,1,1,1,1,1),each=40),ncol=nct),matrix(rep(c(1,1,2,2,2,2,2,2),each=20),ncol=nct),matrix(rep(c(1,1,2,2,3,3,4,4),each=20),ncol=nct),matrix(rep(c(1,1,1,1,1,1,1,1),each=80),ncol=nct))
colnames(ctype_tru)<-paste0("type",seq_len(nct))
ctype_wilcoxtru<-rbind(matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0),each=20),ncol=choose(nct,2)),
                       matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0),each=20),ncol=choose(nct,2)),
                       matrix(rep(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),each=40),ncol=choose(nct,2)),
                       matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),each=20),ncol=choose(nct,2)),
                       matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0),each=20),ncol=choose(nct,2)),
                       matrix(rep(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),each=80),ncol=choose(nct,2)))
colnames(ctype_wilcoxtru)<-paste0("pairwise-type",seq_len(choose(nct,2)))

## coefficient average MSE ############
msebin_mean<-matrix(NA,nrow = T,ncol = nct)
msebin2_mean<-matrix(NA,nrow = T,ncol = nct)
mselin_mean<-matrix(NA,nrow = T,ncol = nct)
mselin2_mean<-matrix(NA,nrow = T,ncol = nct)

## ARI for all iterations ###################
ARI_gene<-numeric(T)
ARI_lin<-numeric(T)
ARI_wilcoxon<-numeric(T)
ARI_bin<-numeric(T)

## Start for loop ########################
system.time(for (t in 1:T) {
  print(paste("Iteration",t))
  
## Similaute total counts ################################################
cts2 <- matrix(rep(c(rnbinom(ngene*n/2,mu=mu1,size=1/nb.disp),
                     rnbinom(ngene*n/2,mu=mu2,size=1/nb.disp)),8),ncol = 8*n)

theta <- runif(ngene,10,50)
theta2 <- runif(ngene,10,50)
colnames(cts2) <-paste0("cell",seq_len(ncol(cts2)))  

## Simulate allelic ratio pattern ########################
prob1 <- rnorm(ngene,.5,sd) # close to 0.5
prob1[1:20] <- rnorm(10,.9,sd)
prob1[21:40] <- rnorm(10,.75,sd)
prob1[81:100] <- rnorm(10,.3,sd)
prob1[101:120] <- rnorm(10,.3,sd)

prob2 <- rnorm(ngene,.5,sd) # close to 0.5
prob2[1:20] <- rnorm(10,0.8,sd)
prob2[21:40] <- rnorm(10,0.65,sd)
prob2[101:120] <- rnorm(10,0.35,sd)

prob3 <- rnorm(ngene,.5,sd) # close to 0.5
prob3[1:20] <- rnorm(10,0.75,sd)
prob3[101:120] <- rnorm(10,0.2,sd)

prob4 <- rnorm(ngene,.5,sd) # close to 0.5
# plot(prob1,prob2,col=gene_cluster$cluster)
# counts for allele 2
ase.cts3<- matrix(c(
  rbetabinom(prod(dim(cts2))/8,prob=prob1,size=cts2[,1:n],theta=rep(theta,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8, prob=prob1,size=cts2[,(n+1):(2*n)], theta=rep(theta2,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8,prob=prob2,size=cts2[,(2*n+1):(3*n)],theta=rep(theta,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8, prob=prob2,size=cts2[,(3*n+1):(4*n)], theta=rep(theta2,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8,prob=prob3,size=cts2[,(4*n+1):(5*n)],theta=rep(theta,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8, prob=prob3,size=cts2[,(5*n+1):(6*n)], theta=rep(theta2,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8,prob=prob4,size=cts2[,(6*n+1):(7*n)],theta=rep(theta,ncol(cts2)/8)),
  rbetabinom(prod(dim(cts2))/8, prob=prob4,size=cts2[,(7*n+1):(8*n)], theta=rep(theta2,ncol(cts2)/8))),
  ncol = 8*n)
ratio3<-(ase.cts3)/(cts2)
## pseudo allelic ratio for gene clustering ##################
ratio33<-(ase.cts3+1)/(cts2+2)
level<-c("type1", "type2", "type3","type4","type5","type6","type7","type8")
jpeg("simultaion_heatmap.jpeg",width = 5, height = 4,units = "in",res=300)
# anno_df <- data.frame(celltype=rep(level,each=n), row.names=colnames(cts2))
pheatmap(ratio33, cluster_rows = FALSE, cluster_cols = FALSE,
annotation_col=anno_df,show_colnames = F)
dev.off()
# 
# #######################################
# ## Hierarchical clustering ############
# ###############################################
gene_dist<-dist(ratio33,method = "manhattan")
# nbcluster<-cmdscale(gene_dist)
# x <- nbcluster[,1]
# y <- nbcluster[,2]
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
#      main="Metric MDS")
gene_hclust <- hclust(gene_dist,method = "ward.D2")
# plot(gene_hclust, labels = FALSE)
# abline(h = 6, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# cut genes into 5 group
gene_cluster<-cutree(gene_hclust, k = 5) %>%
  # turn the named vector into a tibble
  enframe(name="gene",value = "cluster")
table(gene_cluster$cluster)
gene_cluster[["true"]]<-matrix(c(rep(1,20),rep(2,20),rep(3,40),rep(4,20),rep(5,20),rep(3,80)),ncol = 1)

#######################################
## kmeans ############
###############################################
# set.seed(123)
####select optimal k##########################
# n_clust<-fviz_nbclust(ratio33, hcut, method = "gap_stat")
# n_clust<-n_clust$data
# max_cluster<-as.numeric(n_clust$clusters[which.max(n_clust$y)])
# km.res <- kmeans(ratio33, 5,nstart = 50)
# calinhara(ratio33,gene_cluster$cluster)
# 
# gene_cluster<-data.frame(gene=seq_len(ngene),cluster=km.res$cluster,true=matrix(c(rep(1,10),rep(2,10),rep(3,60),rep(4,10),rep(5,10)),ncol = 1))

ARI_gene[t]<-adjustedRandIndex(gene_cluster$cluster,gene_cluster$true)
print(ARI_gene[t])

ct_bin<-matrix(0,ncol=nct,nrow = ngene)
ct_lin<-matrix(0,ncol=nct,nrow = ngene)
ct_wilcoxon<-matrix(0,ncol = choose(nct,2),nrow = ngene)
msebin<-matrix(0,nrow=nct*n,ncol = ngene)
msebin2<-matrix(0,nrow=nct*n,ncol = ngene)
mselin<-matrix(0,nrow=nct*n,ncol = ngene)
mselin2<-matrix(0,nrow=nct*n,ncol = ngene)
########################################################################
## detect cell type according allelic ratio ############################
#####################################################
for (i in 1:5) {
print(paste("gene pattern",i))
poi<-which(gene_cluster$cluster==i)
truth<-t(ratio3[poi,])
# est<-c(rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),each=n))
# test<-(truth-est)^2

gene1_mean<-as.vector(ratio3[poi,])
poina<-which(!is.na(gene1_mean))
count1_mean<-as.vector(cts2[poi,])

rent<-tibble(ratio=gene1_mean,X=as.factor(rep(c(1,2,3,4,5,6,7,8),each=n*length(poi))),cts=count1_mean)
levels(rent$X) <- c("type1", "type2", "type3","type4","type5","type6","type7","type8")

## Gaussian likelihood #####################
formu <- ratio ~ p(X, pen = "gflasso",refcat = "type1")
fit2 <- glmsmurf(formula = formu, family=gaussian(), data = rent,
                                    pen.weights = "glm.stand", lambda = "cv1se.dev",
                                    control = list(lambda.length = 50L,k=10)) #pen.weights="glm" will fused coefficients more
                                                                               # and more consistent with truth sometimes 
coeflin<-coef(fit2)
coeflin2<-coef_reest(fit2)
coeflin[2:nct]<-coeflin[2:nct]+coeflin[1]
coeflin2[2:nct]<-coeflin2[2:nct]+coeflin2[1]
type_lin=match(coeflin2, unique(coeflin2))
ct_lin[poi,]<-rep(type_lin,each=length(poi))

## Wilcoxon Rank Sum test ###################################
res <- pairwise.wilcox.test(rent$ratio,rent$X, p.adjust.method ="BH")
adj<-as.data.frame(res$p.value)[lower.tri(res$p.value, diag = T)]
sigdif<-which(adj<0.05)
ct_wilcoxon[poi,sigdif]<-1
# a<-which(res$p.value<0.05)
# asso<-cbind(colnames(adj)[ceiling(a/3)],rownames(adj)[ifelse(a%%3==0,3,a%%3)])
# adjm<-matrix(TRUE,nrow = 4,ncol = 4)
# row.names(adjm)<-c("type1", "type2", "type3","type4")
# colnames(adjm)<-c("type1", "type2", "type3","type4")
# for (i in 1:nrow(asso)) {
#   adjm[which(row.names(adjm)==asso[i,1]),which(colnames(adjm)==asso[i,2])]<-FALSE
#   adjm[which(colnames(adjm)==asso[i,2]),which(row.names(adjm)==asso[i,1])]<-FALSE
# }
# gr<-graph_from_adjacency_matrix(adjm,mode = "undirected",diag = F)
# plot(gr)

## Binomial likelihood ###########################
fit = tryCatch({
  glmsmurf(formula = formu, family=binomial(link = "logit"), data = rent, weights = count1_mean[poina],
                  pen.weights = "glm.stand", lambda = "cv1se.dev" ,control = list(lambda.length = 50L,k=10)) #
}, error = function(error_condition) {
  message('Failed determining the maximum value of lambda, using equal weight instead')
   glmsmurf(formula = formu, family=binomial(link = "logit"), data = rent, weights = count1_mean[poina],
                  pen.weights = "eq", lambda = "cv1se.dev" ,control = list(lambda.length = 50L,k=10)) #
})
coefbin<-coef(fit)
coefbin[2:nct]<-coefbin[2:nct]+coefbin[1]
coefbin2<-coef_reest(fit)
coefbin2[2:nct]<-coefbin2[2:nct]+coefbin2[1]
type=match(coefbin, unique(coefbin))
ct_bin[poi,]<-rep(type,each=length(poi))

msebin[,poi]<-(truth-exp(coefbin)/(1 + exp(coefbin)))^2
msebin2[,poi]<-(truth-exp(coefbin2)/(1 + exp(coefbin2)))^2
mselin[,poi]<-(truth-coeflin)^2
mselin2[,poi]<-(truth-coeflin2)^2
# mean(test,na.rm = T)
}


res2 <- pbsapply(1:ngene, function(i) {AMI(ctype_tru[i,],ct_lin[i,])})
res2[which(is.nan(res2))]<-1
ARI_lin[t]<-mean(res2) #0.6514286
reswilcox <- pbsapply(1:ngene, function(i) {AMI(ctype_wilcoxtru[i,], ct_wilcoxon[i,])})
reswilcox[which(is.nan(reswilcox))]<-1
ARI_wilcoxon[t]<-mean(reswilcox) #0.6118919
res <- pbsapply(1:ngene, function(i) {AMI(ctype_tru[i,],ct_bin[i,])})
res[which(is.nan(res))]<-1
ARI_bin[t]<-mean(res)

## calculate cell type specific average MSE #########################
msebin_mean[t,]<-pbsapply(1:nct,function(i){mean(unlist(msebin[which(f==i),]),na.rm=T)})
msebin2_mean[t,]<-pbsapply(1:nct,function(i){mean(unlist(msebin2[which(f==i),]),na.rm=T)})
mselin_mean[t,]<-pbsapply(1:nct,function(i){mean(unlist(mselin[which(f==i),]),na.rm=T)})
mselin2_mean[t,]<-pbsapply(1:nct,function(i){mean(unlist(mselin2[which(f==i),]),na.rm=T)})
# res3 <- pbsapply(1:ngene, function(i) {adjustedRandIndex(ctype_tru[i,],ct_bin3[i,])})
# ARI_bin3[t]<-mean(res3)#0.72
})

#############################################################
## Compare coefficient relative efficiency ###################
box<-cbind(msebin_mean,msebin2_mean,mselin_mean,mselin2_mean)
order<-vector(length = 32)
for (i in 1:8) {
  order[(4*i-3):(4*i)]<-c(i,i+8,i+16,i+24)
}
boxplot(sqrt(box[,order]),ylab="Average MSE",ylim=c(0.2,0.26))
abline(h=0, col=2)

ARI<-data.frame(Gaussian=ARI_lin,Wilcoxon=ARI_wilcoxon,Binomial=ARI_bin)
AMI<-ARI %>% gather(key=method,value=AMI,Gaussian:Binomial)

jpeg("violin.jpeg",width = 5, height = 4,units = "in",res=300) 
p <- ggplot(AMI, aes(x=method, y=AMI,fill=method)) + 
  geom_violin()+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ theme_classic()+labs(title="Violin Plot of AMI  by methods",x="Methods", y = "Adjusted Mutual Information")
p+ geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()



