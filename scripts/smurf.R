library("smurf")
library(dplyr)
library(emdbook)
library(mclust)
library(tibble)
library(pbapply)
library(factoextra)
library(aricode)
n <- 20 # cells per cluster
frac <- .5 # fraction of high counts
ngene <- 100
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
ctype_tru<-rbind(matrix(rep(c(1,1,2,2,3,3,4,4),each=10),ncol=nct),matrix(rep(c(1,1,2,2,3,3,3,3),each=10),ncol=nct),matrix(rep(c(1,1,1,1,1,1,1,1),each=60),ncol=nct),matrix(rep(c(1,1,2,2,2,2,2,2),each=10),ncol=nct),matrix(rep(c(1,1,2,2,3,3,4,4),each=10),ncol=nct))
colnames(ctype_tru)<-paste0("type",seq_len(nct))
ctype_wilcoxtru<-rbind(matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0),each=10),ncol=choose(nct,2)),
                       matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0),each=10),ncol=choose(nct,2)),
                       matrix(rep(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),each=60),ncol=choose(nct,2)),
                       matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),each=10),ncol=choose(nct,2)),
                       matrix(rep(c(0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,0),each=10),ncol=choose(nct,2)))
colnames(ctype_wilcoxtru)<-paste0("pairwise-type",seq_len(choose(nct,2)))

## Extract coefficients #######################
bin<-matrix(0,nrow = k*T,ncol = nct)
# bin3<-matrix(0,nrow = k*T,ncol = nct)
lin<-matrix(0,nrow = k*T,ncol = nct)
truth<-matrix(0,nrow = k*T,ncol = nct)

## ARI for all iterations ###################
ARI_gene<-numeric(T)
ARI_lin<-numeric(T)
ARI_wilcoxon<-numeric(T)
ARI_bin<-numeric(T)
# ARI_bin3<-numeric(T)

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
prob1[1:10] <- rnorm(10,.9,sd)
prob1[11:20] <- rnorm(10,.75,sd)
prob1[81:90] <- rnorm(10,.3,sd)
prob1[91:100] <- rnorm(10,.3,sd)

prob2 <- rnorm(ngene,.5,sd) # close to 0.5
prob2[1:10] <- rnorm(10,0.8,sd)
prob2[11:20] <- rnorm(10,0.65,sd)
prob2[91:100] <- rnorm(10,0.35,sd)

prob3 <- rnorm(ngene,.5,sd) # close to 0.5
prob3[1:10] <- rnorm(10,0.75,sd)
prob3[91:100] <- rnorm(10,0.2,sd)

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
anno_df <- data.frame(f, row.names=colnames(cts2))
pheatmap(ratio33, cluster_rows = FALSE, cluster_cols = FALSE,
annotation_col=anno_df,show_colnames = F)

## decide optimal number of clusters ###########
d_clust <- Mclust(as.matrix(ratio33), G=c(2,4,5,6))
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")

#######################################
## Hierarchical clustering ############
###############################################
gene_dist<-dist(ratio33,method = "euclidean")
nbcluster<-cmdscale(gene_dist)
x <- nbcluster[,1]
y <- nbcluster[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS")
gene_hclust <- hclust(gene_dist)
# plot(gene_hclust, labels = FALSE)
# abline(h = 6, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# cut genes into 5 group
gene_cluster<-cutree(gene_hclust1, k = 5) %>%
  # turn the named vector into a tibble
  enframe(name="gene",value = "cluster")
gene_cluster[["true"]]<-matrix(c(rep(1,10),rep(2,10),rep(3,60),rep(4,10),rep(5,10)),ncol = 1)

#######################################
## kmeans ############
###############################################
# set.seed(123)
####select optimal k##########################
# n_clust<-fviz_nbclust(ratio33, hcut, method = "gap_stat")
# n_clust<-n_clust$data
# max_cluster<-as.numeric(n_clust$clusters[which.max(n_clust$y)])
km.res <- kmeans(ratio33, 5,nstart = 15)

gene_cluster<-data.frame(gene=seq_len(ngene),cluster=km.res$cluster,true=matrix(c(rep(1,10),rep(2,10),rep(3,60),rep(4,10),rep(5,10)),ncol = 1))

ARI_gene[t]<-adjustedRandIndex(gene_cluster$cluster,gene_cluster$true)
print(ARI_gene[t])

ct_bin<-matrix(0,ncol=nct,nrow = ngene)
ct_lin<-matrix(0,ncol=nct,nrow = ngene)
ct_bin3<-matrix(0,ncol=nct,nrow = ngene)
ct_wilcoxon<-matrix(0,ncol = choose(nct,2),nrow = ngene)

########################################################################
## detect cell type according allelic ratio ############################
#####################################################
for (i in 1:5) {
print(paste("gene pattern",i))
poi<-which(gene_cluster$cluster==i)
# gene1_mean<-if(length(poi)==1) ratio3[poi,]else colMeans(ratio3[poi,],na.rm = T)
gene1_mean<-as.vector(ratio3[poi,])
poina<-which(!is.na(gene1_mean))
count1_mean<-as.vector(cts2[poi,])
# truth[(5*t-5+i),]<-sapply(1:nct, function(j){mean(gene1_mean[(n*j*length(poi)-n*length(poi)+1):(n*j*length(poi))],na.rm=T)})

rent<-tibble(ratio=gene1_mean,X=as.factor(rep(c(1,2,3,4,5,6,7,8),each=n*length(poi))),cts=count1_mean)
levels(rent$X) <- c("type1", "type2", "type3","type4","type5","type6","type7","type8")

## Gaussian likelihood #####################
formu <- (ratio-0.5) ~ p(X, pen = "gflasso",refcat = "type8")
fit2 <- glmsmurf(formula = formu, family=gaussian(), data = rent,
                                    pen.weights = "glm.stand", lambda = "cv1se.mse",lambda1 = 0.02,
                                    control = list(lambda.length = 50L,k=10),pen.weights.return=T) #pen.weights="glm" will fused coefficients more
                                                                               # and more consistent with truth sometimes 
print(coef(fit2))
fit2$lambda
print(coef_reest(fit2))
coeflin<-coef_reest(fit2)
coeflin[2:(nct+1)]<-coeflin[2:(nct+1)]+coeflin[1]
type_lin=match(coeflin, unique(coeflin))
ct_lin[poi,]<-rep(type_lin,each=length(poi))
plot_lambda(fit2)
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
fit <- glmsmurf(formula = formu, family=binomial(link = "logit"), data = rent, weights = count1_mean[poina],lambda1 = 0.002,
                       pen.weights = "glm.stand", lambda = "cv1se.mse" ,control = list(lambda.length = 50L,k=10)) #

print(coef(fit))
coefbin<-coef(fit)
coefbin[2:nct]<-coefbin[2:nct]+coefbin[1]
type=match(coefbin, unique(coefbin))
ct_bin[poi,]<-rep(type,each=length(poi))

bin[(5*t-5+i),]=sapply(1:nct, function(j){exp(coefbin[j])/(1 + exp(coefbin[j]))})
# bin3[(5*t-5+i),]=sapply(1:nct, function(j){exp(coefbin3[j])/(1 + exp(coefbin3[j]))})
lin[(5*t-5+i),]=sapply(1:nct, function(j){coeflin[j]})
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
# res3 <- pbsapply(1:ngene, function(i) {adjustedRandIndex(ctype_tru[i,],ct_bin3[i,])})
# ARI_bin3[t]<-mean(res3)#0.72
})

#############################################################
## Compare coefficient relative efficiency ###################
box<-data.frame(bindiff1=bin[,1]-truth[,1],
                lindiff1=lin[,1]-truth[,1],
                bindiff2=bin[,2]-truth[,2],
                lindiff2=lin[,2]-truth[,2],
                bindiff3=bin[,3]-truth[,3],
                lindiff3=lin[,3]-truth[,3],
                bindiff4=bin[,4]-truth[,4],
                lindiff4=lin[,4]-truth[,4],
                bindiff1=bin[,5]-truth[,5],
                lindiff1=lin[,5]-truth[,5],
                bindiff2=bin[,6]-truth[,6],
                lindiff2=lin[,6]-truth[,6],
                bindiff3=bin[,7]-truth[,7],
                lindiff3=lin[,7]-truth[,7],
                bindiff4=bin[,8]-truth[,8],
                lindiff4=lin[,8]-truth[,8])
boxplot(box,ylab="estimated-truth",ylim=c(-0.1,0.1))
abline(h=0, col=2)

ARI<-data.frame(Gaussian=ARI_lin,Wilcoxon=ARI_wilcoxon,Binomial=ARI_bin)
pdf("rplot.pdf") 
boxplot(ARI,ylab="Adjusted mutual information")
boxplot(ARI_gene,ylab="ARI for gene cluster")
dev.off()

res2 <- pbsapply(1:ngene, function(i) {adjustedRandIndex(ctype_tru[i,],ct_lin[i,])})
ARI_lin[t]<-mean(res2) #0.6514286
reswilcox <- pbsapply(1:ngene, function(i) {MARI(ctype_wilcoxtru[i,], ct_wilcoxon[i,])})
ARI_wilcoxon[t]<-mean(reswilcox) #0.6118919
res <- pbsapply(1:ngene, function(i) {MARI(ctype_tru[i,],ct_bin[i,])})
ARI_bin[t]<-mean(res) 


### Do parallel bagging ############################
library(doParallel)  # for parallel backend to foreach
library(foreach)     # for parallel processing with for loops

# Create a parallel socket cluster
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

# Fit model in parallel 
system.time(fit22 <- foreach(
  icount(5), #run 5 iterations
  .packages = c("smurf","tibble"), 
  .combine = cbind
) %dopar% {
  # bootstrap copy of training data
  poi<-which(gene_cluster$cluster==2)
  n1<-length(poi)
  df_s =poi[sample(1:n1,size=n1/2,replace=F)]
  gene1_mean<-as.vector(ratio3[df_s,])
  poina<-which(!is.na(gene1_mean))
  count1_mean<-as.vector(cts2[df_s,])
  rent<-tibble(ratio=gene1_mean,X=as.factor(rep(c(1,2,3,4,5,6,7,8),each=n*length(df_s))),cts=count1_mean)
  levels(rent$X) <- c("type1", "type2", "type3","type4","type5","type6","type7","type8")
  
  # fit fused lasso model to bootstrap copy
  formu <- ratio ~ p(X, pen = "gflasso")
  fit <- glmsmurf(formula = formu, family=binomial(link = "logit"), data = rent, weights = count1_mean[poina],
                  pen.weights = "glm.stand", lambda = "cv1se.mse" ,control = list(lambda.length = 50L,k=10)) #
  coefbin<-coef(fit)
  coefbin[2:nct]<-coefbin[2:nct]+coefbin[1]
  match(coefbin, unique(coefbin))
})
