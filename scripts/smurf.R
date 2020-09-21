library("smurf")
library(dplyr)
library(igraph)
library(aricode)
library(emdbook)
library(mclust)
library(tibble)
n <- 20 # cells per cluster
frac <- .5 # fraction of high counts
ngene <- 100
nct<-4
# f <- factor(rep(1:4,each=n/4)) # two clusters
f <- factor(rep(1:nct,each=n))
mu1 <- 3
mu2 <- 100
nb.disp <- 1/100
k=5
set.seed(1)
T=20
################################################################################
## Define true cell type label according to below simulations####################
ctype_tru<-matrix(rep(c(1,2,1,2),each=ngene),ncol=nct)
colnames(ctype_tru)<-paste0("type",seq_len(nct))
ctype_wilcoxtru<-matrix(rep(c(1,0,1,1,0,1),each=ngene),ncol=choose(nct,2))
colnames(ctype_wilcoxtru)<-paste0("pairwise-type",seq_len(choose(nct,2)))

## Extract coefficients #######################
bin<-matrix(0,nrow = k*T,ncol = nct)
lin<-matrix(0,nrow = k*T,ncol = nct)
truth<-matrix(0,nrow = k*T,ncol = nct)

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
                     rnbinom(ngene*n/2,mu=mu2,size=1/nb.disp)),4),ncol = 4*n)

theta <- runif(ngene,1,50)
theta2 <- runif(ngene,1,50)
colnames(cts2) <-paste0("cell",seq_len(ncol(cts2)))  

## Simulate allelic ratio pattern ########################
prob1 <- rnorm(ngene,.5,.05) # close to 0.5
prob1[1:10] <- rnorm(10,.7,.05)
prob1[11:20] <- rnorm(10,.3,.05)

prob2 <- rnorm(ngene,.55,.05) # close to 0.5
prob2[51:60] <- rnorm(10,.7,.05)
prob2[61:70] <- rnorm(10,.3,.05)

## Simulate allelic count ###############################
ase.cts3 <- matrix(c(
  rbetabinom(prod(dim(cts2))/4,prob=prob1,size=cts2[,1:n],theta=rep(theta,ncol(cts2)/4)),
  rbetabinom(prod(dim(cts2))/4, prob=prob2,size=cts2[,(n+1):(2*n)], theta=rep(theta,ncol(cts2)/4)),
  rbetabinom(prod(dim(cts2))/4,prob=prob1,size=cts2[,(2*n+1):(3*n)],theta=rep(theta2,ncol(cts2)/4)),
  rbetabinom(prod(dim(cts2))/4, prob=prob2,size=cts2[,(3*n+1):(4*n)], theta=rep(theta2,ncol(cts2)/4))),ncol = 4*n)
colnames(ase.cts3) <- paste0("cell",seq_len(ncol(cts2)))  

ratio3<-(ase.cts3)/(cts2)
## pseudo allelic ratio for gene clustering ##################
ratio33<-(ase.cts3+1)/(cts2+2)
# ctype_tru<-rbind(matrix(rep(c(1,2,1,2),each=20),ncol=nct),matrix(rep(c(1,1,1,1),each=30),ncol=nct),matrix(rep(c(1,2,1,2),each=20),ncol=nct),matrix(rep(c(1,1,1,1),each=30),ncol=nct))
# colnames(ctype_tru)<-paste0("type",seq_len(nct))
# ctype_wilcoxtru<-rbind(matrix(rep(c(1,0,1,1,0,1),each=20),ncol=choose(nct,2)),matrix(rep(c(0,0,0,0,0,0),each=30),ncol=choose(nct,2)),matrix(rep(c(1,0,1,1,0,1),each=20),ncol=choose(nct,2)),matrix(rep(c(0,0,0,0,0,0),each=30),ncol=choose(nct,2)))
# colnames(ctype_wilcoxtru)<-paste0("pairwise-type",seq_len(choose(nct,2)))

#######################################
## Hierarchical clustering ############
###############################################
# gene_dist<-dist(ratio33,method = "manhattan")
# gene_hclust <- hclust(gene_dist, method = "ward.D2")
# # plot(gene_hclust, labels = FALSE)
# # abline(h = 6, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# # cut genes into 5 group 
# gene_cluster<-cutree(gene_hclust, k = 5) %>% 
#   # turn the named vector into a tibble
#   enframe(name="gene",value = "cluster")
# gene_cluster[["true"]]<-matrix(c(rep(1,10),rep(2,10),rep(3,30),rep(4,10),rep(5,10),rep(3,30)),ncol = 1)

#######################################
## kmeans ############
###############################################
set.seed(123)
####select optimal k##########################
# fviz_nbclust(ratio33, kmeans, method = "wss")
km.res <- kmeans(ratio33, 5,nstart = 15)
gene_cluster<-data.frame(gene=seq_len(ngene),cluster=km.res$cluster,true=matrix(c(rep(1,10),rep(2,10),rep(3,30),rep(4,10),rep(5,10),rep(3,30)),ncol = 1))

ARI_gene[t]<-adjustedRandIndex(gene_cluster$cluster,gene_cluster$true)
print(ARI_gene[t])

ct_bin<-matrix(0,ncol=nct,nrow = ngene)
ct_lin<-matrix(0,ncol=nct,nrow = ngene)
ct_wilcoxon<-matrix(0,ncol = choose(nct,2),nrow = ngene)

########################################################################
## detect cell type according allelic ratio ############################
#####################################################
for (i in 1:5) {
print(paste("gene pattern",i))
poi<-which(gene_cluster$cluster==i)
gene1_mean<-colMeans(ratio3[poi,],na.rm = T)
count1_mean<-colMeans(cts2[poi,])
truth[(5*t-5+i),]<-sapply(1:nct, function(j){mean(gene1_mean[(n*j-n+1):(n*j)],na.rm=T)})

rent<-tibble(ratio=gene1_mean,X=as.factor(rep(c(1,2,3,4),each=n)),cts=count1_mean)
levels(rent$X) <- c("type1", "type2", "type3","type4")
group_by(rent, X) %>%
  summarise(
    count = n(),
    median = median(ratio, na.rm = TRUE),
    IQR = IQR(ratio, na.rm = TRUE)
  )
## Gaussian likelihood #####################
formu <- ratio ~ p(X, pen = "gflasso", refcat = "type1")
fit2 <- glmsmurf(formula = formu, family=gaussian(), data = rent,
                                    pen.weights = "glm", lambda = "cv1se.mse",
                                    control = list(lambda.length = 50L,k=10),pen.weights.return=T) #pen.weights="glm" will fused coefficients more
                                                                               # and more consistent with truth sometimes 
print(coef(fit2))
coeflin<-coef(fit2)
coeflin[2:4]<-coeflin[2:4]+coeflin[1]
type_lin=match(coeflin, unique(coeflin))
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
fit <- glmsmurf(formula = formu, family=binomial(link = "logit"), data = rent, weights = count1_mean,
                       pen.weights = "glm", lambda = "cv1se.mse",
                       control = list(lambda.length = 50L,k=10),pen.weights.return=T) #
print(coef(fit))
coefbin<-coef(fit)
coefbin[2:4]<-coefbin[2:4]+coefbin[1]
type=match(coefbin, unique(coefbin))
ct_bin[poi,]<-rep(type,each=length(poi))
# plot(fit)
# plot_reest(fit)
# summary(fit)

bin[(5*t-5+i),]=sapply(1:nct, function(j){exp(coefbin[j])/(1 + exp(coefbin[j]))})
lin[(5*t-5+i),]=sapply(1:nct, function(j){coeflin[j]})
}


res2 <- pbsapply(1:ngene, function(i) {adjustedRandIndex(ctype_tru[i,],ct_lin[i,])})
ARI_lin[t]<-mean(res2) #0.6514286
reswilcox <- pbsapply(1:ngene, function(i) {adjustedRandIndex(ctype_wilcoxtru[i,], ct_wilcoxon[i,])})
ARI_wilcoxon[t]<-mean(reswilcox) #0.6118919
res <- pbsapply(1:ngene, function(i) {adjustedRandIndex(ctype_tru[i,],ct_bin[i,])})
ARI_bin[t]<-mean(res) #0.72
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
                lindiff4=lin[,4]-truth[,4])
boxplot(box,ylab="estimated-truth")
abline(h=0, col=2)

ARI<-data.frame(Gaussian=ARI_lin,Wilcoxon=ARI_wilcoxon,Binomial=ARI_bin)
boxplot(ARI,ylab="Adjusted rand index")

boxplot(ARI_gene,ylab="Adjusted rand index")
