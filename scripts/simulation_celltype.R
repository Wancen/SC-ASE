library(emdbook)
library(tibble)
n <- 80 # cells per cluster
ngene <- 100
f <- factor(rep(1:4,each=n/4)) # two clusters
mu1 <- 5
mu2 <- 100
nb.disp <- 1/100
set.seed(1)

########################################################################
# total counts are not same ################################################
cts_t3 <- matrix(
  rnbinom(ngene*n*9/10,
          mu=mu1,
          size=1/nb.disp),
  ncol=n*9/10)

cts_f1 <- matrix(
  rnbinom(ngene*n/10,
          mu=mu2,
          size=1/nb.disp),
  ncol=n/10)
cts2<-cbind(cts_t3,cts_f1)
theta <- runif(ngene,1,50)
colnames(cts2) <-paste0("cell",seq_len(ncol(cts2)))  

#### Different gene with same allelic ratio #######################
# create complex simulations for first two cluster here:
prob1 <- rnorm(ngene,.5,.05) # close to 0.5
prob1[1:10] <- rnorm(10,.75,.05)
prob1[11:20] <- rnorm(10,.25,.05)

# counts for allele 2
ase.cts_t2 <- matrix(
  rbetabinom(prod(dim(cts2))/2, 
             prob=prob1,
             size=cts2[,1:(n/2)], 
             theta=rep(theta,ncol(cts2)/2)),
  nrow=ngene)
colnames(ase.cts_t2) <- paste0("cell",seq_len(ncol(ase.cts_t2)))

# create complex simulations for next two cluster here:
prob2 <- rnorm(ngene,.5,.05) # close to 0.5
prob2[51:60] <- rnorm(10,.75,.05)
prob2[61:70] <- rnorm(10,.25,.05)

# counts for allele 2
ase.cts_f2 <- matrix(
  rbetabinom(prod(dim(cts2))/2, 
             prob=prob2,
             size=cts2[,(n/2+1):n], 
             theta=rep(theta,ncol(cts2)/2)),
  nrow=ngene)
colnames(ase.cts_f2) <- seq_len(ncol(ase.cts_f2))

ase.cts3<-cbind(ase.cts_t2,ase.cts_f2)
colnames(ase.cts3) <-paste0("cell",seq_len(ncol(cts2)))  
ratio3<-ase.cts3/cts2
library(pheatmap)
anno_df <- data.frame(f, row.names=colnames(cts2))
# graphics.off()
# jpeg("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/fusion idea/AI but not DAI.jpg")
pheatmap(ase.cts3, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col=anno_df,labels_col = seq_len(80))
pheatmap(ratio3, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col=anno_df,labels_col = seq_len(80))
# dev.off()
#### same gene with different allelic ratio #######################
# create complex simulations for first two cluster here:
prob1 <- rnorm(ngene,.5,.05) # close to 0.5
prob1[1:10] <- rnorm(10,.7,.05)
prob1[11:20] <- rnorm(10,.3,.05)

# counts for allele 2
ase.cts_t2 <- matrix(
  rbetabinom(prod(dim(cts2))/2, 
             prob=prob1,
             size=cts2[,1:(n/2)], 
             theta=rep(theta,ncol(cts2)/2)),
  nrow=ngene)
colnames(ase.cts_t2) <- seq_len(ncol(ase.cts_t2))

# create complex simulations for next two cluster here:
prob2 <- rnorm(ngene,.5,.05) # close to 0.5
prob2[51:60] <- rnorm(10,.4,.04)
prob2[61:70] <- rnorm(10,.6,.04)

# counts for allele 2
ase.cts_f2 <- matrix(
  rbetabinom(prod(dim(cts2))/2, 
             prob=prob2,
             size=cts2[,(n/2+1):n], 
             theta=rep(theta,ncol(cts2)/2)),
  nrow=ngene)
colnames(ase.cts_f2) <- seq_len(ncol(ase.cts_f2))

ase.cts4<-cbind(ase.cts_t2,ase.cts_f2)
colnames(ase.cts4) <-paste0("cell",seq_len(ncol(cts2)))  
ratio4<-ase.cts4/cts2
library(pheatmap)
anno_df <- data.frame(f, row.names=colnames(cts2))
# graphics.off()
# jpeg("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/fusion idea/diff ratio same gene.jpg")
pheatmap(ase.cts4, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col=anno_df,labels_col = seq_len(80))
pheatmap(ratio4, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col=anno_df,labels_col = seq_len(80))
pheatmap(cts2, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col=anno_df,labels_col = seq_len(80))
# dev.off()


