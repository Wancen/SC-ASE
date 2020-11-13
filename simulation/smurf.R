cmd_args=commandArgs(TRUE)

n <- as.numeric(cmd_args[1]) # cells per cell type
cnt <- as.numeric(cmd_args[2]) # the mean "high" total count
out <- cmd_args[3]

source("/proj/milovelab/mu/SC-ASE/simulation/cluster.R")
source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
library("smurf")
library(emdbook)
library(mclust)
library(pbapply)
library(aricode)
library(pheatmap)

ngene<-100
nct<-8
x <- factor(rep(1:nct,each=n))
mu1 <- 5
nb.disp <- 1/100
ncl<-2 #2 gene cluster. AI & no AI

ans <- pbsapply(1:50, function(i) {
  set.seed(i)
  
  # total count
  cts <- matrix(rep(c(rnbinom(ngene*n/2,mu=mu1,size=1/nb.disp),
                     rnbinom(ngene*n/2,mu=cnt,size=1/nb.disp)),nct),ncol = nct*n)
  colnames(cts)<-paste0("cell",1:(nct*n))
  
  p.vec <- (3 + rep(c(seq(from=-2,to=2,length.out=nct/2),rep(0,nct/2)),each=2))/6 
  p <- rep(p.vec, each=n*nct*ngene/length(p.vec)) # true prob
  nclcell<-nct*n*ngene/ncl # #cells within cluster
  ase.cts <- rbind(matrix(rbetabinom(nclcell, prob=head(p,nclcell), size=cts[1:(ngene/ncl),], theta=10),ncol = nct*n),
                 matrix(rbetabinom(nclcell, prob=tail(p,nclcell), size=cts[(ngene/ncl+1):100,], theta=10),ncol = nct*n))# obs counts

  ratio<-(ase.cts)/(cts)
  ratio_pseudo<-(ase.cts+1)/(cts+2) ## pseudo allelic ratio for gene clustering 
  # level<-paste0(rep("type",nct),1:nct) # pheatmap of ratio
  # anno_df <- data.frame(celltype=rep(level,each=n), row.names=colnames(ratio_pseudo))
  # pheatmap(ratio_pseudo, cluster_rows = FALSE, cluster_cols = FALSE,annotation_col=anno_df,show_colnames = F)

  cluster<-genecluster(ratio_pseudo,nct=nct,G=ncl) #return gene cluster 
  acl<-ARI(cluster,rep(c(0,1),each=ngene/ncl))
  
  ## modeling
  out<-list()
  for (j in 1:ncl) {

    poi<-which(cluster==unique(factor(cluster))[j]) # gene position
    r<-as.vector(ratio[poi,])
    poi2<-which(!is.na(r)) # remove missing value
    size<-as.vector(cts[poi,])
    data=data.frame(x=rep(x,each=length(poi)),r)
    f <- r ~ p(x, pen="gflasso", refcat="1") # formula
    t <- system.time(fit<-fusedlasso(formula=f,model="binomial",data,size,misspoi=poi2))[[3]] # saving the elapsed time
    t2 <- system.time(fit2<-fusedlasso(formula=f,model="gaussian",data,size))[[3]] # saving the elapsed time
    co <- coef(fit)
    co <- co + c(0,rep(co[1],nct-1))
    a <- adjustedRandIndex(factor(p.vec[(nct*j-nct+1):(nct*j)]), factor(co))
    co <- coef(fit2)
    co <- co + c(0,rep(co[1],nct-1))
    a2 <- adjustedRandIndex(factor(p.vec[(nct*j-nct+1):(nct*j)]), factor(co))
    out[[j]]=c(acl,a,a2,t,t2)
    }
    out
  }, cl=6)

ans <- do.call(cbind, ans)

# save the results as a data.frame
dat <- data.frame(type=rep(c("bin","gau"),each=ncol(ans)/2),
                  ARI_cl=as.vector(t(ans[1,])),
                  ARI_AI=as.vector(t(ans[2:3,c(TRUE,FALSE)])),
                  ARI_NAI=as.vector(t(ans[2:3,c(FALSE,TRUE)])),
                  time1=as.vector(t(ans[4:5,c(TRUE,FALSE)])),
                  time2=as.vector(t(ans[4:5,c(FALSE,TRUE)])),
                  n=n,
                  cnt=cnt)

# write out as a table
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")








