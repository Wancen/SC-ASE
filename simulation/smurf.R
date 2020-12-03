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

#n<-10
#cnt<-50
ngene<-60
nct<-8
x <- factor(rep(1:nct,each=n))
mu1 <- 5
nb.disp <- 1/100
ncl<-3 #3 gene cluster. AI & no AI
step1<-2 #First AI step [0-4]
step2<-1 #second AI step [0-4]

ans2 <- pbsapply(1:50, function(i) {
  set.seed(i)
  
  # total count
  cts <- matrix(rep(c(rnbinom(ngene*n/2,mu=mu1,size=1/nb.disp),
                     rnbinom(ngene*n/2,mu=cnt,size=1/nb.disp)),nct),ncol = nct*n)
  colnames(cts)<-paste0("cell",1:(nct*n))
  
  p.vec <- (4 + rep(c(seq(from=-step1,to=step1,length.out=nct/2),rep(0,nct/2),seq(from=-step2,to=step2,length.out=nct/2)),each=2))/8 
  p <- rep(p.vec, each=n*nct*ngene/length(p.vec)) # true prob
  nclgene<-ngene/ncl #number genes within cluster
  nclcell<-nct*n*nclgene #number elements within cluster
  ase.cts<-pblapply(1:ncl,function(m) {
    matrix(rbetabinom(nclcell, prob=p[(nclcell*m-nclcell+1):(nclcell*m)], size=cts[(m*nclgene-nclgene+1):(m*nclgene),], theta=10),ncol = nct*n)})
  ase.cts<-do.call(rbind,ase.cts)
  ratio<-(ase.cts)/(cts)
  ratio_pseudo<-(ase.cts+1)/(cts+2) ## pseudo allelic ratio for gene clustering 
  level<-paste0(rep("type",nct),1:nct) # pheatmap of ratio
  anno_df <- data.frame(celltype=rep(level,each=n), row.names=colnames(ratio_pseudo))
  pheatmap(ratio_pseudo, cluster_rows = FALSE, cluster_cols = FALSE,annotation_col=anno_df,show_colnames = F,
           color = colorRampPalette(colors = c("blue","white","red"))(100))

  cluster<-genecluster(ratio_pseudo,nct=nct,G=3) #return gene cluster 
  acl<-ARI(cluster,rep(c(0,1,2),each=ngene/ncl))
  acl
  # modeling
  out<-list()
  for (j in 1:ncl) {

    poi<-which(cluster==unique(factor(cluster))[j]) # gene position
    #poi<-(20*j-20+1):(20*j) # gene position
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
  }, cl=4)

ans2 <- do.call(cbind, ans2)
# a<-pbsapply(1:ncol(ans),function(i){adjustedRandIndex(ans[,i],rep(c(0,1),each=ngene/ncl))})
# save the results as a data.frame
dat <- data.frame(type=rep(c("bin","gau"),each=ncol(ans)),
                  ARI_cl=as.vector(t(ans2[1,])),
                  ARI=as.vector(t(ans2[2:3,])),
                  cl=rep(c("largeAI","NAI","smallAI"),n=ncol(ans)/ncl*2),
                  time=as.vector(t(ans2[4:5,])),
                  n=n,
                  cnt=cnt)

# write out as a table
# write.table(dat, file="out.csv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")








