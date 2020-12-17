# cmd_args=commandArgs(TRUE)
# 
# ngenecl <- as.numeric(cmd_args[1]) # cells per cell type
# out <- cmd_args[2]

source("/proj/milovelab/mu/SC-ASE/simulation/cluster.R")
source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
library("smurf")
library(emdbook)
library(mclust)
library(pbapply)
library(aricode)
library(pheatmap)


ngenecl<-80
n<-10
cnt<-50
ncl<-4 #4 gene cluster. large AI,NO AI, consistent AI, small AI
ngene<-ncl*ngenecl
nct<-8
x <- factor(rep(1:nct,each=n))
mu1 <- 5
nb.disp <- 1/100
ncl<-4 #4 gene cluster. large AI,NO AI, consistent AI, small AI
step1<-4 #First AI step [0-4]

ans <- pbsapply(1:10, function(i) {
  set.seed(i)
  
  # total count
  cts <- matrix(rep(c(rnbinom(ngene*n/2,mu=mu1,size=1/nb.disp),
                     rnbinom(ngene*n/2,mu=cnt,size=1/nb.disp)),nct),ncol = nct*n)
  colnames(cts)<-paste0("cell",1:(nct*n))
  
  p.vec <- (5 + rep(c(seq(from=-step1,to=step1,length.out=nct/2),rep(0,nct/2),rep(2,nct/2),seq(from=3.5,to=4.5,length.out=nct/2)),each=2))/10 
  p <- rep(p.vec, each=n*nct*ngene/length(p.vec)) # true prob
  nclgene<-ngene/ncl #number genes within cluster
  nclcell<-nct*n*nclgene #number elements within cluster
  ase.cts<-lapply(1:ncl,function(m) {
    matrix(rbetabinom(nclcell, prob=p[(nclcell*m-nclcell+1):(nclcell*m)], size=cts[(m*nclgene-nclgene+1):(m*nclgene),], theta=10),ncol = nct*n)})
  ase.cts<-do.call(rbind,ase.cts)
  ratio<-(ase.cts)/(cts)
  ratio_pseudo<-(ase.cts+1)/(cts+2) ## pseudo allelic ratio for gene clustering 
  level<-paste0(rep("type",nct),1:nct) # pheatmap of ratio
  # anno_df <- data.frame(celltype=rep(level,each=n), row.names=colnames(ratio_pseudo))
  # pheatmap(ratio_pseudo, cluster_rows = FALSE, cluster_cols = FALSE,annotation_col=anno_df,show_colnames = F,
  #          color = colorRampPalette(colors = c("blue","white","red"))(100))
  
  cluster<-genecluster(ratio_pseudo,nct=nct,G=ncl) #return gene cluster 
  mcl<-adjustedRandIndex(cluster,rep(1:ncl,each=ngene/ncl))
  # modeling
  out<-list()
  for (j in 1:ncl) {
    # poi<-which(cluster==unique(factor(cluster))[j]) # gene position
    poi<-(ngenecl*j-ngenecl+1):(ngenecl*j) # gene position
    r<-as.vector(ratio[poi,])
    size<-as.vector(cts[poi,])
    data=data.frame(x=rep(x,each=length(poi)),ratio=r,cts=size)
    f <- ratio ~ p(x, pen="gflasso", refcat="1") # formula
    t <- system.time(fit<-fusedlasso(formula=f,model="binomial",data=data,ncores=1))[[3]] # saving the elapsed time
    t2 <- system.time(fit2<-fusedlasso(formula=f,model="gaussian",data,ncores=1))[[3]] # saving the elapsed time
    co <- coef(fit)
    co <- co + c(0,rep(co[1],nct-1))
    a <- adjustedRandIndex(factor(p.vec[(nct*j-nct+1):(nct*j)]), factor(co))
    co2 <- coef(fit2)
    co2 <- co2 + c(0,rep(co2[1],nct-1))
    a2 <- adjustedRandIndex(factor(p.vec[(nct*j-nct+1):(nct*j)]), factor(co2))
    t3 <- system.time(fit3<-wilcox(data,nct,method="holm"))[[3]]
    a3<-adjustedRandIndex(factor(p.vec[(nct*j-nct+1):(nct*j)]), factor(fit3))
    out[[j]]=c(a,a2,a3,t,t2,t3)
   }
  out
},cl=10)

ans <- do.call(cbind, ans)
ans2<-rbind(ans,rep(1:200,each=4))
id<-which(!is.numeric(ans2[3,]))
# a<-pbsapply(1:ncol(ans),function(i){adjustedRandIndex(ans[,i],rep(c(0,1),each=ngene/ncl))})
# save the results as a data.frame
dat2 <- data.frame(type=rep(c("bin","gau","wilcoxon"),each=ncol(ans)),
                  ARI_mcl=as.vector(t(ans[1,])),
                  ARI=as.vector(t(ans[2:4,])),
                  cl=rep(c("large AI gap","NAI","consisAI","small AI gap"),n=ncol(ans)/ncl*3),
                  time=as.vector(t(ans[5:7,])))

# write out as a table
# write.table(dat, file="/proj/milovelab/mu/SC-ASE/simulation/csv/sim2.csv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
write.table(dat2, file="/proj/milovelab/mu/SC-ASE/simulation/csv/sim2_80.csv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")








