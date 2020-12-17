cmd_args=commandArgs(TRUE)

n <- as.numeric(cmd_args[1]) # cells per cell type
cnt <- as.numeric(cmd_args[2]) # the mean "high" total count
k <-as.numeric(cmd_args[3])
out <- cmd_args[4]

k<-40
n <- 400 # cells per cell type
cnt <- 100# the mean "high" total count
source("/proj/milovelab/mu/SC-ASE/simulation/cluster.R")
source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
library(smurf)
library(pbapply)
library(mclust)
library(emdbook)

ans <- pbsapply(1:200, function(i) {
# ans <- for(i in 1:10) {
  set.seed(i)
  k <- k # number of cell types
  low_count <- 5 # the mean "low" total count
  mean_total_count <- rep(rep(c(low_count, cnt),each=n/2), times=k) # total count
  size <- rpois(n * k, mean_total_count)
  size[size == 0] <- 1
  p.vec <- (3 + rep(seq(from=-2,to=2,length.out=k/2),each=2))/6
  # p.vec <- (5 + rep(seq(from=3.5,to=4.5,length.out=k/2),each=2))/10
  p <- rep(p.vec, each=n) # true prob
  #y <- rbinom(k*n, prob=p, size=size) # obs counts
  y <- rbetabinom(k*n, prob=p, size=size, theta=2) # obs counts
  r <- y/size # ratio
  x <- factor(rep(1:k,each=n)) # cell type dummy
  f <- r ~ p(x, pen="gflasso", refcat="1") # formula
  dat=data.frame(x=x,ratio=r,cts=size)
  # t <- system.time(fit<-fusedlasso(formula=f,model="binomial",data=dat,ncores=1))[[3]] # binomial
  # t2 <- system.time(fit2<-fusedlasso(formula=f,model="gaussian",dat,ncores=1,lambda.length=50L))[[3]] # gaussian
  # co <- coef(fit)
  # co <- co + c(0,rep(co[1],k-1))
  # a <- adjustedRandIndex(factor(p.vec), factor(co))
  # co <- coef(fit2)
  # co <- co + c(0,rep(co[1],k-1))
  # a2 <- adjustedRandIndex(factor(p.vec), factor(co))
  t3<-system.time(fit3<-wilcox(dat,nct=k,method = "BH",threshold=0.2))[[3]]
  a3<-adjustedRandIndex(factor(p.vec), factor(fit3))
  fit4<-wilcox_adj(dat,nct=k,k=5,lambda = 0.05,threshold = c(0.01,0.05,0.1,0.2,0.3))
  fit5<-wilcox_adj(dat,nct=k,k=5,lambda = 0.1,threshold = c(0.01,0.05,0.1,0.2,0.3))
  fit6<-wilcox_adj(dat,nct=k,k=5,lambda = 0.5,threshold = c(0.01,0.05,0.1,0.2,0.3))
  fit7<-wilcox_adj(dat,nct=k,k=5,lambda = 1,threshold = c(0.01,0.05,0.1,0.2,0.3))
  a4<-adjustedRandIndex(factor(p.vec), factor(fit4))
  a5<-adjustedRandIndex(factor(p.vec), factor(fit5))
  a6<-adjustedRandIndex(factor(p.vec), factor(fit6))
  a7<-adjustedRandIndex(factor(p.vec), factor(fit7))
  out <- c(a3,a4,a5,a6,a7)
  out
}, cl=10)
# }

# save the results as a data.frame
dat <- data.frame(type=rep(c("bin","gau","Wilcoxon"),each=ncol(ans)),
                  ARI=as.vector(t(ans[1:3,])),
                  time=as.vector(t(ans[4:6,])),
                  n=n,
                  cnt=cnt,
                  k=k)

# write out as a table
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
