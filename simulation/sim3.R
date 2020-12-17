cmd_args=commandArgs(TRUE)

n <- as.numeric(cmd_args[1]) # cells per cell type
cnt <- as.numeric(cmd_args[2]) # the mean "high" total count
k <-as.numeric(cmd_args[3])
out <- cmd_args[4]


n <- 200 # cells per cell type
cnt <- 100# the mean "high" total count
source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
library(smurf)
library(pbapply)
library(mclust)
library(emdbook)

ans <- pbsapply(1:10, function(i) {
  set.seed(i)
  k <- k # number of cell types
  low_count <- 5 # the mean "low" total count
  mean_total_count <- rep(rep(c(low_count, cnt),each=n/2), times=k) # total count
  size <- rpois(n * k, mean_total_count)
  size[size == 0] <- 1
  p.vec <- (3 + rep(seq(from=-2,to=2,length.out=5),each=k/5))/6
  # p.vec <- (5 + rep(seq(from=3.5,to=4.5,length.out=k/2),each=2))/10
  p <- rep(p.vec, each=n) # true prob
  #y <- rbinom(k*n, prob=p, size=size) # obs counts
  y <- rbetabinom(k*n, prob=p, size=size, theta=2) # obs counts
  r <- y/size # ratio
  x <- factor(rep(1:k,each=n)) # cell type dummy
  f <- r ~ p(x, pen="gflasso", refcat="1") # formula
  dat=data.frame(x=x,ratio=r,cts=size)
  t <- system.time(fit<-fusedlasso(formula=f,model="binomial",data=dat,ncores=1))[[3]] # binomial
  t2 <- system.time(fit2<-fusedlasso(formula=f,model="gaussian",dat,ncores=1,lambda.length=50L))[[3]] # gaussian
  co <- coef(fit)
  co <- co + c(0,rep(co[1],k-1))
  a <- adjustedRandIndex(factor(p.vec), factor(co))
  co <- coef(fit2)
  co <- co + c(0,rep(co[1],k-1))
  a2 <- adjustedRandIndex(factor(p.vec), factor(co))
  t3<-system.time(fit3<-wilcox(dat,nct=k))[[3]]
  a3<-adjustedRandIndex(factor(p.vec), factor(fit3))
  out <- c(a,a2,a3,t,t2,t3)
  out
},cl=4)
# }


# save the results as a data.frame
dat <- data.frame(type=rep(c("bin","gau","Wilcoxon"),each=ncol(ans)),
                  ARI=as.vector(t(ans[1:3,])),
                  time=as.vector(t(ans[4:6,])),
                  n=n,
                  cnt=cnt)

# write out as a table
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
