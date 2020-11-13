cmd_args=commandArgs(TRUE)

n <- as.numeric(cmd_args[1]) # cells per cell type
cnt <- as.numeric(cmd_args[2]) # the mean "high" total count
out <- cmd_args[3]

library(smurf)
library(pbapply)
library(mclust)
library(emdbook)

ans <- pbsapply(1:200, function(i) {

  set.seed(i)
  k <- 10 # number of cell types
  low_count <- 5 # the mean "low" total count
  mean_total_count <- rep(rep(c(low_count, cnt),each=n/2), times=k) # total count
  size <- rpois(n * k, mean_total_count)
  size[size == 0] <- 1
  p.vec <- (3 + rep(seq(from=-2,to=2,length.out=k/2),each=2))/6
  p <- rep(p.vec, each=n) # true prob
  #y <- rbinom(k*n, prob=p, size=size) # obs counts
  y <- rbetabinom(k*n, prob=p, size=size, theta=10) # obs counts
  r <- y/size # ratio
  x <- factor(rep(1:k,each=n)) # cell type dummy
  f <- r ~ p(x, pen="gflasso", refcat="1") # formula

  # binomial
  t <- system.time({

    # need to use tryCatch to avoid lambda.max errors
    try1 <- tryCatch({
      fit <- glmsmurf(formula=f, family=binomial(link="logit"), data=data.frame(x,r),
                      weights=size, pen.weights="glm.stand", lambda="cv1se.dev", 
                      control=list(lambda.length=20L, k=5, ncores=1));
      TRUE
      }, error=function(e) FALSE)
    
  })[[3]] # saving the elapsed time
  # gaussian
  t2 <- system.time({
    
    try2 <- tryCatch({
      fit2 <- glmsmurf(formula=f, family=gaussian(), data=data.frame(x,r),
                       pen.weights="glm.stand", lambda="cv1se.dev", 
                       control=list(lambda.length=20L, k=5, ncores=1));
      TRUE
      }, error=function(e) FALSE)
    
  })[[3]]
  
  if (try1 & try2) {
    co <- coef(fit)
    co <- co + c(0,rep(co[1],k-1))
    a <- adjustedRandIndex(factor(p.vec), factor(co))
    co <- coef(fit2)
    co <- co + c(0,rep(co[1],k-1))
    a2 <- adjustedRandIndex(factor(p.vec), factor(co))
    out <- c(a,a2,t,t2)
  } else {
    out <- NULL
  }
  out
}, cl=6)

# dealing with the tryCatch errors...
if (!is.matrix(ans)) {
  ans <- do.call(cbind, ans)
}

# save the results as a data.frame
dat <- data.frame(type=rep(c("bin","gau"),each=ncol(ans)),
                  ARI=as.vector(t(ans[1:2,])),
                  time=as.vector(t(ans[3:4,])),
                  n=n,
                  cnt=cnt)

# write out as a table
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
