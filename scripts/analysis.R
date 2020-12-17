load("/proj/milovelab/mu/data/postprocess.rda")
source("/proj/milovelab/mu/SC-ASE/simulation/cluster.R")
source("/proj/milovelab/mu/SC-ASE/simulation/fusedlasso.R")
library(tidyverse)
library(pheatmap)
library(doParallel)
library(boot)
library(pbapply)
nct=10
ncores=4
cluster<-genecluster(ratio_psedo,nct,G=c(4,8,12,16,20,24))
table(cluster)

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     ...)
train.index <- createDataPartition(dat$x, p = .8, list = FALSE)
train <- dat[ train.index,]
test  <- dat[-train.index,]

foldInds <- createFolds(dat$x, k=5, list=TRUE, returnTrain=T)
lapply(foldInds, function(ii) table(dat$x[ii])) ## verify stratification

a<-rep(0,5)
for(i in 1:5){
  t3<-system.time(fit3<-wilcox(dat[foldInds[[i]],],nct=k,method = "BH",threshold=0.05))[[3]]
  label<-tibble(type=factor(seq_along(1:k)),par=fit3)
  test<-dat[-foldInds[[i]],]
  test2<-test %>% left_join(label,by=c("x"="type"))
  test2<-test2 %>% group_by(par) %>% mutate(grpmean=mean(ratio))
  a[i]<-sum(pbsapply(1:nrow(test2), function(i){(test2$ratio[i]-test2$grpmean[i])^2}))
}
sum(a)


a3<-adjustedRandIndex(factor(p.vec), factor(fit3))
t4<-system.time(fit3<-wilcox(dat,nct=k,method = "BH",threshold=0.2))[[3]]
a4<-adjustedRandIndex(factor(p.vec), factor(fit3))

ans<-as_tibble(ans)
row.names(ans)<-c("a","a1","a2","a3","t1","t2","t3")
ans2<-cbind(t(ans),id=rep(1:100,each=4))

fit <- glmsmurf(formula=f, family=gaussian(), data=dat,
                pen.weights="glm.stand", lambda="cv1se.dss", 
                control=list(lambda.length=20L, k=5))
fit <- glmsmurf(formula=f, family=binomial(link="logit"), data=dat,
                weights=dat$cts[which(!is.na(dat$ratio))], pen.weights="glm.stand", lambda="cv1se.dev", 
                control=list(lambda.length=25L, k=5))
dat %>% group_by(x) %>% summarise(wm=weighted.mean(ratio,cts),m=mean(ratio,na.rm = T)) %>% data.frame()
co2<-fit$coefficients
co2 <- co2 + c(0,rep(co2[1],k-1))
co2

read.table("/proj/milovelab/mu/SC-ASE/simulation/csv/",sep = ",")
write.table(dat, file=out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")