fusedlasso<-function(formula,model="gaussian",data,lambda="cv1se.dev",k=5,...){
  misspoi<-which(!is.na(data$ratio))
  if(model=="binomial"){
      # need to use tryCatch to avoid lambda.max errors
      try1 <- tryCatch({
        fit <- glmsmurf(formula=formula, family=binomial(link="logit"), data=data,
                        weights=data$cts[misspoi], pen.weights="glm.stand", lambda=lambda, 
                        control=list(lambda.length=25L, k=k,...));
        TRUE
      }, error=function(e) {
        message("Failed determining the maximum of lambda, run gaussian model instead")
        fit <-glmsmurf(formula=formula, family=gaussian(), data=data,
                pen.weights="glm.stand", lambda=lambda,control=list(k=k,...))});
    }
  if(model=="gaussian"){
    # need to use tryCatch to avoid lambda.max errors
    try1 <- tryCatch({
       fit <- glmsmurf(formula=formula, family=gaussian(), data=data,
                        pen.weights="glm.stand", lambda=lambda, 
                        control=list(k=k,...));
       TRUE
    }, error=function(e) {
      message("Failed determining the maximum of lambda, using standardized adaptive GAM weight instead")
      fit <-glmsmurf(formula=formula, family=gaussian(), data=data,
                     pen.weights="gam.stand", lambda=lambda,control=list(k=k,...))});
  }
  return(fit)
}

boot_fusedlasso<-function(formula,data,indices,model, lambda1,lambda2,...){
  data_b<-data[indices,]
  poi<-which(!is.na(data_b$ratio))
  if(model=="binomial"){
      fit1 <- glmsmurf(formula=formula, family=binomial(link="logit"), data=data_b,
                      weights=data_b$cts[poi], pen.weights="glm.stand", lambda=lambda1,...);
      co <- coef_reest(fit1)
      co <- co + c(0,rep(co[1],nct-1))
      co <-1/(1+exp(-co))
  }
  if(model=="gaussian"){
    fit2 <- glmsmurf(formula=formula, family=gaussian(), data=data_b,
                    pen.weights="glm.stand", lambda=lambda2,...)
    co <- coef_reest(fit2)
    co <- co + c(0,rep(co[1],nct-1))
  }
  return(co)
}

## Wilcoxon Rank Sum test ###################################
wilcox<-function(data,nct,method="holm",threshold=0.05,...){
res <- pairwise.wilcox.test(data$ratio,data$x, p.adjust.method =method,...)
# res <- pairwise.wilcox.test(dat$ratio,dat$x, p.adjust.method ="holm")
adj<-as.data.frame(res$p.value)[lower.tri(res$p.value, diag = T)]
b<- matrix(0, nct, nct)
b[lower.tri(b, diag=FALSE)]=adj
b2<-b+t(b)
diag(b2)<-1
bb<-ifelse(b2<threshold,1,0) #Binarize p-value to be seen as dismilarity matrix
# bb<-ifelse(b2<0.7,1,0) #Binarize p-value to be seen as dismilarity matrix
clust<-hclust(as.dist(bb))
# plot(clust)
my.clusters<-cutree(clust,h=0)
return(my.clusters)
}

library(caret)
library(tidyverse)
wilcox_adj<-function(data,nct,k,threshold,lambda,method="BH",...){
  set.seed(123)
  foldInds <- createFolds(data$x, k=k, list=TRUE, returnTrain=T)
  out<-list()
  a<-rep(0,k)
  obj<-sapply (1:length(threshold), function(j){
    for(i in 1:5){
      fit3<-wilcox(data[foldInds[[i]],],nct=nct,method = method,threshold=threshold[j])
      label<-tibble(type=factor(seq_along(1:nct)),par=fit3)
      test<-dat[-foldInds[[i]],]
      test2<-test %>% left_join(label,by=c("x"="type"))
      test2<-test2 %>% group_by(par) %>% mutate(grpmean=mean(ratio))
      a[i]<-sum(pbsapply(1:nrow(test2), function(i){abs(test2$ratio[i]-test2$grpmean[i])}))
    }
    fit<-wilcox(data,nct=nct,method = method,threshold=threshold[j])
    b<-sum(a)+lambda*length(unique(fit))*log(nrow(data))
    out[["cl"]]<-fit
    out[["loss"]]<-b
    return(out)
  }) 
  cl<-as_tibble(do.call(rbind, obj[seq(1,length(obj), by = 2)]))
  loss<-as_tibble(do.call(rbind, obj[seq(2,length(obj), by = 2)]))
  return(cl[which.min(unlist(loss)),])
}


