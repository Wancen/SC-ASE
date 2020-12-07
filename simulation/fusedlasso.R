fusedlasso<-function(formula,model="gaussian",data,lambda="cv1se.dev",k=5,...){
  misspoi<-which(!is.na(data$ratio))
  if(model=="binomial"){
      # need to use tryCatch to avoid lambda.max errors
      try1 <- tryCatch({
        fit <- glmsmurf(formula=formula, family=binomial(link="logit"), data=data,
                        weights=data$cts[misspoi], pen.weights="glm.stand", lambda=lambda, 
                        control=list(lambda.length=20L, k=k,...));
        TRUE
      }, error=function(e) {
        message("Failed determining the maximum of lambda, run gaussian model instead")
        fit <-glmsmurf(formula=formula, family=gaussian(), data=data,
                pen.weights="glm.stand", lambda=lambda,control=list(lambda.length=20L, k=k,...))});
    }
  if(model=="gaussian"){
       fit <- glmsmurf(formula=formula, family=gaussian(), data=data,
                        pen.weights="glm.stand", lambda=lambda, 
                        control=list(lambda.length=20L, k=k,...))
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
wilcox<-function(data,nct){
res <- pairwise.wilcox.test(data$ratio,data$x, p.adjust.method ="BH")
adj<-as.data.frame(res$p.value)[lower.tri(res$p.value, diag = T)]
b<- matrix(0, nct, nct)
b[lower.tri(b, diag=FALSE)]=adj
b2<-b+t(b)
diag(b2)<-1
bb<-ifelse(b2<0.05,1,0) #Binarize p-value to be seen as dismilarity matrix
clust<-hclust(as.dist(bb))
# plot(clust)
my.clusters<-cutree(clust,h=0)
return(my.clusters)
}