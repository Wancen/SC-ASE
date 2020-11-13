fusedlasso<-function(formula,model="gaussian",data,size,misspoi=NULL,lambda="cv1se.dev",k=5){
  if(model=="binomial"){
      # need to use tryCatch to avoid lambda.max errors
      try1 <- tryCatch({
        fit <- glmsmurf(formula=formula, family=binomial(link="logit"), data=data,
                        weights=size[misspoi], pen.weights="glm.stand", lambda=lambda, 
                        control=list(lambda.length=20L, k=k, ncores=1));
        TRUE
      }, error=function(e) {
        message("Failed determining the maximum of lambda, run gaussian model instead")
        fit <-glmsmurf(formula=formula, family=gaussian(), data=data,
                pen.weights="glm.stand", lambda="cv1se.dev",control=list(lambda.length=20L, k=k, ncores=1))});
    }
  if(model=="gaussian"){
       fit <- glmsmurf(formula=formula, family=gaussian(), data=data,
                        pen.weights="glm.stand", lambda="cv1se.dev", 
                        control=list(lambda.length=20L, k=k, ncores=1))
  } 
  return(fit)
}

## Wilcoxon Rank Sum test ###################################
# res <- pairwise.wilcox.test(rent$ratio,rent$X, p.adjust.method ="BH")
# adj<-as.data.frame(res$p.value)[lower.tri(res$p.value, diag = T)]
# sigdif<-which(adj<0.05)
# ct_wilcoxon[poi,sigdif]<-1
# a<-which(res$p.value<0.05)
# asso<-cbind(colnames(adj)[ceiling(a/3)],rownames(adj)[ifelse(a%%3==0,3,a%%3)])
# adjm<-matrix(TRUE,nrow = 4,ncol = 4)
