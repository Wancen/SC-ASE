library(rsample)
library(aricode)
### bagging ######################
# Helper packages
library(doParallel)  # for parallel backend to foreach
library(foreach)     # for parallel processing with for loops
library(pbapply)
library(clue)
# Create a parallel socket cluster
cl <- makeCluster(4) # use 4 workers
registerDoParallel(cl) # register the parallel backend

total_sub<-totalexp[gene_feat_mclust]
n1<-length(gene_feat_mclust)
nct=10

gene_ratio<-as.vector(unlist(ratio[gene_feat_mclust,]))
poina<-which(!is.na(gene_ratio))
gene_total<-as.vector(unlist(total_fil2[gene_feat_mclust,]))
data<-tibble(ratio=gene_ratio,X=factor(rep(celltype,each=length(gene_feat_mclust)),levels = cell_meta_unique),cts=gene_total)
data %>%
  group_by(X) %>%
  summarise(weighted_mean = weighted.mean(ratio, cts,na.rm=T))
formu <- ratio ~ p(X, pen = "gflasso")
start_time <- Sys.time()
fit_total <- glmsmurf(formula = formu, family=binomial(link = "logit"), data = data, weights = gene_total[poina],
                            pen.weights = "glm.stand", lambda = "cv1se.dev",control = list(lambda.length = 50L,k=10)) #
end_time <- Sys.time()
end_time - start_time
stopCluster(cl)
# bootstrap of cells
boot<-bootstraps(data,times = 100,"X")

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
# Fit model in parallel 
start_time <- Sys.time()
fit22 <- foreach(
   i=1:100, 
  .packages = c("smurf","tibble","rsample"), 
  .combine = 'comb',
  .multicombine=TRUE,
  .init=list(list(), list())
) %dopar% {
  data_b<-as.data.frame( boot$splits[[i]] )
  poina<-which(!is.na(data_b$ratio))
  # fit fused lasso model to bootstrap copy
  formu <- ratio ~ p(X, pen = "gflasso")
  fit_boot <- glmsmurf(formula = formu, family=binomial(link = "logit"), data = data_b, weights = data_b$cts[poina],
                  pen.weights = "glm.stand", lambda = fit_total$lambda) #
  coefbin<-coef_reest(fit_boot)
  coefbin[2:nct]<-coefbin[2:nct]+coefbin[1]
  bin_coef=exp(coefbin)/(1+exp(coefbin))
  list(match(coefbin, unique(coefbin)),bin_coef)
}
end_time <- Sys.time()
end_time - start_time
stopCluster(cl)

## consensus partation #################
cluster <- fit22[[1]]
consens_par<-cl_consensus(cl_ensemble(list=lapply(cluster, as.cl_hard_partition)),method = "SM")

class <- max.col(consens_par$.Data)
# a<-diceR::consensus_matrix(t(cluster2))
# row.names(a)<-cell_meta_unique
# colnames(a)<-cell_meta_unique
# dist<-as.dist(a)


# maxlab<-apply(cluster2,1,max)
# nmi<-pbsapply(1:100, function(i) {NMI(unlist(cluster2[i,]),kmean_clue$cluster)})
# nmi<-pbsapply(1:100, function(i) {NMI(unlist(cluster2[i,]),c(1,2,3,4,5,5,6,5,5,5))})
# nmi%*%maxlab/sum(maxlab)
# 
# 
# ANMI_calculation(cluster2,kmean_clue$cluster)
# ANMI_calculation(cluster2,bb)
## coef test ###########################
coef<-fit22[[2]]
coef<-do.call(rbind.data.frame,coef)
colnames(coef)<-cell_meta_unique
boxplot(coef)

# Add a column with your condition for the color
coef2<-coef %>% gather(key=cell,value=coeffici,zy:lateblast)
coef2[]
jpeg(file="boxplot.jpeg",width = 7, height = 5,units = "in",res=450)
coef2 %>%
  mutate(cell = factor(cell,levels = cell_meta_unique) )%>%
  ggplot( aes(x=cell, y=coeffici, fill=cell)) + 
  geom_boxplot()+ theme_classic()+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  ylab("allelic ratio")+
  xlab("class") +
  theme(legend.position="none") +
  xlab("")
dev.off()
## calculate confidence interval  ###############################
coef_est<-coef_reest(fit)
coef_est[2:nct]<-coef_est[2:nct]+coef_est[1]
ratio_est=exp(coef_est)/(1+exp(coef_est))
ratio_est_mat<-matrix(rep(ratio_est,each=100),ncol=10)
bias<-coef-ratio_est_mat
bias_order<-pbsapply(1:10, function(i) {bias[order(bias[,i]),i]})
lower_CI<-ratio_est-bias_order[95,]
upper_CI<-ratio_est-bias_order[5,]

## bias correction estimator ################
bias_boot<-colMeans(coef)-ratio_est
ratio_bc<-ratio_est-bias_boot
