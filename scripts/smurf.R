library("smurf")
library(dplyr)
#######################################
## Hierarchical clustering ############
###############################################
gene_dist<-dist(ratio4)
gene_hclust <- hclust(gene_dist, method = "ward.D2")
plot(gene_hclust, labels = FALSE)
abline(h = 6, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
# cut genes into 3 group 
gene_cluster<-cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe(name="gene",value = "cluster")

table(gene_cluster$cluster)
# ordered by dendrogram order to visualize
library(dendextend)
library(RColorBrewer)
palette(brewer.pal(8, "Dark2"))
dend<-as.dendrogram(gene_hclust)
o.dend <- order.dendrogram(dend)

labels_colors(dend) <- as.integer(gene_cluster$cluster[o.dend])
plot(dend)

## Visualize gene pattern per cluster
# have to transfrom wide format into long format

# extrace average allelic ration per cluster
gene1_mean<-colMeans(ratio4[which(gene_cluster$cluster==5),],na.rm = T)
count1_mean<-colMeans(cts2[which(gene_cluster$cluster==5),])

########################################################################
## detect cell type according allelic ratio ############################
#####################################################
rent<-tibble(ratio=gene1_mean,X=as.factor(rep(c(1,2,3,4),each=20)),cts=count1_mean)
levels(rent$X) <- c("type1", "type2", "type3","type4")
## Gaussian likelihood #####################
formu <- ratio ~ p(X, pen = "gflasso", refcat = "type1")
system.time(munich.fit2 <- glmsmurf(formula = formu, family=gaussian(), data = rent,
                                    pen.weights = "eq", lambda = "cv1se.mse",
                                    control = list(lambda.length = 50L),pen.weights.return=T)) #pen.weights="glm" will fused coefficients more
                                                                               # and more consistent with truth sometimes 
plot(munich.fit2)
plot_reest(munich.fit2)
summary(munich.fit2)
## Binomial likelihood ###########################
system.time(munich.fit <- glmsmurf(formula = formu, family=binomial(link = "logit"), data = rent, weights = count1_mean,
                       pen.weights = "eq", lambda = "cv1se.mse",
                       control = list(lambda.length = 50L),pen.weights.return=T)) #
munich.fit$lambda
plot_lambda(munich.fit)
plot(munich.fit)
plot_reest(munich.fit)
summary(munich.fit)



