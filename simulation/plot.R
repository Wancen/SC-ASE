library(ggplot2)
library(dplyr)
l <- lapply(list.files(path="/proj/milovelab/mu/SC-ASE/simulation/csv", pattern="*.csv", full=TRUE), function(f) {
  read.csv(f, header=FALSE)
})
dat <- do.call(rbind, l)
colnames(dat) <- c("type","ARI","time","n","cnt")
dat %>% group_by(cnt, n) %>% summarize(count=n()) %>% data.frame

ggplot(dat, aes(x=ARI,fill=type)) +
  geom_bar(stat="count", position="dodge", width=.075) +
  facet_grid(n ~ cnt, labeller = label_both)

ggplot(dat, aes(x=type, y=time)) +
  geom_boxplot(outlier.color=NA) +
  geom_jitter(width=.1) +
  facet_grid(n ~ cnt, labeller = label_both)

#############################################################
## Compare coefficient relative efficiency ###################
box<-cbind(msebin_mean,msebin2_mean,mselin_mean,mselin2_mean)
order<-vector(length = 32)
for (i in 1:8) {
  order[(4*i-3):(4*i)]<-c(i,i+8,i+16,i+24)
}
boxplot(sqrt(box[,order]),ylab="Average MSE",ylim=c(0.2,0.26))
abline(h=0, col=2)

ARI<-data.frame(Gaussian=ARI_lin,Wilcoxon=ARI_wilcoxon,Binomial=ARI_bin)
AMI<-ARI %>% gather(key=method,value=AMI,Gaussian:Binomial)

jpeg("violin.jpeg",width = 5, height = 4,units = "in",res=300) 
p <- ggplot(AMI, aes(x=method, y=AMI,fill=method)) + 
  geom_violin()+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ theme_classic()+labs(title="Violin Plot of AMI  by methods",x="Methods", y = "Adjusted Mutual Information")
p+ geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()