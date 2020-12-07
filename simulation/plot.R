### Frame 1 simulation ##############################################
library(ggplot2)
library(tidyverse)
l <- lapply(list.files(path="/proj/milovelab/mu/SC-ASE/simulation/csv", pattern="^n.*csv$", full=TRUE), function(f) {
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
## Frame 2 simulation ###################
colnames(dat) <- c("type","ARI_mcl","ARI","cl","time")

ggplot(dat[which(dat$cl=="NAI"&dat$type=="bin"),], aes(x=ARI_mcl)) +
  geom_bar(stat="count", position="dodge", width=.075) 
  # +facet_grid(n ~ cnt, labeller = label_both)

ggplot(dat, aes(x=ARI,fill=type)) +
  geom_histogram(position = "dodge") +
  facet_wrap(vars(cl), labeller = label_both)

#############################################################
## Real data: Deng ####################### 
load(file = "/proj/milovelab/mu/SC-ASE/data/deng.rda")
load("../../data/postprocess.rda")

## time
time<-as_tibble(do.call(rbind, ans[seq(1,length(ans), by = 4)]))
colnames(time)<-c("bin","gau","binboot","gauboot")
time2<-time %>%cbind(cluster=1:16) %>%gather(key="model",value = "time",bin:gauboot)
time2 %>% group_by(model) %>% summarize(sumtime=sum(time)) %>% data.frame
ggplot(data = time2, mapping = aes(x = model, y=time)) + 
  geom_point()+
  ylim(0, 200)+
  facet_wrap(vars(cluster), labeller = label_both)

## estimator
ratio<-as_tibble(do.call(rbind, ans[seq(2,length(ans), by = 4)]))
colnames(ratio)<-cell_meta_unique
ratio2<-tibble(ratio,cluster=rep(1:16,each=2),model=rep(c("bin","gau"),16))
ratio3<-ratio2 %>% gather(key="type",value = "ratio",zy:lateblast)
ratio3[which(ratio3[,"ratio"]>1),"ratio"]<-1
ratio3[which(ratio3[,"ratio"]<0),"ratio"]<-0

ggplot(data = ratio3, mapping = aes(x = type, y=ratio,col=factor(model))) + 
  geom_point()+
  facet_wrap(vars(cluster), labeller = label_both)+
  scale_x_discrete(name ="Cell types", 
                   limits=cell_meta_unique)+
  theme(axis.text.x = element_text(angle=45))

## confidence interval
ci<-as_tibble(do.call(rbind, ans[seq(3,length(ans), by = 4)]))
colnames(ci)<-cell_meta_unique

ci2<-tibble(ci,cluster=rep(1:16,each=4),bound=rep(c("bin","bin","gau","gau"),16))
ci3<-ci2 %>% gather(key="type",value = "ci",zy:lateblast)
ci3[which(ci3[,"ci"]>1),"ci"]<-1
ci3[which(ci3[,"ci"]<0),"ci"]<-0

ggplot(data = ci3, mapping = aes(x = type, y=ci,col=factor(bound))) + 
  geom_point()+
  facet_wrap(vars(cluster), labeller = label_both)+
  scale_x_discrete(name ="Cell types", 
                   limits=cell_meta_unique)+
  theme(axis.text.x = element_text(angle=45))

## consensus partation ############
consens.par<-as_tibble(do.call(rbind, ans[seq(4,length(ans), by = 4)]))
colnames(consens.par)<-cell_meta_unique
consens.par2<-tibble(consens.par,cluster=rep(1:16,each=2),model=rep(c("bin","gau"),16))
