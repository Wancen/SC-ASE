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
## Frame 2 simulation ###################
dat <- do.call(rbind, l)
colnames(dat) <- c("type","ARI_cl","ARI","cl","time1","n","cnt")
dat %>% group_by(cnt, n) %>% summarize(count=n()) %>% data.frame
ggplot(dat[which(dat$type=="bin"&dat$cnt==50),], aes(x=ARI_cl,fill=type)) +
  geom_bar(stat="count", position="dodge", width=.075) 

ggplot(dat, aes(x=ARI,fill=type)) +
  geom_bar(stat="count", position="dodge", width=.075) +
  facet_grid(cl ~ cnt, labeller = label_both)
