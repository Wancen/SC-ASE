# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(version = "3.11")
# BiocManager::install("splatter")
setwd("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/workflow/simulation")
library(splatter)
library("scater")
# create a default set of parameters for sim
params <- newSplatParams()
params
# params<-setParam(params, "mean.shape", 30)
params <- setParam(params, "lib.loc", 14)
params <- setParam(params, "lib.scale", 0.3)
# Set multiple parameters at once (using a list)
# # Set multiple parameters at once (using a list)
# params <- setParams(params, update = list(nGenes = 8000, mean.rate = 0.5))
sim <- splatSimulate(params,group.prob=c(.5, .5),
                     method="groups",
                     nGenes=103639,
                     batchCells=60,
                     seed=123)
sim <- logNormCounts(sim)
sim <- runPCA(sim)
plotPCA(sim, colour_by = "Group")
# order the cells by their true grouping
sim <- sim[,order(sim$Group)]
# Information about cells
colData(sim)
rowData(sim)
# label the DE genes
mcols(sim)$de <- with(mcols(sim), DEFacGroup2/DEFacGroup1 != 1)

# how many?
table(mcols(sim)$de)

# how many cells?
table(sim$Group)

total <- counts(sim)

hist(log(total[total>5]), xlab="log10 counts > 5")

# simulate 80% 0.5,5% 1, 5% 0, 5% 0.25, 5% 0.75
ref<-c(rep(0.5,16),1,0,0.25,0.75)
p<-matrix(sample(ref,103639*30,replace = T),ncol = 30)
table(p)/(103639*30)
# 0       0.25        0.5       0.75          1 
# 0.04996333 0.05061000 0.79985667 0.04946667 0.05010333 
hist(p)
# get_mat<-function(x,p){
#   rbinom(1,x,p)
# }

mat<-matrix(rbinom(103639*30,total,p),ncol=30)
# idx <- rowSums(total >= 3)>=10
# table(idx)
# a <- as.vector(mat[idx,]/total[idx,])
a<-mat/total
pat<-total-mat
h <- hist(a, breaks = 100, plot=FALSE)
h$counts=h$counts/sum(h$counts)
plot(h)
h <- hist(stem_salmon_ratio, breaks = 100, plot=FALSE)

library(polyester)
simulate_experiment_countmat(fasta="mat.transcripts.fa",
                    readmat = mat,
                    outdir="mat",
                    paired=FALSE,
                    strand_specific=FALSE,
                    seed=1)
simulate_experiment_countmat(fasta="pat.transcripts.fa",
                             readmat = pat,
                             outdir="pat",
                             paired=FALSE,
                             strand_specific=FALSE,
                             seed=1)
library(Biostrings)
library(R.utils)
for (allel in c("mat","pat")) {
for (i in 1:30) {
  cat(i)
  ii<-if(i<10) paste0("0",i) else i
reads1<-readDNAStringSet(sprintf("%s/sample_%s.fasta",allel,ii))
set.seed(1)
idx<-sample(length(reads1))
reads1<-reads1[idx]
file1<-sprintf("%s/shuffle_%s.fasta",allel,ii)
writeXStringSet(reads1,file=file1)
gzip(file1,overwrite=T)
}
}

