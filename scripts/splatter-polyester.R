# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# BiocManager::install("splatter")
setwd("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/GRA/workflow/simulation")
library(splatter)
# create a default set of parameters for sim
params <- newSplatParams()
# simulate two groups, 20 cells, 1000 "genes"
sim <- splatSimulate(params,
                     method="single",
                     nGenes=103639,
                     batchCells=30,
                     seed=123)
total <- counts(sim)
# simulate 80% 0.5,5% 1, 5% 0, 5% 0.25, 5% 0.75
ref<-c(rep(0.5,16),1,0,0.25,0.75)
p<-matrix(sample(ref,103639*30,replace = T),ncol = 30)
table(p)/(103639*30)
# 0       0.25        0.5       0.75          1 
# 0.04996333 0.05061000 0.79985667 0.04946667 0.05010333 
get_mat<-function(x,p){
  rbinom(1,x,p)
}
mat<-total
for(i in 1:10000){
  for (j in 1:30) {
    mat[i,j]<-rbinom(1,total[i,j],p[i,j])
  }
}
pat<-total-mat

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

