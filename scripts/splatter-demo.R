library(splatter)
# create a default set of parameters for sim
params <- newSplatParams()
# simulate two groups, 20 cells, 1000 "genes"
sim <- splatSimulate(params,
                     group.prob=c(.5, .5),
                     method="groups",
                     nGenes=1000,
                     batchCells=20,
                     seed=123)

# order the cells by their true grouping
sim <- sim[,order(sim$Group)]

# label the DE genes
mcols(sim)$de <- with(mcols(sim), DEFacGroup2/DEFacGroup1 != 1)

# how many?
table(mcols(sim)$de)

# how many cells?
table(sim$Group)

# look at counts
assays(sim)[["counts"]][1:5,1:5]

# total counts
colSums(assays(sim)[["counts"]])

# total counts over the expected library size
plot(sim$ExpLibSize, colSums(assays(sim)[["counts"]]));abline(0,1)

# heatmap
library(rafalib)
imagemat(log10(assays(sim)[["counts"]]+1))
