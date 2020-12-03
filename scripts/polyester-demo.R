library(polyester)
# length of simulated txps
len <- 2000
# number of txps
n <- 500
# make some fake txps
txps <- replicate(n, {
  paste(sample(c("A","C","G","T"), len, TRUE), collapse="")
})
names(txps) <- paste0("txp-",1:n)
out <- as.vector( rbind(paste0(">",names(txps)), txps) )
# write them as FASTA format
write(out, file="transcripts.fa")

# make some simulated counts
sim_counts <- rnbinom(n=n, size=1/5, mu=10)

# make sure we don't blow up computer
stopifnot(sum(sim_counts) < 30e6)

# the NB dispersion
disp <- .01
# some simulated fold changes (here, all null)
fold_changes <- matrix(c(rep(1,n),rep(1,n)),ncol=2)

# counting number of transcripts...
# grep '^>' gencode.v32.transcripts.fa | wc -l
# 227462

# write the paired-end reads
simulate_experiment(fasta="transcripts.fa",
                    outdir="test",
                    num_reps=c(1,1),
                    reads_per_transcript=sim_counts,
                    size=1/disp,
                    fold_changes=fold_changes,
                    fraglen=300,
                    fragsd=50,
                    readlen=100,
                    seed=123)

# look at the simulated counts
load("test/sim_counts_matrix.rda")
counts_matrix
plot(sim_counts, counts_matrix[,1])
points(sim_counts, counts_matrix[,2], col="red")
abline(0,1)
