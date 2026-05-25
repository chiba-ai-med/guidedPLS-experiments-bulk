source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

# PLS on RBH orthologs
resgPLS <- svd(t(scaled_X_At_RBH) %*% scaled_X_Pj_RBH)
score_At <- resgPLS$u[, seq(4)]
score_Pj <- resgPLS$v[, seq(4)]

# Output
save(resgPLS, score_At, score_Pj, file=outfile)
