source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

# PCA (At)
res_pca_At <- svd(scaled_X_At)

# PCA (Pj)
res_pca_Pj <- svd(scaled_X_Pj)

# Score (samples in rows)
score_At <- t(scaled_X_At) %*% res_pca_At$u[, seq(4)]
score_Pj <- t(scaled_X_Pj) %*% res_pca_Pj$u[, seq(4)]

# Output
save(res_pca_At, res_pca_Pj, score_At, score_Pj, file=outfile)
