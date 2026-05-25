source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

# guided-PCA: SVD of X %*% Y (genes x samples %*% samples x factors = genes x factors)
scaled_Y_At_all <- scale(Y_At_all, center=TRUE, scale=TRUE)
scaled_Y_Pj_all <- scale(Y_Pj_all, center=TRUE, scale=TRUE)

resgPCA_At <- svd(scaled_X_At %*% scaled_Y_At_all)
resgPCA_Pj <- svd(scaled_X_Pj %*% scaled_Y_Pj_all)

# Score: samples x components (t(X) %*% u = samples x genes %*% genes x components)
k_at <- min(4, ncol(resgPCA_At$u))
k_pj <- min(4, ncol(resgPCA_Pj$u))
score_At <- t(scaled_X_At) %*% resgPCA_At$u[, seq(k_at)]
score_Pj <- t(scaled_X_Pj) %*% resgPCA_Pj$u[, seq(k_pj)]

# Output
save(resgPCA_At, resgPCA_Pj, score_At, score_Pj, file=outfile)
