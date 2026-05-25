source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

# guided-PLS
k <- min(4, ncol(Y_At_all), ncol(Y_Pj_all))
resgPLS <- .guidedPLS(t(scaled_X_At), t(scaled_X_Pj),
    Y_At_all, Y_Pj_all, k, cortest=TRUE)
score_At <- resgPLS$score1
score_Pj <- resgPLS$score2

# Output
save(resgPLS, score_At, score_Pj, file=outfile)
