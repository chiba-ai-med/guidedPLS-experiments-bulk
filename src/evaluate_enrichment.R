source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
method <- args[2]
outfile <- args[3]

load(paste0("data/", sample, ".RData"))
infile <- paste0("output/", sample, "/", method, ".RData")
load(infile)

# Extract loadings
if(method == "guided_pls"){
    loadings_at <- resgPLS$loading1
    loadings_pj <- resgPLS$loading2
    rownames(loadings_at) <- rownames(X_At)
    rownames(loadings_pj) <- rownames(X_Pj)
}
if(method == "guided_pca"){
    loadings_at <- resgPCA_At$u
    loadings_pj <- resgPCA_Pj$u
}
if(method == "pca"){
    loadings_at <- res_pca_At$u[, seq(4)]
    loadings_pj <- res_pca_Pj$u[, seq(4)]
}
if(method == "rbh_pls"){
    loadings_at <- scaled_X_At_RBH %*% resgPLS$u[, seq(4)]
    loadings_pj <- scaled_X_Pj_RBH %*% resgPLS$v[, seq(4)]
    rownames(loadings_at) <- rownames(scaled_X_At_RBH)
    rownames(loadings_pj) <- rownames(scaled_X_Pj_RBH)
}

# Binarization (top 500 positive and negative)
loadings_at_pos <- apply(loadings_at, 2, function(x){
    out <- rep(0, length(x))
    names(out) <- names(x)
    out[which(rank(-x) <= 500)] <- 1
    out
})
loadings_at_neg <- apply(loadings_at, 2, function(x){
    out <- rep(0, length(x))
    names(out) <- names(x)
    out[which(rank(x) <= 500)] <- 1
    out
})
loadings_pj_pos <- apply(loadings_pj, 2, function(x){
    out <- rep(0, length(x))
    names(out) <- names(x)
    out[which(rank(-x) <= 500)] <- 1
    out
})
loadings_pj_neg <- apply(loadings_pj, 2, function(x){
    out <- rep(0, length(x))
    names(out) <- names(x)
    out[which(rank(x) <= 500)] <- 1
    out
})

# Output directory
outdir <- paste0("output/", sample, "/enrichment/", method)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Enrichment Analysis
GOType <- c("BP", "MF", "CC")

# At enrichment
if(method == "rbh_pls"){
    expmat_at <- scaled_X_At_RBH
} else {
    expmat_at <- scaled_X_At
}
GOList_At <- list(At_GO_BP, At_GO_MF, At_GO_CC)
for(i in seq(ncol(loadings_at))){
    for(j in seq(3)){
        outGO_pos <- .EnrichLoadings(loadings_at_pos[,i], GOList_At[[j]], expmat_at, 0.1)
        outGO_neg <- .EnrichLoadings(loadings_at_neg[,i], GOList_At[[j]], expmat_at, 0.1)
        write.table(outGO_pos,
            paste0(outdir, "/At_", GOType[j], "_pos", i, ".txt"),
            quote=FALSE, row.names=FALSE, sep="\t")
        write.table(outGO_neg,
            paste0(outdir, "/At_", GOType[j], "_neg", i, ".txt"),
            quote=FALSE, row.names=FALSE, sep="\t")
    }
}

# Pj enrichment
if(method == "rbh_pls"){
    expmat_pj <- scaled_X_Pj_RBH
} else {
    expmat_pj <- scaled_X_Pj
}
GOList_Pj <- list(Pj_GO_BP, Pj_GO_MF, Pj_GO_CC)
for(i in seq(ncol(loadings_pj))){
    for(j in seq(3)){
        outGO_pos <- .EnrichLoadings(loadings_pj_pos[,i], GOList_Pj[[j]], expmat_pj, 0.1)
        outGO_neg <- .EnrichLoadings(loadings_pj_neg[,i], GOList_Pj[[j]], expmat_pj, 0.1)
        write.table(outGO_pos,
            paste0(outdir, "/Pj_", GOType[j], "_pos", i, ".txt"),
            quote=FALSE, row.names=FALSE, sep="\t")
        write.table(outGO_neg,
            paste0(outdir, "/Pj_", GOType[j], "_neg", i, ".txt"),
            quote=FALSE, row.names=FALSE, sep="\t")
    }
}

# Output
file.create(outfile)
