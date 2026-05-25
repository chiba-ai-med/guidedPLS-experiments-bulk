source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

####################################
# Methods to evaluate
####################################
methods <- c("guided_pls", "rbh_pls", "guided_pca", "pca")
topK_values <- c(100, 200, 500)

####################################
# Gold standard gene set: Procambium markers
####################################
gold_at <- intersect(procambium_At$GENEID, rownames(X_At))
gold_pj <- intersect(procambium_Pj$GENEID, rownames(X_Pj))

results <- list()

for(m in methods){
    infile <- paste0("output/", sample, "/", m, ".RData")
    if(!file.exists(infile)) next
    load(infile)

    # Extract loadings per method
    if(m == "guided_pls"){
        loadings_at <- resgPLS$loading1
        loadings_pj <- resgPLS$loading2
        rownames(loadings_at) <- rownames(X_At)
        rownames(loadings_pj) <- rownames(X_Pj)
    }
    if(m == "guided_pca"){
        loadings_at <- resgPCA_At$u
        loadings_pj <- resgPCA_Pj$u
        rownames(loadings_at) <- rownames(X_At)
        rownames(loadings_pj) <- rownames(X_Pj)
    }
    if(m == "pca"){
        loadings_at <- res_pca_At$u[, seq(4)]
        loadings_pj <- res_pca_Pj$u[, seq(4)]
        rownames(loadings_at) <- rownames(X_At)
        rownames(loadings_pj) <- rownames(X_Pj)
    }
    if(m == "rbh_pls"){
        loadings_at <- scaled_X_At_RBH %*% resgPLS$u[, seq(4)]
        loadings_pj <- scaled_X_Pj_RBH %*% resgPLS$v[, seq(4)]
    }

    # Compute max absolute loading across all components per gene
    max_loading_at <- apply(abs(loadings_at), 1, max)
    max_loading_pj <- apply(abs(loadings_pj), 1, max)

    for(K in topK_values){
        # Top-K genes by max absolute loading
        top_at <- names(sort(max_loading_at, decreasing=TRUE))[seq(min(K, length(max_loading_at)))]
        top_pj <- names(sort(max_loading_pj, decreasing=TRUE))[seq(min(K, length(max_loading_pj)))]

        # Precision, Recall, F1 for At
        tp_at <- length(intersect(top_at, gold_at))
        prec_at <- tp_at / length(top_at)
        rec_at <- tp_at / max(1, length(gold_at))
        f1_at <- ifelse(prec_at + rec_at > 0,
            2 * prec_at * rec_at / (prec_at + rec_at), 0)

        # Precision, Recall, F1 for Pj
        tp_pj <- length(intersect(top_pj, gold_pj))
        prec_pj <- tp_pj / length(top_pj)
        rec_pj <- tp_pj / max(1, length(gold_pj))
        f1_pj <- ifelse(prec_pj + rec_pj > 0,
            2 * prec_pj * rec_pj / (prec_pj + rec_pj), 0)

        results[[length(results)+1]] <- data.frame(
            method=m, species="At", K=K,
            precision=prec_at, recall=rec_at, f1=f1_at,
            gold_size=length(gold_at), tp=tp_at, auroc=NA)
        results[[length(results)+1]] <- data.frame(
            method=m, species="Pj", K=K,
            precision=prec_pj, recall=rec_pj, f1=f1_pj,
            gold_size=length(gold_pj), tp=tp_pj, auroc=NA)
    }

    # AUROC using loading magnitude as score
    # At
    labels_at <- as.numeric(names(max_loading_at) %in% gold_at)
    ord_at <- order(max_loading_at, decreasing=TRUE)
    labels_sorted_at <- labels_at[ord_at]
    tpr_at <- cumsum(labels_sorted_at) / max(1, sum(labels_sorted_at))
    fpr_at <- cumsum(1 - labels_sorted_at) / max(1, sum(1 - labels_sorted_at))
    auroc_at <- sum(diff(fpr_at) * (tpr_at[-1] + tpr_at[-length(tpr_at)]) / 2)

    # Pj
    labels_pj <- as.numeric(names(max_loading_pj) %in% gold_pj)
    ord_pj <- order(max_loading_pj, decreasing=TRUE)
    labels_sorted_pj <- labels_pj[ord_pj]
    tpr_pj <- cumsum(labels_sorted_pj) / max(1, sum(labels_sorted_pj))
    fpr_pj <- cumsum(1 - labels_sorted_pj) / max(1, sum(1 - labels_sorted_pj))
    auroc_pj <- sum(diff(fpr_pj) * (tpr_pj[-1] + tpr_pj[-length(tpr_pj)]) / 2)

    results[[length(results)+1]] <- data.frame(
        method=m, species="At", K=NA,
        precision=NA, recall=NA, f1=NA,
        gold_size=length(gold_at), tp=NA, auroc=auroc_at)
    results[[length(results)+1]] <- data.frame(
        method=m, species="Pj", K=NA,
        precision=NA, recall=NA, f1=NA,
        gold_size=length(gold_pj), tp=NA, auroc=auroc_pj)
}

result_df <- do.call(rbind, results)

# Output
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
write.csv(result_df, outfile, row.names=FALSE)
