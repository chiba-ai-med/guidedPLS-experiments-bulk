source("src/Functions.R")
library("viridisLite")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

methods <- c("guided_pls", "rbh_pls", "guided_pca", "pca")

####################################
# Load DEG results
####################################
.load_deg_table <- function(sample, species, conditions, counts_mat){
    qtable <- c()
    for(cond in conditions){
        filename <- paste0("output/", sample, "/deg/", species, "_", cond, ".RData")
        if(!file.exists(filename)) next
        load(filename)
        pos_vec <- rep(0, nrow(counts_mat))
        neg_vec <- rep(0, nrow(counts_mat))
        names(pos_vec) <- rownames(counts_mat)
        names(neg_vec) <- rownames(counts_mat)
        tmp <- deg@.Data[[1]]
        pos <- head(rownames(tmp)[which(tmp$logFC > 0)], 0.1*nrow(counts_mat))
        neg <- head(rownames(tmp)[which(tmp$logFC < 0)], 0.1*nrow(counts_mat))
        pos_vec[pos] <- 1
        neg_vec[neg] <- 1
        qtable <- cbind(qtable, pos_vec, neg_vec)
    }
    cnames <- unlist(lapply(conditions, function(x) paste0(x, c("+", "-"))))
    colnames(qtable) <- cnames
    qtable
}

qtable_at <- .load_deg_table(sample, "at", deg_conditions, counts_At)
qtable_pj <- .load_deg_table(sample, "pj", deg_conditions, counts_Pj)

####################################
# For each method, compute Jaccard and extract best factor-dim correspondence
####################################
all_results <- list()

for(m in methods){
    infile <- paste0("output/", sample, "/", m, ".RData")
    if(!file.exists(infile)) next
    load(infile)

    # Extract loadings
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
        rownames(loadings_at) <- rownames(scaled_X_At_RBH)
        rownames(loadings_pj) <- rownames(scaled_X_Pj_RBH)
    }

    # Positive/Negative split
    loadings_at_split <- loadings_at[, unlist(lapply(seq(ncol(loadings_at)), function(x) rep(x, 2)))]
    loadings_pj_split <- loadings_pj[, unlist(lapply(seq(ncol(loadings_pj)), function(x) rep(x, 2)))]
    for(i in seq(ncol(loadings_at_split))){
        if(i %% 2 == 1){
            loadings_at_split[,i] <- -loadings_at_split[,i]
            loadings_pj_split[,i] <- -loadings_pj_split[,i]
        }
    }
    colnames(loadings_at_split) <- unlist(lapply(seq(ncol(loadings_at_split)/2),
        function(x) paste0("Dim", x, c("+", "-"))))
    colnames(loadings_pj_split) <- unlist(lapply(seq(ncol(loadings_pj_split)/2),
        function(x) paste0("Dim", x, c("+", "-"))))

    # Binarize top 10%
    .binarize <- function(mat){
        apply(mat, 2, function(x){
            out <- rep(0, length(x))
            names(out) <- names(x)
            out[which(rank(x) <= 0.1*length(x))] <- 1
            out
        })
    }
    bin_at <- .binarize(loadings_at_split)
    bin_pj <- .binarize(loadings_pj_split)

    # Align for RBH
    if(m == "rbh_pls"){
        common_at <- intersect(rownames(qtable_at), rownames(bin_at))
        qt_at <- qtable_at[common_at, ]
        b_at <- bin_at[common_at, ]
        common_pj <- intersect(rownames(qtable_pj), rownames(bin_pj))
        qt_pj <- qtable_pj[common_pj, ]
        b_pj <- bin_pj[common_pj, ]
    } else {
        qt_at <- qtable_at
        b_at <- bin_at
        qt_pj <- qtable_pj
        b_pj <- bin_pj
    }

    # Jaccard (At)
    numer_at <- t(qt_at) %*% b_at
    denom_at <- outer(colSums(qt_at), colSums(b_at), "+") - numer_at
    jaccard_at <- numer_at / denom_at

    # Jaccard (Pj)
    numer_pj <- t(qt_pj) %*% b_pj
    denom_pj <- outer(colSums(qt_pj), colSums(b_pj), "+") - numer_pj
    jaccard_pj <- numer_pj / denom_pj

    # For each DEG condition, find the best matching Dim and its Jaccard
    for(cond in rownames(jaccard_at)){
        best_dim_at <- colnames(jaccard_at)[which.max(jaccard_at[cond, ])]
        best_jac_at <- max(jaccard_at[cond, ])
        best_dim_pj <- colnames(jaccard_pj)[which.max(jaccard_pj[cond, ])]
        best_jac_pj <- max(jaccard_pj[cond, ])

        all_results[[length(all_results)+1]] <- data.frame(
            method=m, deg_condition=cond,
            best_dim_At=best_dim_at, jaccard_At=best_jac_at,
            best_dim_Pj=best_dim_pj, jaccard_Pj=best_jac_pj)
    }

    # Mean Jaccard across all conditions (summary metric)
    mean_jac_at <- mean(apply(jaccard_at, 1, max))
    mean_jac_pj <- mean(apply(jaccard_pj, 1, max))
    all_results[[length(all_results)+1]] <- data.frame(
        method=m, deg_condition="MEAN_BEST",
        best_dim_At="", jaccard_At=mean_jac_at,
        best_dim_Pj="", jaccard_Pj=mean_jac_pj)

    # Factor separation score: how many distinct Dims are used as best match
    unique_dims_at <- length(unique(apply(jaccard_at, 1, which.max)))
    unique_dims_pj <- length(unique(apply(jaccard_pj, 1, which.max)))
    all_results[[length(all_results)+1]] <- data.frame(
        method=m, deg_condition="N_DISTINCT_DIMS",
        best_dim_At=as.character(unique_dims_at), jaccard_At=NA,
        best_dim_Pj=as.character(unique_dims_pj), jaccard_Pj=NA)
}

result_df <- do.call(rbind, all_results)

# Save
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
write.csv(result_df, outfile, row.names=FALSE)

####################################
# Summary plot: Mean best Jaccard by method
####################################
summary_df <- result_df[result_df$deg_condition == "MEAN_BEST", ]
summary_long <- rbind(
    data.frame(method=summary_df$method, species="At", mean_jaccard=summary_df$jaccard_At),
    data.frame(method=summary_df$method, species="Pj", mean_jaccard=summary_df$jaccard_Pj))
summary_long$method <- factor(summary_long$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

g <- ggplot(summary_long, aes(x=method, y=mean_jaccard, fill=method)) +
    geom_bar(stat="identity") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          text = element_text(size=14)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Mean best Jaccard (DEG x Loading)") +
    xlab("Method") +
    ggtitle(paste0("Factor-DEG correspondence: ", sample))

outdir <- paste0("plot/", sample)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
ggsave(file=paste0(outdir, "/jaccard_comparison.png"), plot=g, width=8, height=5)

####################################
# Distinct dims plot
####################################
dims_df <- result_df[result_df$deg_condition == "N_DISTINCT_DIMS", ]
dims_long <- rbind(
    data.frame(method=dims_df$method, species="At",
        n_dims=as.numeric(dims_df$best_dim_At)),
    data.frame(method=dims_df$method, species="Pj",
        n_dims=as.numeric(dims_df$best_dim_Pj)))
dims_long$method <- factor(dims_long$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

g2 <- ggplot(dims_long, aes(x=method, y=n_dims, fill=method)) +
    geom_bar(stat="identity") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          text = element_text(size=14)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Number of distinct best-matching Dims") +
    xlab("Method") +
    ggtitle(paste0("Factor separation: ", sample))

ggsave(file=paste0(outdir, "/factor_separation.png"), plot=g2, width=8, height=5)
