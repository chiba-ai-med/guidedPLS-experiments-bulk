source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

methods <- c("guided_pls", "rbh_pls", "guided_pca", "pca")

####################################
# Build gold-standard gene sets from DEG
# Use strict DEG (FDR < 0.05, |logFC| > 1) as "true response genes"
####################################
.extract_deg_genes <- function(sample, species, condition, direction="both"){
    f <- paste0("output/", sample, "/deg/", species, "_", condition, ".RData")
    if(!file.exists(f)) return(character(0))
    load(f)
    tmp <- deg@.Data[[1]]
    sig <- tmp[tmp$FDR < 0.05 & abs(tmp$logFC) > 1, ]
    if(direction == "up") return(rownames(sig[sig$logFC > 0, ]))
    if(direction == "down") return(rownames(sig[sig$logFC < 0, ]))
    rownames(sig)
}

# Define gold sets per dataset
if(sample == "parasitism1"){
    gold_sets <- list(
        parasm_At = .extract_deg_genes(sample, "at", "parasm"),
        parasm_Pj = .extract_deg_genes(sample, "pj", "parasm"),
        wol_At = .extract_deg_genes(sample, "at", "wol"),
        wol_Pj = .extract_deg_genes(sample, "pj", "wol"),
        time1d_At = .extract_deg_genes(sample, "at", "1d"),
        time1d_Pj = .extract_deg_genes(sample, "pj", "1d"),
        time7d_At = .extract_deg_genes(sample, "at", "7d"),
        time7d_Pj = .extract_deg_genes(sample, "pj", "7d")
    )
}
if(sample == "grafting"){
    gold_sets <- list(
        graft_At = .extract_deg_genes(sample, "at", "graft"),
        graft_Pj = .extract_deg_genes(sample, "pj", "graft"),
        time1d_At = .extract_deg_genes(sample, "at", "1d"),
        time1d_Pj = .extract_deg_genes(sample, "pj", "1d"),
        time7d_At = .extract_deg_genes(sample, "at", "7d"),
        time7d_Pj = .extract_deg_genes(sample, "pj", "7d"),
        time14d_At = .extract_deg_genes(sample, "at", "14d"),
        time14d_Pj = .extract_deg_genes(sample, "pj", "14d")
    )
}

cat("Gold set sizes:\n")
for(nm in names(gold_sets)){
    cat(sprintf("  %s: %d genes\n", nm, length(gold_sets[[nm]])))
}

####################################
# Evaluate each method
####################################
topK_values <- c(100, 200, 500, 1000)
results <- list()

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

    # Max absolute loading per gene across all components
    max_loading_at <- apply(abs(loadings_at), 1, max)
    max_loading_pj <- apply(abs(loadings_pj), 1, max)

    for(gs_name in names(gold_sets)){
        gold <- gold_sets[[gs_name]]
        if(length(gold) == 0) next

        # Determine species
        if(grepl("_At$", gs_name)){
            max_loading <- max_loading_at
        } else {
            max_loading <- max_loading_pj
        }

        # Restrict gold to genes in loading space
        gold_in <- intersect(gold, names(max_loading))
        if(length(gold_in) == 0) next

        for(K in topK_values){
            top_genes <- names(sort(max_loading, decreasing=TRUE))[seq(min(K, length(max_loading)))]
            tp <- length(intersect(top_genes, gold_in))
            prec <- tp / length(top_genes)
            rec <- tp / length(gold_in)
            f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)

            results[[length(results)+1]] <- data.frame(
                method=m, gold_set=gs_name, K=K,
                gold_size=length(gold_in), tp=tp,
                precision=prec, recall=rec, f1=f1, auroc=NA)
        }

        # AUROC
        labels <- as.numeric(names(max_loading) %in% gold_in)
        ord <- order(max_loading, decreasing=TRUE)
        labels_sorted <- labels[ord]
        n_pos <- sum(labels_sorted)
        n_neg <- length(labels_sorted) - n_pos
        if(n_pos > 0 && n_neg > 0){
            tpr <- cumsum(labels_sorted) / n_pos
            fpr <- cumsum(1 - labels_sorted) / n_neg
            auroc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
        } else {
            auroc <- NA
        }

        results[[length(results)+1]] <- data.frame(
            method=m, gold_set=gs_name, K=NA,
            gold_size=length(gold_in), tp=NA,
            precision=NA, recall=NA, f1=NA, auroc=auroc)
    }
}

result_df <- do.call(rbind, results)

# Save
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
write.csv(result_df, outfile, row.names=FALSE)

####################################
# Summary plot: F1@500 by method and gold set
####################################
df_f1 <- result_df[!is.na(result_df$K) & result_df$K == 500, ]
df_f1$method <- factor(df_f1$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

g1 <- ggplot(df_f1, aes(x=gold_set, y=f1, fill=method)) +
    geom_bar(stat="identity", position="dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=11),
          text = element_text(size=13)) +
    scale_fill_brewer(palette="Set2") +
    ylab("F1 @K=500") +
    xlab("Gold standard gene set (DEG-derived)") +
    ggtitle(paste0("DEG recovery F1@500: ", sample))

outdir <- paste0("plot/", sample)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
ggsave(file=paste0(outdir, "/deg_gold_f1.png"), plot=g1, width=12, height=6)

####################################
# AUROC by method and gold set
####################################
df_auroc <- result_df[is.na(result_df$K), ]
if(nrow(df_auroc) > 0 && "auroc" %in% colnames(df_auroc)){
    df_auroc <- df_auroc[!is.na(df_auroc$auroc), ]
    df_auroc$method <- factor(df_auroc$method,
        levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

    g2 <- ggplot(df_auroc, aes(x=gold_set, y=auroc, fill=method)) +
        geom_bar(stat="identity", position="dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=11),
              text = element_text(size=13)) +
        scale_fill_brewer(palette="Set2") +
        ylab("AUROC") +
        xlab("Gold standard gene set (DEG-derived)") +
        geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
        ggtitle(paste0("DEG recovery AUROC: ", sample))

    ggsave(file=paste0(outdir, "/deg_gold_auroc.png"), plot=g2, width=12, height=6)
}

####################################
# Mean F1@500 across all gold sets per method
####################################
mean_f1 <- aggregate(f1 ~ method, data=df_f1, FUN=mean)
mean_f1$method <- factor(mean_f1$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

g3 <- ggplot(mean_f1, aes(x=method, y=f1, fill=method)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          text = element_text(size=14)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Mean F1 @K=500 across DEG gold sets") +
    xlab("Method") +
    ggtitle(paste0("Overall DEG recovery: ", sample))

ggsave(file=paste0(outdir, "/deg_gold_mean_f1.png"), plot=g3, width=6, height=5)
