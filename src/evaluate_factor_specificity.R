source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

methods <- c("guided_pls", "rbh_pls", "guided_pca", "pca")
topK <- 200

####################################
# Build DEG gene sets per condition
####################################
.extract_deg_genes <- function(sample, species, condition){
    f <- paste0("output/", sample, "/deg/", species, "_", condition, ".RData")
    if(!file.exists(f)) return(character(0))
    load(f)
    tmp <- deg@.Data[[1]]
    rownames(tmp[tmp$FDR < 0.05 & abs(tmp$logFC) > 1, ])
}

deg_sets_at <- lapply(deg_conditions, function(c) .extract_deg_genes(sample, "at", c))
names(deg_sets_at) <- deg_conditions
deg_sets_pj <- lapply(deg_conditions, function(c) .extract_deg_genes(sample, "pj", c))
names(deg_sets_pj) <- deg_conditions

cat("DEG set sizes:\n")
for(cond in deg_conditions){
    cat(sprintf("  %s: At=%d, Pj=%d\n", cond,
        length(deg_sets_at[[cond]]), length(deg_sets_pj[[cond]])))
}

####################################
# Per-dimension enrichment for each method
####################################
.fisher_enrichment <- function(top_genes, deg_genes, all_genes){
    top_genes <- intersect(top_genes, all_genes)
    deg_genes <- intersect(deg_genes, all_genes)
    a <- length(intersect(top_genes, deg_genes))
    b <- length(top_genes) - a
    c <- length(deg_genes) - a
    d <- length(all_genes) - a - b - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow=2), alternative="greater")
    c(neg_log10_p = -log10(ft$p.value + 1e-300),
      odds_ratio = as.numeric(ft$estimate),
      n_overlap = a, n_top = length(top_genes), n_deg = length(deg_genes))
}

results <- list()

for(m in methods){
    infile <- paste0("output/", sample, "/", m, ".RData")
    if(!file.exists(infile)) next
    load(infile)

    # Extract per-dimension loadings
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

    ndim <- ncol(loadings_at)

    for(d in seq(ndim)){
        # Top genes per dimension
        top_at <- names(sort(abs(loadings_at[, d]), decreasing=TRUE))[
            seq(min(topK, nrow(loadings_at)))]
        top_pj <- names(sort(abs(loadings_pj[, d]), decreasing=TRUE))[
            seq(min(topK, nrow(loadings_pj)))]

        for(cond in deg_conditions){
            # At enrichment
            res_at <- .fisher_enrichment(top_at, deg_sets_at[[cond]], rownames(X_At))
            results[[length(results)+1]] <- data.frame(
                method=m, dim=d, species="At", condition=cond,
                neg_log10_p=res_at["neg_log10_p"],
                odds_ratio=res_at["odds_ratio"],
                n_overlap=res_at["n_overlap"],
                n_top=res_at["n_top"],
                n_deg=res_at["n_deg"],
                row.names=NULL)

            # Pj enrichment
            res_pj <- .fisher_enrichment(top_pj, deg_sets_pj[[cond]], rownames(X_Pj))
            results[[length(results)+1]] <- data.frame(
                method=m, dim=d, species="Pj", condition=cond,
                neg_log10_p=res_pj["neg_log10_p"],
                odds_ratio=res_pj["odds_ratio"],
                n_overlap=res_pj["n_overlap"],
                n_top=res_pj["n_top"],
                n_deg=res_pj["n_deg"],
                row.names=NULL)
        }
    }
}

result_df <- do.call(rbind, results)

####################################
# Compute specificity metrics per method
####################################
cat("\n=== Factor Specificity Summary ===\n")

specificity_results <- list()

for(m in unique(result_df$method)){
    for(sp in c("At", "Pj")){
        sub <- result_df[result_df$method == m & result_df$species == sp, ]
        ndim <- max(sub$dim)

        # Build enrichment matrix: dim × condition
        enrich_mat <- matrix(0, nrow=ndim, ncol=length(deg_conditions))
        rownames(enrich_mat) <- paste0("Dim", seq(ndim))
        colnames(enrich_mat) <- deg_conditions

        for(d in seq(ndim)){
            for(j in seq_along(deg_conditions)){
                val <- sub$neg_log10_p[sub$dim == d & sub$condition == deg_conditions[j]]
                if(length(val) > 0) enrich_mat[d, j] <- val
            }
        }

        # Per-dimension: best condition and specificity
        best_cond <- apply(enrich_mat, 1, function(x) colnames(enrich_mat)[which.max(x)])
        max_val <- apply(enrich_mat, 1, max)

        # Specificity: max / (max + 2nd_max) — 1.0 = perfectly specific
        specificity <- apply(enrich_mat, 1, function(x){
            sorted <- sort(x, decreasing=TRUE)
            if(sorted[1] == 0) return(0)
            sorted[1] / (sorted[1] + sorted[2])
        })

        # N distinct conditions captured
        n_distinct <- length(unique(best_cond))

        # Cross-species: does At-dim-d and Pj-dim-d map to same condition?
        cat(sprintf("\n%s %s:\n", m, sp))
        for(d in seq(ndim)){
            cat(sprintf("  Dim%d -> %s (score=%.1f, specificity=%.2f)\n",
                d, best_cond[d], max_val[d], specificity[d]))
        }
        cat(sprintf("  N_DISTINCT conditions: %d/%d\n", n_distinct, ndim))
        cat(sprintf("  Mean specificity: %.3f\n", mean(specificity)))

        specificity_results[[length(specificity_results)+1]] <- data.frame(
            method=m, species=sp,
            n_distinct=n_distinct, n_dim=ndim,
            mean_specificity=mean(specificity),
            mean_max_score=mean(max_val))
    }
}

spec_df <- do.call(rbind, specificity_results)
cat("\n=== Specificity comparison ===\n")
print(spec_df)

####################################
# Cross-species concordance: do At and Pj dims map to same factors?
####################################
cat("\n=== Cross-species dimension concordance ===\n")
concordance_results <- list()

for(m in unique(result_df$method)){
    sub_at <- result_df[result_df$method == m & result_df$species == "At", ]
    sub_pj <- result_df[result_df$method == m & result_df$species == "Pj", ]
    ndim <- max(sub_at$dim)

    n_concordant <- 0
    for(d in seq(ndim)){
        best_at <- deg_conditions[which.max(
            sub_at$neg_log10_p[sub_at$dim == d])]
        best_pj <- deg_conditions[which.max(
            sub_pj$neg_log10_p[sub_pj$dim == d])]
        match <- best_at == best_pj
        if(match) n_concordant <- n_concordant + 1
        cat(sprintf("  %s Dim%d: At->%s, Pj->%s %s\n",
            m, d, best_at, best_pj, ifelse(match, "MATCH", "")))
    }
    concordance_results[[length(concordance_results)+1]] <- data.frame(
        method=m, n_concordant=n_concordant, n_dim=ndim,
        concordance_rate=n_concordant/ndim)
}

conc_df <- do.call(rbind, concordance_results)
cat("\n")
print(conc_df)

####################################
# Save
####################################
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
write.csv(result_df, outfile, row.names=FALSE)
write.csv(spec_df,
    paste0(dirname(outfile), "/factor_specificity_summary.csv"), row.names=FALSE)
write.csv(conc_df,
    paste0(dirname(outfile), "/factor_concordance.csv"), row.names=FALSE)

####################################
# Plots
####################################
library("viridisLite")
outdir <- paste0("plot/", sample)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# 1. Heatmaps: dimension × condition enrichment per method per species
for(m in unique(result_df$method)){
    for(sp in c("At", "Pj")){
        sub <- result_df[result_df$method == m & result_df$species == sp, ]
        sub$dim_label <- paste0("Dim", sub$dim)
        sub$dim_label <- factor(sub$dim_label,
            levels=paste0("Dim", seq(max(sub$dim))))
        sub$condition <- factor(sub$condition, levels=deg_conditions)

        g <- ggplot(sub, aes(x=condition, y=dim_label, fill=neg_log10_p)) +
            geom_tile(color="white") +
            geom_text(aes(label=round(neg_log10_p, 1)), size=3.5) +
            scale_fill_gradientn(colours=viridis(100), name="-log10(p)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle=45, hjust=1, size=11),
                  axis.text.y = element_text(size=11),
                  text = element_text(size=13)) +
            xlab("DEG condition") + ylab("Dimension") +
            ggtitle(paste0(m, " ", sp, ": dimension-factor enrichment"))

        ggsave(file=paste0(outdir, "/factor_spec_", m, "_", sp, ".png"),
            plot=g, width=max(6, length(deg_conditions)*0.8), height=4)
    }
}

# 2. Summary bar: mean specificity by method
spec_df$method <- factor(spec_df$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

g_spec <- ggplot(spec_df, aes(x=method, y=mean_specificity, fill=method)) +
    geom_bar(stat="identity") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Mean factor specificity") +
    xlab("") +
    ggtitle(paste0("Factor specificity: ", sample))

ggsave(file=paste0(outdir, "/factor_specificity.png"),
    plot=g_spec, width=8, height=5)

# 3. Cross-species concordance bar
conc_df$method <- factor(conc_df$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

g_conc <- ggplot(conc_df, aes(x=method, y=concordance_rate, fill=method)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Cross-species concordance rate") +
    xlab("") +
    ylim(0, 1) +
    ggtitle(paste0("Cross-species dimension concordance: ", sample))

ggsave(file=paste0(outdir, "/factor_concordance.png"),
    plot=g_conc, width=6, height=5)
