source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

####################################
# Load guided-PLS and RBH-PLS results
####################################
load(paste0("output/", sample, "/guided_pls.RData"))
gpls_loading_at <- resgPLS$loading1
gpls_loading_pj <- resgPLS$loading2
rownames(gpls_loading_at) <- rownames(X_At)
rownames(gpls_loading_pj) <- rownames(X_Pj)

load(paste0("output/", sample, "/rbh_pls.RData"))
rbh_loading_at <- scaled_X_At_RBH %*% resgPLS$u[, seq(4)]
rbh_loading_pj <- scaled_X_Pj_RBH %*% resgPLS$v[, seq(4)]
rownames(rbh_loading_at) <- rownames(scaled_X_At_RBH)
rownames(rbh_loading_pj) <- rownames(scaled_X_Pj_RBH)

####################################
# Top genes by max absolute loading
####################################
topK <- 500

gpls_max_at <- apply(abs(gpls_loading_at), 1, max)
gpls_max_pj <- apply(abs(gpls_loading_pj), 1, max)
gpls_top_at <- names(sort(gpls_max_at, decreasing=TRUE))[seq(topK)]
gpls_top_pj <- names(sort(gpls_max_pj, decreasing=TRUE))[seq(topK)]

rbh_max_at <- apply(abs(rbh_loading_at), 1, max)
rbh_max_pj <- apply(abs(rbh_loading_pj), 1, max)
rbh_top_at <- names(sort(rbh_max_at, decreasing=TRUE))[seq(min(topK, length(rbh_max_at)))]
rbh_top_pj <- names(sort(rbh_max_pj, decreasing=TRUE))[seq(min(topK, length(rbh_max_pj)))]

####################################
# Classify guided-PLS top genes
####################################
# Ortholog genes (present in RBH table)
ortho_genes_at <- orthtable$V1[orthtable$V1 %in% rownames(X_At)]
ortho_genes_pj <- orthtable$V2[orthtable$V2 %in% rownames(X_Pj)]

# guided-PLS top genes: ortholog vs non-ortholog
gpls_top_at_ortho <- intersect(gpls_top_at, ortho_genes_at)
gpls_top_at_nonortho <- setdiff(gpls_top_at, ortho_genes_at)
gpls_top_pj_ortho <- intersect(gpls_top_pj, ortho_genes_pj)
gpls_top_pj_nonortho <- setdiff(gpls_top_pj, ortho_genes_pj)

cat("=== Gene classification (top", topK, ") ===\n")
cat(sprintf("guided-PLS At: %d ortholog, %d non-ortholog (total %d)\n",
    length(gpls_top_at_ortho), length(gpls_top_at_nonortho), topK))
cat(sprintf("guided-PLS Pj: %d ortholog, %d non-ortholog (total %d)\n",
    length(gpls_top_pj_ortho), length(gpls_top_pj_nonortho), topK))
cat(sprintf("RBH-PLS At: %d genes (all orthologs by definition)\n", length(rbh_top_at)))
cat(sprintf("RBH-PLS Pj: %d genes (all orthologs by definition)\n", length(rbh_top_pj)))

# Overlap between guided-PLS ortholog genes and RBH-PLS top genes
overlap_at <- intersect(gpls_top_at_ortho, rbh_top_at)
overlap_pj <- intersect(gpls_top_pj_ortho, rbh_top_pj)
cat(sprintf("\nOverlap (guided-PLS ortho ∩ RBH-PLS): At=%d, Pj=%d\n",
    length(overlap_at), length(overlap_pj)))

####################################
# DEG validation: are non-ortholog genes biologically meaningful?
####################################
.get_deg_genes <- function(sample, species, condition){
    f <- paste0("output/", sample, "/deg/", species, "_", condition, ".RData")
    if(!file.exists(f)) return(character(0))
    load(f)
    tmp <- deg@.Data[[1]]
    rownames(tmp[tmp$FDR < 0.05 & abs(tmp$logFC) > 1, ])
}

# Get DEG union across all conditions
if(sample == "parasitism1"){
    conditions <- c("1d", "3d", "7d", "wol", "parasm")
} else {
    conditions <- c("0d", "1d", "2d", "3d", "4d", "5d", "6d", "7d", "14d", "graft")
}

deg_union_at <- unique(unlist(lapply(conditions, function(c) .get_deg_genes(sample, "at", c))))
deg_union_pj <- unique(unlist(lapply(conditions, function(c) .get_deg_genes(sample, "pj", c))))

cat(sprintf("\nDEG union size: At=%d, Pj=%d\n", length(deg_union_at), length(deg_union_pj)))

# DEG enrichment in each gene category
results <- list()

.deg_overlap <- function(gene_set, deg_set, all_genes, label){
    gene_set <- intersect(gene_set, all_genes)
    deg_in <- intersect(gene_set, deg_set)
    deg_rate <- length(deg_in) / max(1, length(gene_set))
    # Fisher's exact test
    a <- length(deg_in)
    b <- length(gene_set) - a
    c <- length(deg_set) - a
    d <- length(all_genes) - a - b - c
    ft <- fisher.test(matrix(c(a, b, c, d), nrow=2), alternative="greater")
    data.frame(category=label, n_genes=length(gene_set),
        n_deg=length(deg_in), deg_rate=deg_rate,
        odds_ratio=ft$estimate, pvalue=ft$p.value)
}

# At
results[[1]] <- .deg_overlap(gpls_top_at_ortho, deg_union_at, rownames(X_At),
    "gPLS_ortho_At")
results[[2]] <- .deg_overlap(gpls_top_at_nonortho, deg_union_at, rownames(X_At),
    "gPLS_nonortho_At")
results[[3]] <- .deg_overlap(rbh_top_at, deg_union_at, rownames(X_At),
    "RBH_At")

# Pj
results[[4]] <- .deg_overlap(gpls_top_pj_ortho, deg_union_pj, rownames(X_Pj),
    "gPLS_ortho_Pj")
results[[5]] <- .deg_overlap(gpls_top_pj_nonortho, deg_union_pj, rownames(X_Pj),
    "gPLS_nonortho_Pj")
results[[6]] <- .deg_overlap(rbh_top_pj, deg_union_pj, rownames(X_Pj),
    "RBH_Pj")

result_df <- do.call(rbind, results)
cat("\n=== DEG enrichment by gene category ===\n")
print(result_df)

####################################
# GO enrichment for non-ortholog genes (guided-PLS unique discovery)
####################################
outdir <- paste0("output/", sample, "/discovery")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

GOType <- c("BP", "MF", "CC")
GOList_At <- list(At_GO_BP, At_GO_MF, At_GO_CC)
GOList_Pj <- list(Pj_GO_BP, Pj_GO_MF, Pj_GO_CC)

# Binarize non-ortholog genes
bin_nonortho_at <- rep(0, nrow(X_At))
names(bin_nonortho_at) <- rownames(X_At)
bin_nonortho_at[gpls_top_at_nonortho] <- 1

bin_nonortho_pj <- rep(0, nrow(X_Pj))
names(bin_nonortho_pj) <- rownames(X_Pj)
bin_nonortho_pj[gpls_top_pj_nonortho] <- 1

for(j in seq(3)){
    outGO_at <- .EnrichLoadings(bin_nonortho_at, GOList_At[[j]], X_At, 0.1)
    outGO_pj <- .EnrichLoadings(bin_nonortho_pj, GOList_Pj[[j]], X_Pj, 0.1)
    write.table(outGO_at,
        paste0(outdir, "/nonortho_At_", GOType[j], ".txt"),
        quote=FALSE, row.names=FALSE, sep="\t")
    write.table(outGO_pj,
        paste0(outdir, "/nonortho_Pj_", GOType[j], ".txt"),
        quote=FALSE, row.names=FALSE, sep="\t")
}

####################################
# Save results
####################################
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
write.csv(result_df, outfile, row.names=FALSE)

####################################
# Plots
####################################
outdir_plot <- paste0("plot/", sample)

# 1. Venn-like bar: gene composition of guided-PLS top genes
comp_df <- data.frame(
    species = c("At", "At", "At", "Pj", "Pj", "Pj"),
    category = c("gPLS ortholog\n(shared w/ RBH-PLS)",
                 "gPLS non-ortholog\n(unique to gPLS)",
                 "RBH-PLS only",
                 "gPLS ortholog\n(shared w/ RBH-PLS)",
                 "gPLS non-ortholog\n(unique to gPLS)",
                 "RBH-PLS only"),
    count = c(
        length(gpls_top_at_ortho),
        length(gpls_top_at_nonortho),
        length(setdiff(rbh_top_at, gpls_top_at)),
        length(gpls_top_pj_ortho),
        length(gpls_top_pj_nonortho),
        length(setdiff(rbh_top_pj, gpls_top_pj))
    )
)
comp_df$category <- factor(comp_df$category,
    levels=c("gPLS ortholog\n(shared w/ RBH-PLS)",
             "gPLS non-ortholog\n(unique to gPLS)",
             "RBH-PLS only"))

g1 <- ggplot(comp_df, aes(x=category, y=count, fill=category)) +
    geom_bar(stat="identity") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13),
          legend.position="none") +
    scale_fill_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb")) +
    ylab(paste0("Number of top-", topK, " genes")) +
    xlab("") +
    ggtitle(paste0("Gene discovery composition: ", sample))

ggsave(file=paste0(outdir_plot, "/discovery_composition.png"),
    plot=g1, width=10, height=5)

# 2. DEG enrichment rate comparison
result_df$category <- factor(result_df$category,
    levels=c("gPLS_ortho_At", "gPLS_nonortho_At", "RBH_At",
             "gPLS_ortho_Pj", "gPLS_nonortho_Pj", "RBH_Pj"))
result_df$species <- ifelse(grepl("At$", result_df$category), "At", "Pj")
result_df$type <- gsub("_(At|Pj)$", "", result_df$category)
result_df$type <- factor(result_df$type,
    levels=c("gPLS_ortho", "gPLS_nonortho", "RBH"))

g2 <- ggplot(result_df, aes(x=type, y=deg_rate, fill=type)) +
    geom_bar(stat="identity") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13)) +
    scale_fill_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb")) +
    ylab("Fraction of genes that are DEG") +
    xlab("") +
    ggtitle(paste0("DEG enrichment by category: ", sample))

ggsave(file=paste0(outdir_plot, "/discovery_deg_rate.png"),
    plot=g2, width=8, height=5)

# 3. Significance (-log10 pvalue)
result_df$neg_log10_p <- -log10(result_df$pvalue + 1e-300)

g3 <- ggplot(result_df, aes(x=type, y=neg_log10_p, fill=type)) +
    geom_bar(stat="identity") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13)) +
    scale_fill_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb")) +
    ylab("-log10(p-value) Fisher's exact test") +
    xlab("") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
    ggtitle(paste0("DEG enrichment significance: ", sample))

ggsave(file=paste0(outdir_plot, "/discovery_significance.png"),
    plot=g3, width=8, height=5)
