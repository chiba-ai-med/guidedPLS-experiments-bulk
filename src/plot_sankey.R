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
    rownames(loadings_at) <- rownames(X_At)
    rownames(loadings_pj) <- rownames(X_Pj)
}
if(method == "pca"){
    loadings_at <- res_pca_At$u[, seq(min(4, ncol(res_pca_At$u)))]
    loadings_pj <- res_pca_Pj$u[, seq(min(4, ncol(res_pca_Pj$u)))]
    rownames(loadings_at) <- rownames(X_At)
    rownames(loadings_pj) <- rownames(X_Pj)
}
if(method == "rbh_pls"){
    loadings_at <- scaled_X_At_RBH %*% resgPLS$u[, seq(4)]
    loadings_pj <- scaled_X_Pj_RBH %*% resgPLS$v[, seq(4)]
    rownames(loadings_at) <- rownames(scaled_X_At_RBH)
    rownames(loadings_pj) <- rownames(scaled_X_Pj_RBH)
}

####################################
# Time-resolved loading analysis
# Track how top genes' DEG status changes across timepoints
####################################
if(sample == "parasitism1"){
    timepoints <- c("1d", "3d", "7d")
} else {
    timepoints <- c("0d", "1d", "2d", "3d", "4d", "5d", "6d", "7d", "14d")
}

topK <- 200
n_comp <- min(ncol(loadings_at), 4)

outdir <- paste0("plot/", sample, "/sankey/", method)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

for(comp in seq(n_comp)){
    top_genes_at <- names(sort(abs(loadings_at[, comp]), decreasing=TRUE))[seq(topK)]

    # Build transition matrix across timepoints
    categories <- c("up", "down", "ns")
    gene_status <- matrix("ns", nrow=length(top_genes_at), ncol=length(timepoints))
    rownames(gene_status) <- top_genes_at
    colnames(gene_status) <- timepoints

    for(t_idx in seq_along(timepoints)){
        tp <- timepoints[t_idx]
        f <- paste0("output/", sample, "/deg/at_", tp, ".RData")
        if(!file.exists(f)) next
        load(f)
        tmp <- deg@.Data[[1]]
        for(gene in top_genes_at){
            if(gene %in% rownames(tmp)){
                lfc <- tmp[gene, "logFC"]
                fdr <- tmp[gene, "FDR"]
                if(fdr < 0.05){
                    gene_status[gene, t_idx] <- ifelse(lfc > 0, "up", "down")
                }
            }
        }
    }

    # Build alluvial-style data for ggplot
    if(length(timepoints) >= 2){
        alluvial_data <- data.frame()
        for(g_idx in seq(nrow(gene_status))){
            for(t_idx in seq_along(timepoints)){
                alluvial_data <- rbind(alluvial_data, data.frame(
                    gene = rownames(gene_status)[g_idx],
                    timepoint = timepoints[t_idx],
                    status = gene_status[g_idx, t_idx],
                    stringsAsFactors = FALSE))
            }
        }

        alluvial_data$timepoint <- factor(alluvial_data$timepoint, levels=timepoints)
        alluvial_data$status <- factor(alluvial_data$status, levels=categories)

        # Count transitions
        summary_data <- data.frame()
        for(t_idx in seq_along(timepoints)){
            tp <- timepoints[t_idx]
            counts <- table(gene_status[, t_idx])
            for(cat in names(counts)){
                summary_data <- rbind(summary_data, data.frame(
                    timepoint=tp, status=cat, count=as.numeric(counts[cat])))
            }
        }
        summary_data$timepoint <- factor(summary_data$timepoint, levels=timepoints)
        summary_data$status <- factor(summary_data$status, levels=categories)

        # Stacked bar chart showing DEG status over time
        status_colors <- c("up"="red", "down"="blue", "ns"="grey70")

        g <- ggplot(summary_data, aes(x=timepoint, y=count, fill=status)) +
            geom_bar(stat="identity", position="stack") +
            scale_fill_manual(values=status_colors) +
            theme_minimal() +
            theme(text = element_text(size=16)) +
            ylab(paste0("Top ", topK, " genes (Component ", comp, ")")) +
            xlab("Timepoint") +
            ggtitle(paste0(method, " - At Component ", comp, " loading genes"))

        pngfile <- paste0(outdir, "/at_comp", comp, ".png")
        ggsave(file=pngfile, plot=g, width=8, height=5)
    }
}

file.create(outfile)
