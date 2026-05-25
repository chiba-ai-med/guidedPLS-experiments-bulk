source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

load(paste0("data/", sample, ".RData"))

methods <- c("guided_pls", "rbh_pls")
# guided-PCA and PCA are single-species, so alignment is not applicable

####################################
# Define shared condition labels for At and Pj
####################################
if(sample == "parasitism1"){
    # Shared conditions: timepoint × genotype for parasitized samples
    # At: Pj.Col.1d = wt+Pj 1d, Pj.wol.1d = wol+Pj 1d, etc.
    # Pj: Pj.Col.1d = Pj parasitizing wt 1d, Pj.wol.1d = Pj parasitizing wol 1d
    shared_conditions <- c("1d_wt", "3d_wt", "7d_wt", "1d_wol", "3d_wol", "7d_wol")

    # Assign condition labels to At samples
    cond_At <- rep(NA, ncol(X_At))
    names(cond_At) <- colnames(X_At)
    cond_At[grep("Pj\\.Col\\.1d", names(cond_At))] <- "1d_wt"
    cond_At[grep("Pj\\.Col\\.3d", names(cond_At))] <- "3d_wt"
    cond_At[grep("Pj\\.Col\\.7d", names(cond_At))] <- "7d_wt"
    cond_At[grep("Pj\\.wol\\.1d", names(cond_At))] <- "1d_wol"
    cond_At[grep("Pj\\.wol\\.3d", names(cond_At))] <- "3d_wol"
    cond_At[grep("Pj\\.wol\\.7d", names(cond_At))] <- "7d_wol"

    # Assign condition labels to Pj samples
    cond_Pj <- rep(NA, ncol(X_Pj))
    names(cond_Pj) <- colnames(X_Pj)
    cond_Pj[grep("Pj\\.Col\\.1d", names(cond_Pj))] <- "1d_wt"
    cond_Pj[grep("Pj\\.Col\\.3d", names(cond_Pj))] <- "3d_wt"
    cond_Pj[grep("Pj\\.Col\\.7d", names(cond_Pj))] <- "7d_wt"
    cond_Pj[grep("Pj\\.wol\\.1d", names(cond_Pj))] <- "1d_wol"
    cond_Pj[grep("Pj\\.wol\\.3d", names(cond_Pj))] <- "3d_wol"
    cond_Pj[grep("Pj\\.wol\\.7d", names(cond_Pj))] <- "7d_wol"
}

if(sample == "grafting"){
    # Shared: both organisms grafted together at each timepoint
    shared_conditions <- c("1d", "2d", "3d", "4d", "5d", "6d", "7d", "14d")

    cond_At <- rep(NA, ncol(X_At))
    names(cond_At) <- colnames(X_At)
    cond_At[grep("Gra01bo", names(cond_At))] <- "1d"
    cond_At[grep("Gra02bo", names(cond_At))] <- "2d"
    cond_At[grep("Gra03bo", names(cond_At))] <- "3d"
    cond_At[grep("Gra04bo", names(cond_At))] <- "4d"
    cond_At[grep("Gra05bo", names(cond_At))] <- "5d"
    cond_At[grep("Gra06bo", names(cond_At))] <- "6d"
    cond_At[grep("Gra07bo", names(cond_At))] <- "7d"
    cond_At[grep("Gra14bo", names(cond_At))] <- "14d"

    cond_Pj <- rep(NA, ncol(X_Pj))
    names(cond_Pj) <- colnames(X_Pj)
    cond_Pj[grep("Gra01bo", names(cond_Pj))] <- "1d"
    cond_Pj[grep("Gra02bo", names(cond_Pj))] <- "2d"
    cond_Pj[grep("Gra03bo", names(cond_Pj))] <- "3d"
    cond_Pj[grep("Gra04bo", names(cond_Pj))] <- "4d"
    cond_Pj[grep("Gra05bo", names(cond_Pj))] <- "5d"
    cond_Pj[grep("Gra06bo", names(cond_Pj))] <- "6d"
    cond_Pj[grep("Gra07bo", names(cond_Pj))] <- "7d"
    cond_Pj[grep("Gra14bo", names(cond_Pj))] <- "14d"
}

####################################
# Evaluate alignment for each method
####################################
results <- list()

for(m in methods){
    infile <- paste0("output/", sample, "/", m, ".RData")
    if(!file.exists(infile)) next
    load(infile)

    # Scores are samples x components
    s_At <- score_At
    s_Pj <- score_Pj

    # Filter to shared-condition samples only
    at_shared <- which(!is.na(cond_At))
    pj_shared <- which(!is.na(cond_Pj))

    s_At_shared <- s_At[at_shared, , drop=FALSE]
    s_Pj_shared <- s_Pj[pj_shared, , drop=FALSE]
    c_At_shared <- cond_At[at_shared]
    c_Pj_shared <- cond_Pj[pj_shared]

    # Metric 1: Mean cross-species distance for same-condition pairs
    # For each condition, compute centroid distance between At and Pj
    same_cond_dist <- c()
    diff_cond_dist <- c()

    centroids_At <- list()
    centroids_Pj <- list()
    for(cond in shared_conditions){
        idx_at <- which(c_At_shared == cond)
        idx_pj <- which(c_Pj_shared == cond)
        if(length(idx_at) > 0 && length(idx_pj) > 0){
            cent_at <- colMeans(s_At_shared[idx_at, , drop=FALSE])
            cent_pj <- colMeans(s_Pj_shared[idx_pj, , drop=FALSE])
            centroids_At[[cond]] <- cent_at
            centroids_Pj[[cond]] <- cent_pj
            same_cond_dist <- c(same_cond_dist, sqrt(sum((cent_at - cent_pj)^2)))
        }
    }

    # Cross-condition distances (different conditions)
    conds_with_data <- names(centroids_At)
    if(length(conds_with_data) >= 2){
        for(i in seq(length(conds_with_data)-1)){
            for(j in (i+1):length(conds_with_data)){
                c1 <- conds_with_data[i]
                c2 <- conds_with_data[j]
                d_at_pj <- sqrt(sum((centroids_At[[c1]] - centroids_Pj[[c2]])^2))
                diff_cond_dist <- c(diff_cond_dist, d_at_pj)
                d_pj_at <- sqrt(sum((centroids_Pj[[c1]] - centroids_At[[c2]])^2))
                diff_cond_dist <- c(diff_cond_dist, d_pj_at)
            }
        }
    }

    mean_same <- mean(same_cond_dist)
    mean_diff <- mean(diff_cond_dist)
    alignment_ratio <- mean_same / mean_diff  # Lower = better alignment

    # Metric 2: Silhouette-like score
    # For each shared sample, is it closer to its same-condition cross-species centroid
    # than to different-condition cross-species centroids?
    correct <- 0
    total <- 0
    for(cond in conds_with_data){
        idx_at <- which(c_At_shared == cond)
        for(i in idx_at){
            d_same <- sqrt(sum((s_At_shared[i,] - centroids_Pj[[cond]])^2))
            d_others <- sapply(setdiff(conds_with_data, cond), function(c2){
                sqrt(sum((s_At_shared[i,] - centroids_Pj[[c2]])^2))
            })
            if(d_same < mean(d_others)) correct <- correct + 1
            total <- total + 1
        }
        idx_pj <- which(c_Pj_shared == cond)
        for(i in idx_pj){
            d_same <- sqrt(sum((s_Pj_shared[i,] - centroids_At[[cond]])^2))
            d_others <- sapply(setdiff(conds_with_data, cond), function(c2){
                sqrt(sum((s_Pj_shared[i,] - centroids_At[[c2]])^2))
            })
            if(d_same < mean(d_others)) correct <- correct + 1
            total <- total + 1
        }
    }
    alignment_accuracy <- correct / total

    results[[length(results)+1]] <- data.frame(
        method=m,
        mean_same_cond_dist=mean_same,
        mean_diff_cond_dist=mean_diff,
        alignment_ratio=alignment_ratio,
        alignment_accuracy=alignment_accuracy)
}

result_df <- do.call(rbind, results)

# Save
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
write.csv(result_df, outfile, row.names=FALSE)

# Plot
result_df$method <- factor(result_df$method, levels=c("guided_pls", "rbh_pls"))

g1 <- ggplot(result_df, aes(x=method, y=alignment_accuracy, fill=method)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(text = element_text(size=14)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Cross-species alignment accuracy") +
    xlab("Method") +
    ylim(0, 1) +
    ggtitle(paste0("Species alignment: ", sample))

outdir <- paste0("plot/", sample)
ggsave(file=paste0(outdir, "/alignment_accuracy.png"), plot=g1, width=6, height=5)

g2 <- ggplot(result_df, aes(x=method, y=alignment_ratio, fill=method)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(text = element_text(size=14)) +
    scale_fill_brewer(palette="Set2") +
    ylab("Same/Different condition distance ratio\n(lower = better)") +
    xlab("Method") +
    ggtitle(paste0("Distance ratio: ", sample))

ggsave(file=paste0(outdir, "/alignment_ratio.png"), plot=g2, width=6, height=5)
