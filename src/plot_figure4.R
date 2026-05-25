source("src/Functions.R")

outdir <- "plot/Figure/3"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

####################################
# Figure 3: Quantitative comparison across methods
# (a) Factor specificity score
# (b) Cross-species concordance rate
####################################
samples <- c("parasitism1", "grafting")

# Collect specificity summaries
spec_list <- list()
conc_list <- list()

for(sample in samples){
    s <- read.csv(paste0("output/", sample, "/evaluation/factor_specificity_summary.csv"))
    s$sample <- sample
    spec_list[[sample]] <- s

    c <- read.csv(paste0("output/", sample, "/evaluation/factor_concordance.csv"))
    c$sample <- sample
    conc_list[[sample]] <- c
}

spec_df <- do.call(rbind, spec_list)
conc_df <- do.call(rbind, conc_list)

# Method labels
method_levels <- c("guided_pls", "rbh_pls", "guided_pca", "pca")
method_labels <- c("guided-PLS", "RBH-PLS", "guided-PCA", "PCA")

spec_df$method <- factor(spec_df$method, levels=method_levels, labels=method_labels)
conc_df$method <- factor(conc_df$method, levels=method_levels, labels=method_labels)

# (a) Factor specificity — average over species
spec_avg <- aggregate(mean_specificity ~ method + sample, data=spec_df, FUN=mean)

g_spec <- ggplot(spec_avg, aes(x=method, y=mean_specificity, fill=method)) +
    geom_bar(stat="identity", width=0.7) +
    facet_wrap(~ sample) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13),
          legend.position="none",
          strip.text = element_text(size=12, face="bold")) +
    scale_fill_brewer(palette="Set2") +
    ylab("Factor specificity") +
    xlab("") +
    geom_hline(yintercept=0.5, linetype="dashed", color="grey50", linewidth=0.4)

# (b) Cross-species concordance
# For RBH-PLS grafting, all dims map to 0d → concordance=1.0 is misleading
# Add annotation
conc_df$note <- ""
conc_df$note[conc_df$method == "RBH-PLS" & conc_df$sample == "grafting"] <- "*"

g_conc <- ggplot(conc_df, aes(x=method, y=concordance_rate, fill=method)) +
    geom_bar(stat="identity", width=0.7) +
    geom_text(aes(label=note), vjust=-0.3, size=5) +
    facet_wrap(~ sample) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=30, hjust=1, size=11),
          text = element_text(size=13),
          legend.position="none",
          strip.text = element_text(size=12, face="bold")) +
    scale_fill_brewer(palette="Set2") +
    ylab("Cross-species concordance") +
    xlab("") +
    ylim(0, 1.1)

fig4 <- g_spec / g_conc +
    plot_annotation(
        title = "Factor specificity and cross-species concordance",
        tag_levels = "a",
        theme = theme(plot.title = element_text(size=14, face="bold")))

ggsave(file=paste0(outdir, "/Figure3.png"), plot=fig4, width=10, height=9, dpi=300)
ggsave(file=paste0(outdir, "/Figure3.pdf"), plot=fig4, width=10, height=9)
