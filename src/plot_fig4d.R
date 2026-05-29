source("src/Functions.R")

outdir <- "plot/Figures/main"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

samples <- c("parasitism1", "grafting")

# Collect data
spec_list <- list()
conc_list <- list()

for (sample in samples) {
    s <- read.csv(paste0("output/", sample, "/evaluation/factor_specificity_summary.csv"))
    s$sample <- sample
    spec_list[[sample]] <- s

    c <- read.csv(paste0("output/", sample, "/evaluation/factor_concordance.csv"))
    c$sample <- sample
    conc_list[[sample]] <- c
}

spec_df <- do.call(rbind, spec_list)
conc_df <- do.call(rbind, conc_list)

# Method and sample labels
method_levels <- c("guided_pls", "rbh_pls", "guided_pca", "pca")
method_labels <- c("guided-PLS", "RBH-PLS", "guided-PCA", "PCA")

spec_df$method <- factor(spec_df$method, levels = method_levels, labels = method_labels)
conc_df$method <- factor(conc_df$method, levels = method_levels, labels = method_labels)

# Common theme (no x-axis text, Dark2 palette)
common_theme <- theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(size = 16),
          legend.position = "none")

# Factor specificity — average over species
spec_avg <- aggregate(mean_specificity ~ method + sample, data = spec_df, FUN = mean)

# Save individual plots
for (sample in samples) {
    # Factor specificity
    d_spec <- spec_avg[spec_avg$sample == sample, ]
    g <- ggplot(d_spec, aes(x = method, y = mean_specificity, fill = method)) +
        geom_bar(stat = "identity", width = 0.7) +
        common_theme +
        scale_fill_brewer(palette = "Dark2") +
        ylab("Factor specificity") +
        xlab("") +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.4)

    ggsave(file = paste0(outdir, "/Fig4D_factor_specificity_", sample, ".png"),
        plot = g, width = 5, height = 4, dpi = 300)
    ggsave(file = paste0(outdir, "/Fig4D_factor_specificity_", sample, ".pdf"),
        plot = g, width = 5, height = 4)

    # Cross-species concordance
    d_conc <- conc_df[conc_df$sample == sample, ]
    d_conc$note <- ""
    if (sample == "grafting") {
        d_conc$note[d_conc$method == "RBH-PLS"] <- "*"
    }

    g <- ggplot(d_conc, aes(x = method, y = concordance_rate, fill = method)) +
        geom_bar(stat = "identity", width = 0.7) +
        geom_text(aes(label = note), vjust = -0.3, size = 5) +
        common_theme +
        scale_fill_brewer(palette = "Dark2") +
        ylab("Cross-species concordance") +
        xlab("") +
        ylim(0, 1.1)

    ggsave(file = paste0(outdir, "/Fig4D_concordance_", sample, ".png"),
        plot = g, width = 5, height = 4, dpi = 300)
    ggsave(file = paste0(outdir, "/Fig4D_concordance_", sample, ".pdf"),
        plot = g, width = 5, height = 4)
}

# Legend as separate horizontal file
g_legend <- ggplot(spec_avg, aes(x = method, y = mean_specificity, fill = method)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Dark2", name = "") +
    theme_minimal() +
    theme(legend.direction = "horizontal",
          legend.text = element_text(size = 20),
          legend.key.size = unit(1.2, "cm"))

legend_grob <- cowplot::get_legend(g_legend + theme(legend.position = "bottom"))

ggsave(file = paste0(outdir, "/Fig4D_legend.png"),
    plot = legend_grob, width = 10, height = 1.5, dpi = 300)
ggsave(file = paste0(outdir, "/Fig4D_legend.pdf"),
    plot = legend_grob, width = 10, height = 1.5)
