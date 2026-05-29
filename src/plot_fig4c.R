source("src/Functions.R")
library("viridisLite")

outdir <- "plot/Figures/main"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

samples <- c("parasitism1", "grafting")

for (sample in samples) {
    df <- read.csv(paste0("output/", sample, "/evaluation/factor_specificity.csv"))
    sub <- df[df$method == "guided_pls", ]

    # Rename "parasm" to "parasitism" in condition labels (before factor conversion)
    sub$condition <- as.character(sub$condition)
    sub$condition[sub$condition == "parasm"] <- "parasitism"

    for (sp in c("At", "Pj")) {
        d <- sub[sub$species == sp, ]
        d$dim_label <- paste0("Dim", d$dim)
        d$dim_label <- factor(d$dim_label, levels = paste0("Dim", seq(max(d$dim))))
        d$condition <- factor(d$condition, levels = unique(d$condition))

        n_cond <- length(unique(d$condition))
        # Scale width proportionally so tile size is consistent across datasets
        fig_width <- 1.4 * n_cond + 2
        fig_height <- 6

        g <- ggplot(d, aes(x = condition, y = dim_label, fill = neg_log10_p)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = round(neg_log10_p, 1)), size = 9, color = "white") +
            scale_fill_gradientn(colours = viridis(100), name = "-log10(p)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
                  axis.text.y = element_text(size = 28),
                  legend.position = "none") +
            xlab("") + ylab("")

        fname <- paste0("Fig4C_heatmap_", sample, "_", sp)
        ggsave(file = paste0(outdir, "/", fname, ".png"),
            plot = g, width = fig_width, height = fig_height, dpi = 300)
        ggsave(file = paste0(outdir, "/", fname, ".pdf"),
            plot = g, width = fig_width, height = fig_height)
    }
}

# Legend as separate file
g_legend <- ggplot(data.frame(x = 1, y = 1, z = c(0, 5, 10)),
        aes(x = x, y = y, fill = z)) +
    geom_tile() +
    scale_fill_gradientn(colours = viridis(100), name = "-log10(p)") +
    theme_minimal() +
    theme(legend.direction = "horizontal",
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.key.width = unit(2, "cm"),
          legend.key.height = unit(0.6, "cm"))

legend_grob <- cowplot::get_legend(g_legend + theme(legend.position = "bottom"))

ggsave(file = paste0(outdir, "/Fig4C_legend.png"),
    plot = legend_grob, width = 8, height = 1.5, dpi = 300)
ggsave(file = paste0(outdir, "/Fig4C_legend.pdf"),
    plot = legend_grob, width = 8, height = 1.5)
