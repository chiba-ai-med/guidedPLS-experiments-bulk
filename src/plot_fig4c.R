source("src/Functions.R")
library("viridisLite")

outdir <- "plot/Figures/main"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

samples <- c("parasitism1", "grafting")
panels <- list()
panel_idx <- 0

for (sample in samples) {
    df <- read.csv(paste0("output/", sample, "/evaluation/factor_specificity.csv"))
    sub <- df[df$method == "guided_pls", ]

    for (sp in c("At", "Pj")) {
        panel_idx <- panel_idx + 1
        d <- sub[sub$species == sp, ]
        d$dim_label <- paste0("Dim", d$dim)
        d$dim_label <- factor(d$dim_label, levels = paste0("Dim", seq(max(d$dim))))
        d$condition <- factor(d$condition, levels = unique(d$condition))

        panels[[panel_idx]] <- ggplot(d, aes(x = condition, y = dim_label, fill = neg_log10_p)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = round(neg_log10_p, 1)), size = 6.4, color = "white") +
            scale_fill_gradientn(colours = viridis(100), name = "-log10(p)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                  axis.text.y = element_text(size = 14),
                  legend.key.height = unit(0.4, "cm")) +
            xlab("") + ylab("") +
            ggtitle(paste0(sample, " \u2014 ", sp))
    }
}

fig4c <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
    plot_annotation(tag_levels = "a")

ggsave(file = paste0(outdir, "/Fig4C_dimension_factor_heatmap.png"),
    plot = fig4c, width = 14, height = 8, dpi = 300)
ggsave(file = paste0(outdir, "/Fig4C_dimension_factor_heatmap.pdf"),
    plot = fig4c, width = 14, height = 8)
