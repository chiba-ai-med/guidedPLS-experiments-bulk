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

        g <- ggplot(d, aes(x = condition, y = dim_label, fill = neg_log10_p)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = round(neg_log10_p, 1)), size = 6.4, color = "white") +
            scale_fill_gradientn(colours = viridis(100), name = "-log10(p)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                  axis.text.y = element_text(size = 14),
                  legend.key.height = unit(0.4, "cm"),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 14)) +
            xlab("") + ylab("")

        fname <- paste0("Fig4C_heatmap_", sample, "_", sp)
        ggsave(file = paste0(outdir, "/", fname, ".png"),
            plot = g, width = 7, height = 4, dpi = 300)
        ggsave(file = paste0(outdir, "/", fname, ".pdf"),
            plot = g, width = 7, height = 4)
    }
}
