source("src/Functions.R")
library("viridisLite")

outdir <- "plot/Figure/2"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

####################################
# Figure 2: Dimension x Condition enrichment heatmaps (guided-PLS only)
# (a) parasitism1 At, (b) parasitism1 Pj
# (c) grafting At, (d) grafting Pj
####################################
samples <- c("parasitism1", "grafting")
panels <- list()
panel_idx <- 0

for(sample in samples){
    df <- read.csv(paste0("output/", sample, "/evaluation/factor_specificity.csv"))
    sub <- df[df$method == "guided_pls", ]

    for(sp in c("At", "Pj")){
        panel_idx <- panel_idx + 1
        d <- sub[sub$species == sp, ]
        d$dim_label <- paste0("Dim", d$dim)
        d$dim_label <- factor(d$dim_label, levels=paste0("Dim", seq(max(d$dim))))
        d$condition <- factor(d$condition, levels=unique(d$condition))

        panels[[panel_idx]] <- ggplot(d, aes(x=condition, y=dim_label, fill=neg_log10_p)) +
            geom_tile(color="white", linewidth=0.5) +
            geom_text(aes(label=round(neg_log10_p, 1)), size=3.2, color="white") +
            scale_fill_gradientn(colours=viridis(100), name="-log10(p)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle=45, hjust=1, size=10),
                  axis.text.y = element_text(size=10),
                  plot.title = element_text(size=12, face="bold"),
                  legend.key.height = unit(0.4, "cm")) +
            xlab("") + ylab("") +
            ggtitle(paste0(sample, " — ", sp))
    }
}

fig3 <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
    plot_annotation(
        title = "guided-PLS: dimension–factor correspondence",
        tag_levels = "a",
        theme = theme(plot.title = element_text(size=14, face="bold")))

ggsave(file=paste0(outdir, "/Figure2.png"), plot=fig3, width=14, height=8, dpi=300)
ggsave(file=paste0(outdir, "/Figure2.pdf"), plot=fig3, width=14, height=8)
