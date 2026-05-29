source("src/Functions.R")

# Fig4A: parasitism1 (guided_pls, rbh_pls) + legend
# Fig4B: grafting (guided_pls, rbh_pls) + legend

panels <- list(
    list(panel = "Fig4A", sample = "parasitism1", pt_size = 9),
    list(panel = "Fig4B", sample = "grafting", pt_size = 7)
)
methods <- c("guided_pls", "rbh_pls")

outdir <- "plot/Figures/main"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (p in panels) {
    sample <- p$sample
    panel <- p$panel
    load(paste0("data/", sample, ".RData"))

    for (method in methods) {
        infile <- paste0("output/", sample, "/", method, ".RData")
        load(infile)

        if (ncol(score_At) < 2) next

        # Build data frames
        df_At <- data.frame(
            x = score_At[, 1], y = score_At[, 2],
            condition = names(label_At),
            color = label_At,
            species = "At",
            stringsAsFactors = FALSE
        )
        df_Pj <- data.frame(
            x = score_Pj[, 1], y = score_Pj[, 2],
            condition = names(label_Pj),
            color = label_Pj,
            species = "Pj",
            stringsAsFactors = FALSE
        )
        df <- rbind(df_At, df_Pj)

        # Color mapping: condition -> color
        color_map <- unique(df[, c("condition", "color")])
        color_vec <- setNames(color_map$color, color_map$condition)

        g <- ggplot(df, aes(x = x, y = y, color = condition, shape = species)) +
            geom_point(size = p$pt_size) +
            scale_color_manual(values = color_vec) +
            scale_shape_manual(values = c("At" = 16, "Pj" = 18)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
            xlab("Component-1") + ylab("Component-2") +
            theme_minimal() +
            theme(
                text = element_text(size = 16),
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 18),
                legend.position = "none",
                panel.grid.minor = element_blank()
            )

        for (ext in c("png", "pdf")) {
            outfile <- paste0(outdir, "/", panel, "_", method, ".", ext)
            ggsave(file = outfile, plot = g, width = 7, height = 7, dpi = 300)
        }
    }

    # Legend as separate file
    # Combine all conditions with their shapes
    df_legend <- data.frame(
        condition = c(unique(names(label_At)), unique(names(label_Pj))),
        color = c(label_At[!duplicated(names(label_At))],
                  label_Pj[!duplicated(names(label_Pj))]),
        species = c(rep("At", length(unique(names(label_At)))),
                     rep("Pj", length(unique(names(label_Pj))))),
        stringsAsFactors = FALSE
    )
    color_vec_leg <- setNames(df_legend$color, df_legend$condition)

    g_legend <- ggplot(df_legend, aes(x = 1, y = 1,
            color = condition, shape = species)) +
        geom_point(size = 5) +
        scale_color_manual(values = color_vec_leg, name = "") +
        scale_shape_manual(values = c("At" = 16, "Pj" = 18), name = "") +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)),
               shape = guide_legend(ncol = 1, override.aes = list(size = 4))) +
        theme_minimal() +
        theme(legend.text = element_text(size = 16),
              legend.key.size = unit(1, "cm"))

    legend_grob <- cowplot::get_legend(
        g_legend + theme(legend.position = "right"))

    for (ext in c("png", "pdf")) {
        outfile <- paste0(outdir, "/", panel, "_legend.", ext)
        ggsave(file = outfile, plot = legend_grob,
            width = 5, height = 10, dpi = 300)
    }
}
