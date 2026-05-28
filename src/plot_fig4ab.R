source("src/Functions.R")

# Fig4A: parasitism1 (guided_pls, rbh_pls) + legend
# Fig4B: grafting (guided_pls, rbh_pls) + legend

panels <- list(
    list(panel = "Fig4A", sample = "parasitism1"),
    list(panel = "Fig4B", sample = "grafting")
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

        x <- 1
        y <- 2
        xlim <- range(score_At[, x], score_Pj[, x])
        ylim <- range(score_At[, y], score_Pj[, y])

        # Scatter plot (no arrows, no legend)
        for (ext in c("png", "pdf")) {
            outfile <- paste0(outdir, "/", panel, "_", method, ".", ext)
            if (ext == "png") {
                png(file = outfile, width = 700, height = 700)
            } else {
                pdf(file = outfile, width = 7, height = 7)
            }
            par(mar = c(4, 4, 2, 2))

            plot(score_At[, c(x, y)],
                col = label_At, xlim = xlim, ylim = ylim,
                xlab = paste0("Component-", x),
                ylab = paste0("Component-", y),
                pch = 20, cex = 3.5, bty = "n",
                main = "")
            par(new = TRUE)
            plot(score_Pj[, c(x, y)],
                col = label_Pj, xlim = xlim, ylim = ylim,
                pch = 18, cex = 3.5, bty = "n", ann = FALSE,
                xaxt = "n", yaxt = "n")
            segments(xlim[1], 0, xlim[2], 0, lty = 2)
            segments(0, ylim[1], 0, ylim[2], lty = 2)

            dev.off()
        }
    }

    # Legend as separate file
    labels <- c(label_At, label_Pj)
    pchs <- c(rep(20, length = length(unique(label_At))),
        rep(18, length = length(unique(label_Pj))))

    for (ext in c("png", "pdf")) {
        outfile <- paste0(outdir, "/", panel, "_legend.", ext)
        if (ext == "png") {
            png(file = outfile, width = 400, height = 700)
        } else {
            pdf(file = outfile, width = 4, height = 7)
        }
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend("center",
            col = labels[unique(names(labels))],
            legend = unique(names(labels)),
            pch = pchs, cex = 1.5, bg = "transparent",
            bty = "n")
        dev.off()
    }
}
