source("src/Functions.R")

# Score scatter plots for Fig4B (no arrows, no legend)
samples <- c("grafting", "parasitism1")
methods <- c("guided_pls", "rbh_pls", "guided_pca", "pca")

outdir <- "plot/Figures/main"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (sample in samples) {
    load(paste0("data/", sample, ".RData"))
    for (method in methods) {
        infile <- paste0("output/", sample, "/", method, ".RData")
        load(infile)

        if (ncol(score_At) < 2) next

        x <- 1
        y <- 2
        xlim <- range(score_At[, x], score_Pj[, x])
        ylim <- range(score_At[, y], score_Pj[, y])

        outfile <- paste0(outdir, "/Fig4B_", method, "_", sample, ".png")
        png(file = outfile, width = 700, height = 700)
        par(mar = c(4, 4, 2, 2))

        plot(score_At[, c(x, y)],
            col = label_At, xlim = xlim, ylim = ylim,
            xlab = paste0("Component-", x),
            ylab = paste0("Component-", y),
            pch = 20, cex = 3.5, bty = "n",
            main = paste0(method, " (", sample, ")"))
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
