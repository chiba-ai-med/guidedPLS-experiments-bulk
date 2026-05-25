source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
method <- args[2]
outfile <- args[3]

load(paste0("data/", sample, ".RData"))
infile <- paste0("output/", sample, "/", method, ".RData")
load(infile)

# Score correlation vectors for biplot arrows
score_cor1 <- NULL
score_cor2 <- NULL
scale1 <- 1
scale2 <- 1

if(method == "guided_pls"){
    score_cor1 <- resgPLS$score_cor1
    score_cor2 <- resgPLS$score_cor2
    scale1 <- 1/7
    scale2 <- 1/7
}
if(method == "guided_pca"){
    scaled_Y_At_all <- scale(Y_At_all, center=TRUE, scale=TRUE)
    scaled_Y_Pj_all <- scale(Y_Pj_all, center=TRUE, scale=TRUE)
    score_cor1 <- t(resgPCA_At$u) %*% scaled_X_At %*% scaled_Y_At_all
    score_cor2 <- t(resgPCA_Pj$u) %*% scaled_X_Pj %*% scaled_Y_Pj_all
    scale1 <- 1/15
    scale2 <- 1/15
}

# Output directory
outdir <- paste0("plot/", sample, "/scatter/", method)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if(ncol(score_At) >= 2){
    dimcombn <- combn(ncol(score_At), 2)
    labels <- c(label_At, label_Pj)
    pchs <- c(rep(20, length=length(unique(label_At))),
        rep(18, length=length(unique(label_Pj))))

    for(i in seq(ncol(dimcombn))){
        x <- dimcombn[1,i]
        y <- dimcombn[2,i]
        xlim <- c(min(score_At[,x], score_Pj[,x]),
            max(score_At[,x], score_Pj[,x]))
        ylim <- c(min(score_At[,y], score_Pj[,y]),
            max(score_At[,y], score_Pj[,y]))

        outfile_plot <- paste0(outdir, "/", x, "_", y, ".png")
        png(file=outfile_plot, width=900, height=700)
        par(mar = c(3,3,3,12))

        if(!is.null(score_cor1)){
            # Empty plot
            plot(score_At[,c(x,y)],
                col=rgb(1,1,1), xlim=xlim, ylim=ylim,
                xlab="", ylab="",
                pch=20, cex=3.5, bty="n")
            # Arrows (At: solid, Pj: dashed)
            arrows(0,0,
                (score_cor1[,c(x,y)]*scale1)[,1],
                (score_cor1[,c(x,y)]*scale1)[,2],
                col=colorAt2[seq(nrow(score_cor1))], lty=1)
            arrows(0,0,
                (score_cor2[,c(x,y)]*scale2)[,1],
                (score_cor2[,c(x,y)]*scale2)[,2],
                col=colorPj2[seq(nrow(score_cor2))], lty=2)
            text((score_cor1[,c(x,y)]*scale1)[,1],
                (score_cor1[,c(x,y)]*scale1)[,2],
                label=paste0(rownames(score_cor1), "(At)"),
                col=colorAt2[seq(nrow(score_cor1))], cex=1.5)
            text((score_cor2[,c(x,y)]*scale2)[,1],
                (score_cor2[,c(x,y)]*scale2)[,2],
                label=paste0(rownames(score_cor2), "(Pj)"),
                col=colorPj2[seq(nrow(score_cor2))], cex=1.5)
            par(new=TRUE)
        }

        # Points
        plot(score_At[,c(x,y)],
            col=label_At, xlim=xlim, ylim=ylim,
            xlab=paste0("Component-", x), ylab=paste0("Component-", y),
            pch=20, cex=3.5, bty="n")
        par(new=TRUE)
        plot(score_Pj[,c(x,y)],
            col=label_Pj, xlim=xlim, ylim=ylim,
            pch=18, cex=3.5, bty="n", ann=FALSE)
        par(xpd = TRUE)
        legend(par()$usr[2], par()$usr[4],
            col = labels[unique(names(labels))],
            legend = unique(names(labels)),
            pch = pchs, cex = 1.5, bg = "transparent")
        segments(xlim[1], 0, 0, 0, lty=2)
        segments(xlim[2], 0, 0, 0, lty=2)
        segments(0, ylim[1], 0, 0, lty=2)
        segments(0, ylim[2], 0, 0, lty=2)
        dev.off()
    }
}else{
    # 1D case
    xlim <- c(min(score_At[,1], score_Pj[,1]),
        max(score_At[,1], score_Pj[,1]))
    labels <- c(label_At, label_Pj)

    outfile_plot <- paste0(outdir, "/1_1.png")
    png(file=outfile_plot, width=1200, height=300)
    par(mar = c(3,3,3,12))
    plot(cbind(score_At[,1], 1),
        xlab="Component-1", xlim=xlim, ylim=c(0.5, 2.5),
        pch=20, cex=2, bty="n", col=label_At)
    par(new=TRUE)
    plot(cbind(score_Pj[,1], 2),
        xlab="", xlim=xlim, ylim=c(0.5, 2.5),
        pch=18, cex=2, bty="n", col=label_Pj, ann=FALSE)
    dev.off()
}

# Output
file.create(outfile)
