source("src/Functions.R")
library("viridisLite")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
method <- args[2]
outfile <- args[3]

load(paste0("data/", sample, ".RData"))

####################################
# Load DEG results
####################################
.load_deg_table <- function(sample, species, conditions, counts_mat){
    qtable <- c()
    for(cond in conditions){
        filename <- paste0("output/", sample, "/deg/", species, "_", cond, ".RData")
        if(!file.exists(filename)) next
        load(filename)
        pos_vec <- rep(0, length=nrow(counts_mat))
        neg_vec <- rep(0, length=nrow(counts_mat))
        names(pos_vec) <- rownames(counts_mat)
        names(neg_vec) <- rownames(counts_mat)
        tmp <- deg@.Data[[1]]
        pos <- head(rownames(tmp)[which(tmp$logFC > 0)], 0.1*nrow(counts_mat))
        neg <- head(rownames(tmp)[which(tmp$logFC < 0)], 0.1*nrow(counts_mat))
        pos_vec[pos] <- 1
        neg_vec[neg] <- 1
        qtable <- cbind(qtable, pos_vec, neg_vec)
    }
    cnames <- unlist(lapply(conditions, function(x) paste0(x, c("+", "-"))))
    colnames(qtable) <- cnames
    qtable
}

qtable_at <- .load_deg_table(sample, "at", deg_conditions, counts_At)
qtable_pj <- .load_deg_table(sample, "pj", deg_conditions, counts_Pj)

# Add procambium markers
pos_at <- rep(0, nrow(qtable_at)); names(pos_at) <- rownames(counts_At)
neg_at <- rep(0, nrow(qtable_at)); names(neg_at) <- rownames(counts_At)
pos <- intersect(procambium_At$GENEID[which(procambium_At$avg_logFC > 0)], rownames(qtable_at))
neg <- intersect(procambium_At$GENEID[which(procambium_At$avg_logFC < 0)], rownames(qtable_at))
pos_at[pos] <- 1; neg_at[neg] <- 1
qtable_at <- cbind(qtable_at, "procambium+"=pos_at, "procambium-"=neg_at)

pos_pj <- rep(0, nrow(qtable_pj)); names(pos_pj) <- rownames(counts_Pj)
neg_pj <- rep(0, nrow(qtable_pj)); names(neg_pj) <- rownames(counts_Pj)
pos <- intersect(procambium_Pj$GENEID[which(procambium_Pj$avg_logFC > 0)], rownames(qtable_pj))
neg <- intersect(procambium_Pj$GENEID[which(procambium_Pj$avg_logFC < 0)], rownames(qtable_pj))
pos_pj[pos] <- 1; neg_pj[neg] <- 1
qtable_pj <- cbind(qtable_pj, "procambium+"=pos_pj, "procambium-"=neg_pj)

####################################
# Load method loadings
####################################
infile <- paste0("output/", sample, "/", method, ".RData")
load(infile)

if(method == "guided_pls"){
    loadings_at <- resgPLS$loading1
    loadings_pj <- resgPLS$loading2
    rownames(loadings_at) <- rownames(X_At)
    rownames(loadings_pj) <- rownames(X_Pj)
}
if(method == "guided_pca"){
    loadings_at <- resgPCA_At$u
    loadings_pj <- resgPCA_Pj$u
    rownames(loadings_at) <- rownames(X_At)
    rownames(loadings_pj) <- rownames(X_Pj)
}
if(method == "pca"){
    loadings_at <- res_pca_At$u[, seq(4)]
    loadings_pj <- res_pca_Pj$u[, seq(4)]
    rownames(loadings_at) <- rownames(X_At)
    rownames(loadings_pj) <- rownames(X_Pj)
}
if(method == "rbh_pls"){
    loadings_at <- scaled_X_At_RBH %*% resgPLS$u[, seq(4)]
    loadings_pj <- scaled_X_Pj_RBH %*% resgPLS$v[, seq(4)]
    rownames(loadings_at) <- rownames(scaled_X_At_RBH)
    rownames(loadings_pj) <- rownames(scaled_X_Pj_RBH)
    # Align gene sets
    common_at <- intersect(rownames(qtable_at), rownames(loadings_at))
    loadings_at <- loadings_at[common_at, ]
    qtable_at <- qtable_at[common_at, ]
    common_pj <- intersect(rownames(qtable_pj), rownames(loadings_pj))
    loadings_pj <- loadings_pj[common_pj, ]
    qtable_pj <- qtable_pj[common_pj, ]
}

####################################
# Positive/Negative split
####################################
loadings_at <- loadings_at[, unlist(lapply(seq(ncol(loadings_at)), function(x) rep(x, 2)))]
loadings_pj <- loadings_pj[, unlist(lapply(seq(ncol(loadings_pj)), function(x) rep(x, 2)))]
for(i in seq(ncol(loadings_at))){
    if(i %% 2 == 1){
        loadings_at[,i] <- -loadings_at[,i]
        loadings_pj[,i] <- -loadings_pj[,i]
    }
}
colnames(loadings_at) <- unlist(lapply(seq(ncol(loadings_at)/2),
    function(x) paste0("Dim", x, c("+", "-"))))
colnames(loadings_pj) <- unlist(lapply(seq(ncol(loadings_pj)/2),
    function(x) paste0("Dim", x, c("+", "-"))))

####################################
# Binarization
####################################
.binarize_loading <- function(loadings_mat){
    apply(loadings_mat, 2, function(x){
        out <- x
        out[which(rank(x) <= 0.1*length(x))] <- 1
        out[which(rank(x) > 0.1*length(x))] <- 0
        out
    })
}
loadings_at <- .binarize_loading(loadings_at)
loadings_pj <- .binarize_loading(loadings_pj)

####################################
# Jaccard index
####################################
numer_at <- t(qtable_at) %*% loadings_at
denom_at <- outer(colSums(qtable_at), colSums(loadings_at), "+") - numer_at
jaccard_at <- numer_at / denom_at

numer_pj <- t(qtable_pj) %*% loadings_pj
denom_pj <- outer(colSums(qtable_pj), colSums(loadings_pj), "+") - numer_pj
jaccard_pj <- numer_pj / denom_pj

####################################
# Plot
####################################
gdata_at <- melt(jaccard_at)
gdata_pj <- melt(jaccard_pj)
colnames(gdata_at)[3] <- "jaccard"
colnames(gdata_pj)[3] <- "jaccard"
gdata <- rbind(gdata_at, gdata_pj)
gdata$species <- c(rep("At", nrow(gdata_at)), rep("Pj", nrow(gdata_pj)))

g <- ggplot(gdata, aes(x=Var2, y=Var1, fill = jaccard))
g <- g + geom_tile()
g <- g + facet_wrap(~ species)
g <- g + scale_fill_gradientn(colours = viridis(100))
g <- g + theme(text = element_text(size=24))
g <- g + theme(axis.title.x = element_blank())
g <- g + theme(axis.title.y = element_blank())
g <- g + theme(axis.text.x = element_text(size=18, angle=90, hjust=1))
g <- g + ggtitle(method)

ggsave(file=outfile, plot=g, width=15, height=12)
