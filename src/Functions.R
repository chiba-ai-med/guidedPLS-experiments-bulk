library("guidedPLS")
library("edgeR")
library("RColorBrewer")
library("ggplot2")
library("tagcloud")
library("reshape2")
library("patchwork")
library("GO.db")
library("irlba")

####################################
# guidedPLS function
####################################
.guidedPLS <- function(X1, X2, Y1, Y2, k, cortest=FALSE){
    # Auto scaling
    YX1 <- scale(t(Y1) %*% X1, center=TRUE, scale=TRUE)
    YX2 <- scale(t(Y2) %*% X2, center=TRUE, scale=TRUE)
    # guided PLS
    res <- irlba(t(YX1) %*% YX2, k)
    # Loading
    loading1 <- as.matrix(res$u[,seq(k)])
    loading2 <- as.matrix(res$v[,seq(k)])
    # Score
    score1 <- X1 %*% loading1
    score2 <- X2 %*% loading2
    # Smaller Score
    score_cor1 <- YX1 %*% loading1
    score_cor2 <- YX2 %*% loading2
    if(cortest){
        # Correlation Coefficient
        cor1 <- matrix(0, nrow=ncol(YX1), ncol=k)
        cor2 <- matrix(0, nrow=ncol(YX2), ncol=k)
        for(i in seq(k)){
            cor1[,i] <- apply(YX1, 2, function(x){
                cor(x, score_cor2[, i])
            })
            cor2[,i] <- apply(YX2, 2, function(x){
                cor(x, score_cor1[, i])
            })
        }
        # P-value / Q-value
        pval1 <- matrix(0, nrow=ncol(YX1), ncol=k)
        pval2 <- matrix(0, nrow=ncol(YX2), ncol=k)
        qval1 <- matrix(0, nrow=ncol(YX1), ncol=k)
        qval2 <- matrix(0, nrow=ncol(YX2), ncol=k)
        for(i in seq(k)){
            pval1[,i] <- apply(YX1, 2, function(x){
                cor.test(x, score_cor2[, i])$p.value
            })
            pval2[,i] <- apply(YX2, 2, function(x){
                cor.test(x, score_cor1[, i])$p.value
            })
            qval1[,i] <- p.adjust(pval1[,i], "BH")
            qval2[,i] <- p.adjust(pval2[,i], "BH")
        }
        # Output
        list(res=res, loading1=loading1, loading2=loading2,
            score1=score1, score2=score2,
            score_cor1=score_cor1, score_cor2=score_cor2,
            cor1=cor1, cor2=cor2,
            pval1=pval1, pval2=pval2, qval1=qval1, qval2=qval2)
    }else{
        # Output
        list(res=res, loading1=loading1, loading2=loading2,
            score1=score1, score2=score2,
            score_cor1=score_cor1, score_cor2=score_cor2)
    }
}

####################################
# Setting for Plot (Parasitism1)
####################################
colorAt <- c(
    # 1d: red
    rgb(1,0,0), # wt root
    brewer.pal(8, "Dark2")[4], # wol root
    brewer.pal(8, "Dark2")[2], # wt w Pj
    brewer.pal(8, "Accent")[3], # wol w Pj
    # 3d: green
    rgb(0,1,0), # wt root
    brewer.pal(8, "Dark2")[1], # wol root
    brewer.pal(8, "Dark2")[5], # wt w Pj
    brewer.pal(8, "Accent")[1], # wol w Pj
    # 7d: blue
    rgb(0,0,1), # wt root
    brewer.pal(8, "Dark2")[3], # wol root
    brewer.pal(8, "Accent")[2], # wt w Pj
    brewer.pal(8, "Accent")[5] # wol w Pj
)

names(colorAt) <- c(
    # 1d: red
    "At_wt_root_1d",
    "At_wol_root_1d",
    "At_wt_1d",
    "At_wol_1d",
    # 3d: green
    "At_wt_root_3d",
    "At_wol_root_3d",
    "At_wt_3d",
    "At_wol_3d",
    # 7d: blue
    "At_wt_root_7d",
    "At_wol_root_7d",
    "At_wt_7d",
    "At_wol_7d")

colorAt2 <- c(
    rgb(1,0,0), # 1d: red
    rgb(0,1,0), # 3d: green
    rgb(0,0,1), # 7d: blue
    brewer.pal(8, "Dark2")[1], # wt
    brewer.pal(8, "Accent")[3], # wol
    brewer.pal(11, "RdBu")[1], # wPj
    brewer.pal(11, "RdBu")[11] # woPj
)

names(colorAt2) <- c(
    "1d",
    "3d",
    "7d",
    "wt",
    "wol",
    "wPj",
    "woPj")

colorPj <- c(
    # 1d: red
    rgb(1,0,0), # wt root
    brewer.pal(8, "Dark2")[2], # wt w Pj
    brewer.pal(8, "Accent")[3], # wol w Pj
    # 3d: green
    rgb(0,1,0), # wt root
    brewer.pal(8, "Dark2")[5], # wt w Pj
    brewer.pal(8, "Accent")[1], # wol w Pj
    # 7d: blue
    rgb(0,0,1), # wt root
    brewer.pal(8, "Accent")[2], # wt w Pj
    brewer.pal(8, "Accent")[5] # wol w Pj
)

names(colorPj) <- c(
    # 1d: red
    "Pj_root_1d",
    "Pj_wt_1d",
    "Pj_wol_1d",
    # 3d: green
    "Pj_root_3d",
    "Pj_wt_3d",
    "Pj_wol_3d",
    # 7d: blue
    "Pj_root_7d",
    "Pj_wt_7d",
    "Pj_wol_7d")

colorPj2 <- c(
    rgb(1,0,0), # 1d: red
    rgb(0,1,0), # 3d: green
    rgb(0,0,1), # 7d: blue
    brewer.pal(8, "Dark2")[1], # wt
    brewer.pal(8, "Accent")[3], # wol
    brewer.pal(11, "RdBu")[1], # wAt
    brewer.pal(11, "RdBu")[11] # woAt
)

names(colorPj2) <- c(
    "1d",
    "3d",
    "7d",
    "wt",
    "wol",
    "wAt",
    "woAt")

####################################
# Setting for Enrichment Analysis
####################################
.EnrichLoadings <- function(loading, gotable, expmat, qval=0.1){
    target_at <- unlist(sapply(rownames(expmat), function(x){
        which(gotable$GENEID == x)
    }))
    gotable <- gotable[target_at, ]
    uniq.goid <- unique(gotable$GOID)
    pos_position <- which(loading == 1)
    k <- length(pos_position)
    # Enrichment Analysis in each GO Term
    out <- do.call("rbind",
            lapply(uniq.goid, function(u){
                all_geneid <- gotable$GENEID[which(gotable$GOID == u)]
                sig_geneid <- intersect(all_geneid, rownames(expmat)[pos_position])
                if(length(sig_geneid) != 0){
                    x <- length(sig_geneid)
                    m <- length(all_geneid)
                    n <- length(loading) - m
                    p <- dhyper(x, m, n, k)
                    go <- unlist(gotable[which(gotable$GOID == u)[1], 3:4])
                    data.frame(GOID=u, Definition=go[1], GOTerm=go[2],
                        GeneID=paste(sig_geneid, collapse="|"), Pvalue=p)
                }else{
                    NULL
                }
            }))
    out <- data.frame(out, Qvalue=p.adjust(out$Pvalue, "BH"))
    # Sort
    out <- out[order(out$Qvalue), ]
    # Thresholding
    out[which(out$Qvalue < qval), ]
}
