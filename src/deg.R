source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
species <- args[2]
degs <- args[3]
outfile <- args[4]

load(paste0("data/", sample, ".RData"))

# Select species
if(species == "at"){
    counts <- counts_At
}
if(species == "pj"){
    counts <- counts_Pj
}

# Group Index
if(sample == "parasitism1"){
    # Parasitism1 dataset
    if((species == "at") && (degs == "1d")){
        index1 <- which(Y_At_time[,1] == 1)
    }
    if((species == "at") && (degs == "3d")){
        index1 <- which(Y_At_time[,2] == 1)
    }
    if((species == "at") && (degs == "7d")){
        index1 <- which(Y_At_time[,3] == 1)
    }
    if((species == "at") && (degs == "wol")){
        index1 <- which(Y_At_wol[,1] == 1)
    }
    if((species == "at") && (degs == "parasm")){
        index1 <- which(Y_At_parasm[,1] == 1)
    }
    if((species == "pj") && (degs == "1d")){
        index1 <- which(Y_Pj_time[,1] == 1)
    }
    if((species == "pj") && (degs == "3d")){
        index1 <- which(Y_Pj_time[,2] == 1)
    }
    if((species == "pj") && (degs == "7d")){
        index1 <- which(Y_Pj_time[,3] == 1)
    }
    if((species == "pj") && (degs == "wol")){
        index1 <- which(Y_Pj_wol[,1] == 1)
    }
    if((species == "pj") && (degs == "parasm")){
        index1 <- which(Y_Pj_parasm[,1] == 1)
    }
}

if(sample == "grafting"){
    # Grafting dataset
    time_map <- c("0d"=1, "1d"=2, "2d"=3, "3d"=4, "4d"=5,
                  "5d"=6, "6d"=7, "7d"=8, "14d"=9)
    if(species == "at"){
        if(degs == "graft"){
            index1 <- which(Y_At_graft[,1] == 1)
        } else {
            idx <- time_map[degs]
            index1 <- which(Y_At_time[,idx] == 1)
        }
    }
    if(species == "pj"){
        if(degs == "graft"){
            index1 <- which(Y_Pj_graft[,1] == 1)
        } else {
            idx <- time_map[degs]
            index1 <- which(Y_Pj_time[,idx] == 1)
        }
    }
}

index2 <- setdiff(seq(ncol(counts)), index1)

# edgeR
group <- rep("A", length=ncol(counts))
group[index1] <- "B"
group <- factor(group)
design <- model.matrix(~ group)
d <- DGEList(counts = counts, group = group)
d <- calcNormFactors(d)
d <- estimateDisp(d, design)
fit <- glmFit(d, design)
lrt <- glmLRT(fit, coef = 2)
deg <- topTags(lrt, n=nrow(lrt))

# Save
save(deg, file=outfile)
