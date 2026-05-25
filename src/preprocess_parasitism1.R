source("src/Functions.R")

####################################
# Gene Expression Matrix
####################################
# Log conversion
read.csv("data/naist/filtered_At_tpm.csv", header=TRUE, row.names=1) |>
    as.matrix() |> (`+`)(1) |> log10() -> X_At
read.csv("data/naist/filtered_Pj_tpm.csv", header=TRUE, row.names=1) |>
    as.matrix() |> (`+`)(1) |> log10() -> X_Pj

# High Variance Genes
targetAt <- which(rank(1/apply(X_At, 1, var)) <= 5000)
targetPj <- which(rank(1/apply(X_Pj, 1, var)) <= 5000)

X_At <- X_At[targetAt, ]
X_Pj <- X_Pj[targetPj, ]

# Centering
scaled_X_At <- t(scale(t(X_At), center=TRUE, scale=FALSE))
scaled_X_Pj <- t(scale(t(X_Pj), center=TRUE, scale=FALSE))

####################################
# RBH
####################################
orthtable <- unique(read.table("data/naist/Tosend20211025/TAIR10_Pjv1_RBH2.txt")[,1:2])

targetRBH <- unlist(apply(orthtable, 1, function(x){
    position_at <- which(x[1] == rownames(X_At))
    position_pj <- which(x[2] == rownames(X_Pj))
    if((length(position_at) == 1) && (length(position_pj) == 1)){
        c(rownames(X_At)[position_at], rownames(X_Pj)[position_pj])
    }
}))
targetRBH <- as.matrix(targetRBH)
dim(targetRBH) <- c(2, nrow(targetRBH)/2)
targetRBH <- t(targetRBH)

X_At_RBH <- X_At[targetRBH[,1], ]
X_Pj_RBH <- X_Pj[targetRBH[,2], ]

scaled_X_At_RBH <- t(scale(t(X_At_RBH), center=TRUE, scale=FALSE))
scaled_X_Pj_RBH <- t(scale(t(X_Pj_RBH), center=TRUE, scale=FALSE))

####################################
# Count table for DEGs
####################################
counts_all <- read.table("data/naist/counts_all.txt", header=TRUE)

counts_At <- counts_all[30338:nrow(counts_all), c(1:9, 22:ncol(counts_all))]
counts_Pj <- counts_all[seq(30337), c(10:21, 22:47)]

counts_At <- counts_At[intersect(names(targetAt), rownames(counts_At)), ]
counts_Pj <- counts_Pj[intersect(names(targetPj), rownames(counts_Pj)), ]

####################################
# Procambium gene set
####################################
procambium_At <- unique(read.csv("data/naist/Tosend20211025/Wendrich2020/Procambium.csv"))
procambium_Pj <- unique(merge(orthtable, procambium_At, by.x="V1", by.y="gene_AT"))
colnames(procambium_At)[1:2] <- c("GENEID", "SYMBOL")
colnames(procambium_Pj)[1:3] <- c("GENEID_At", "GENEID", "SYMBOL")

####################################
# Gene Ontology
####################################
At_GO <- unique(read.delim("data/naist/At_GO.list", header=FALSE, sep="\t"))
Pj_GO <- unique(read.delim("data/naist/Pj_GO.list", header=FALSE, sep="\t"))
colnames(At_GO) <- c("GENEID", "GOID")
colnames(Pj_GO) <- c("GENEID", "GOID")

gotable <- select(GO.db, columns=columns(GO.db), keytype="GOID", keys=keys(GO.db))
gotable_BP <- gotable[which(gotable$ONTOLOGY == "BP"), ]
gotable_MF <- gotable[which(gotable$ONTOLOGY == "MF"), ]
gotable_CC <- gotable[which(gotable$ONTOLOGY == "CC"), ]

At_GO_BP <- merge(gotable_BP, At_GO, by="GOID")
At_GO_MF <- merge(gotable_MF, At_GO, by="GOID")
At_GO_CC <- merge(gotable_CC, At_GO, by="GOID")

Pj_GO_BP <- merge(gotable_BP, Pj_GO, by="GOID")
Pj_GO_MF <- merge(gotable_MF, Pj_GO, by="GOID")
Pj_GO_CC <- merge(gotable_CC, Pj_GO, by="GOID")

####################################
# Indicator Matrix for At
####################################
Y_At_time <- matrix(0, nrow=ncol(X_At), ncol=3)
Y_At_wol <- matrix(0, nrow=ncol(X_At), ncol=2)
Y_At_parasm <- matrix(0, nrow=ncol(X_At), ncol=2)

rownames(Y_At_time) <- colnames(X_At)
rownames(Y_At_wol) <- colnames(X_At)
rownames(Y_At_parasm) <- colnames(X_At)

colnames(Y_At_time) <- c("1d", "3d", "7d")
colnames(Y_At_wol) <- c("wt", "wol")
colnames(Y_At_parasm) <- c("wPj", "woPj")

Y_At_time[grep("1d", rownames(Y_At_time)), 1] <- 1
Y_At_time[grep("3d", rownames(Y_At_time)), 2] <- 1
Y_At_time[grep("7d", rownames(Y_At_time)), 3] <- 1

targewol <- grep("wol", rownames(Y_At_wol))
Y_At_wol[targewol, 1] <- 1
Y_At_wol[setdiff(seq(nrow(Y_At_wol)), targewol), 2] <- 1

targetParasm <- grep("^Pj", rownames(Y_At_parasm))
Y_At_parasm[targetParasm, 1] <- 1
Y_At_parasm[setdiff(seq(nrow(Y_At_parasm)), targetParasm), 2] <- 1

Y_At_all <- cbind(Y_At_time, Y_At_wol, Y_At_parasm)

####################################
# Indicator Matrix for Pj
####################################
Y_Pj_time <- matrix(0, nrow=ncol(X_Pj), ncol=3)
Y_Pj_wol <- matrix(0, nrow=ncol(X_Pj), ncol=2)
Y_Pj_parasm <- matrix(0, nrow=ncol(X_Pj), ncol=2)

rownames(Y_Pj_time) <- colnames(X_Pj)
rownames(Y_Pj_wol) <- colnames(X_Pj)
rownames(Y_Pj_parasm) <- colnames(X_Pj)

colnames(Y_Pj_time) <- c("1d", "3d", "7d")
colnames(Y_Pj_wol) <- c("wt", "wol")
colnames(Y_Pj_parasm) <- c("wAt", "woAt")

Y_Pj_time[grep("1d", rownames(Y_Pj_time)), 1] <- 1
Y_Pj_time[grep("3d", rownames(Y_Pj_time)), 2] <- 1
Y_Pj_time[grep("7d", rownames(Y_Pj_time)), 3] <- 1

targewol <- grep("wol", rownames(Y_Pj_wol))
Y_Pj_wol[targewol, 1] <- 1
Y_Pj_wol[setdiff(seq(nrow(Y_Pj_wol)), targewol), 2] <- 1

# Pj root samples (first 12) = without At, remaining = with At
Y_Pj_parasm[1:12, 2] <- 1
Y_Pj_parasm[13:38, 1] <- 1

Y_Pj_all <- cbind(Y_Pj_time, Y_Pj_wol, Y_Pj_parasm)

####################################
# Setting for Plot
####################################
label_At <- rep(0, length=ncol(X_At))
names(label_At) <- c(
    rep("At_wt_root_1d", 3),
    rep("At_wt_root_3d", 3),
    rep("At_wt_root_7d", 3),
    rep("At_wt_1d", 4),
    rep("At_wt_3d", 4),
    rep("At_wt_7d", 3),
    rep("At_wol_1d", 4),
    rep("At_wol_3d", 5),
    rep("At_wol_7d", 6),
    rep("At_wol_root_1d", 3),
    rep("At_wol_root_3d", 3),
    rep("At_wol_root_7d", 4))

for(nm in names(colorAt)){
    label_At[which(names(label_At) == nm)] <- colorAt[nm]
}

label_Pj <- rep(0, length=ncol(X_Pj))
names(label_Pj) <- c(
    rep("Pj_root_1d", 4),
    rep("Pj_root_3d", 4),
    rep("Pj_root_7d", 4),
    rep("Pj_wt_1d", 4),
    rep("Pj_wt_3d", 4),
    rep("Pj_wt_7d", 3),
    rep("Pj_wol_1d", 4),
    rep("Pj_wol_3d", 5),
    rep("Pj_wol_7d", 6))

for(nm in names(colorPj)){
    label_Pj[which(names(label_Pj) == nm)] <- colorPj[nm]
}

####################################
# DEG conditions
####################################
deg_conditions <- c("1d", "3d", "7d", "wol", "parasm")

####################################
# Output
####################################
save(
    X_At, X_Pj,
    scaled_X_At, scaled_X_Pj,
    X_At_RBH, X_Pj_RBH,
    scaled_X_At_RBH, scaled_X_Pj_RBH,
    counts_At, counts_Pj,
    procambium_At, procambium_Pj,
    orthtable,
    At_GO_BP, At_GO_MF, At_GO_CC,
    Pj_GO_BP, Pj_GO_MF, Pj_GO_CC,
    Y_At_time, Y_At_wol, Y_At_parasm, Y_At_all,
    Y_Pj_time, Y_Pj_wol, Y_Pj_parasm, Y_Pj_all,
    label_At, label_Pj,
    deg_conditions,
    file="data/parasitism1.RData")
