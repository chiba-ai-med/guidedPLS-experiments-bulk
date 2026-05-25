source("src/Functions.R")

####################################
# Gene Expression Matrix
####################################
# Load Kurotani2020 TPM data
At_all <- read.table("data/naist/Kurotani2020_Datta/At_TPM_Kurotani2020.txt",
    header=TRUE, row.names=1) |> as.matrix()
Pj_all <- read.table("data/naist/Kurotani2020_Datta/Pj-pep_TPM_Kurotani2020.txt",
    header=TRUE, row.names=1) |> as.matrix()

# Select Grafting experiment samples only
# At-relevant: Gra00at (day0 At alone) + Gra01bo-Gra14bo (day1-14 both)
at_cols <- grep("Gra00at|Gra0[1-7]bo|Gra14bo", colnames(At_all))
# Pj-relevant: Gra00pj (day0 Pj alone) + Gra01bo-Gra14bo (day1-14 both)
pj_cols <- grep("Gra00pj|Gra0[1-7]bo|Gra14bo", colnames(Pj_all))

X_At_raw <- At_all[, at_cols]
X_Pj_raw <- Pj_all[, pj_cols]

# Log conversion
X_At <- log10(X_At_raw + 1)
X_Pj <- log10(X_Pj_raw + 1)

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
counts_At_all <- read.table("data/naist/Kurotani2020_Datta/Kurotani_At_count_new.txt",
    header=TRUE, row.names=1)
counts_Pj_all <- read.table("data/naist/Kurotani2020_Datta/Kurotani2020_count_Pj_pep.txt",
    header=TRUE, row.names=1)

counts_At <- counts_At_all[, at_cols]
counts_Pj <- counts_Pj_all[, pj_cols]

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
# Samples: Gra00at{1-4} (day0, alone) + Gra01bo{1-4}...Gra14bo{1-4} (day1-14, grafted)
####################################
n_At <- ncol(X_At)
cn_At <- colnames(X_At)

# Time indicator (10 timepoints: 0d,1d,2d,3d,4d,5d,6d,7d,14d)
time_labels <- c("0d", "1d", "2d", "3d", "4d", "5d", "6d", "7d", "14d")
Y_At_time <- matrix(0, nrow=n_At, ncol=length(time_labels))
rownames(Y_At_time) <- cn_At
colnames(Y_At_time) <- time_labels

Y_At_time[grep("Gra00at", cn_At), "0d"] <- 1
Y_At_time[grep("Gra01bo", cn_At), "1d"] <- 1
Y_At_time[grep("Gra02bo", cn_At), "2d"] <- 1
Y_At_time[grep("Gra03bo", cn_At), "3d"] <- 1
Y_At_time[grep("Gra04bo", cn_At), "4d"] <- 1
Y_At_time[grep("Gra05bo", cn_At), "5d"] <- 1
Y_At_time[grep("Gra06bo", cn_At), "6d"] <- 1
Y_At_time[grep("Gra07bo", cn_At), "7d"] <- 1
Y_At_time[grep("Gra14bo", cn_At), "14d"] <- 1

# Graft indicator (alone vs grafted)
Y_At_graft <- matrix(0, nrow=n_At, ncol=2)
rownames(Y_At_graft) <- cn_At
colnames(Y_At_graft) <- c("grafted", "alone")

Y_At_graft[grep("Gra00at", cn_At), "alone"] <- 1
Y_At_graft[grep("bo", cn_At), "grafted"] <- 1

Y_At_all <- cbind(Y_At_time, Y_At_graft)

####################################
# Indicator Matrix for Pj
# Samples: Gra00pj{1-4} (day0, alone) + Gra01bo{1-4}...Gra14bo{1-4} (day1-14, grafted)
####################################
n_Pj <- ncol(X_Pj)
cn_Pj <- colnames(X_Pj)

Y_Pj_time <- matrix(0, nrow=n_Pj, ncol=length(time_labels))
rownames(Y_Pj_time) <- cn_Pj
colnames(Y_Pj_time) <- time_labels

Y_Pj_time[grep("Gra00pj", cn_Pj), "0d"] <- 1
Y_Pj_time[grep("Gra01bo", cn_Pj), "1d"] <- 1
Y_Pj_time[grep("Gra02bo", cn_Pj), "2d"] <- 1
Y_Pj_time[grep("Gra03bo", cn_Pj), "3d"] <- 1
Y_Pj_time[grep("Gra04bo", cn_Pj), "4d"] <- 1
Y_Pj_time[grep("Gra05bo", cn_Pj), "5d"] <- 1
Y_Pj_time[grep("Gra06bo", cn_Pj), "6d"] <- 1
Y_Pj_time[grep("Gra07bo", cn_Pj), "7d"] <- 1
Y_Pj_time[grep("Gra14bo", cn_Pj), "14d"] <- 1

Y_Pj_graft <- matrix(0, nrow=n_Pj, ncol=2)
rownames(Y_Pj_graft) <- cn_Pj
colnames(Y_Pj_graft) <- c("grafted", "alone")

Y_Pj_graft[grep("Gra00pj", cn_Pj), "alone"] <- 1
Y_Pj_graft[grep("bo", cn_Pj), "grafted"] <- 1

Y_Pj_all <- cbind(Y_Pj_time, Y_Pj_graft)

####################################
# Setting for Plot
####################################
# Color scheme for grafting: timepoint-based gradient
colorGraft_At <- c(
    brewer.pal(9, "YlOrRd")[1],  # 0d alone
    brewer.pal(9, "YlOrRd")[2],  # 1d
    brewer.pal(9, "YlOrRd")[3],  # 2d
    brewer.pal(9, "YlOrRd")[4],  # 3d
    brewer.pal(9, "YlOrRd")[5],  # 4d
    brewer.pal(9, "YlOrRd")[6],  # 5d
    brewer.pal(9, "YlOrRd")[7],  # 6d
    brewer.pal(9, "YlOrRd")[8],  # 7d
    brewer.pal(9, "YlOrRd")[9]   # 14d
)
names(colorGraft_At) <- c("At_0d", "At_1d", "At_2d", "At_3d", "At_4d",
    "At_5d", "At_6d", "At_7d", "At_14d")

colorGraft_Pj <- c(
    brewer.pal(9, "YlGnBu")[1],
    brewer.pal(9, "YlGnBu")[2],
    brewer.pal(9, "YlGnBu")[3],
    brewer.pal(9, "YlGnBu")[4],
    brewer.pal(9, "YlGnBu")[5],
    brewer.pal(9, "YlGnBu")[6],
    brewer.pal(9, "YlGnBu")[7],
    brewer.pal(9, "YlGnBu")[8],
    brewer.pal(9, "YlGnBu")[9]
)
names(colorGraft_Pj) <- c("Pj_0d", "Pj_1d", "Pj_2d", "Pj_3d", "Pj_4d",
    "Pj_5d", "Pj_6d", "Pj_7d", "Pj_14d")

label_At <- rep(0, length=n_At)
names(label_At) <- c(
    rep("At_0d", 4),
    rep("At_1d", 4), rep("At_2d", 4), rep("At_3d", 4), rep("At_4d", 4),
    rep("At_5d", 4), rep("At_6d", 4), rep("At_7d", 4), rep("At_14d", 4))
for(nm in names(colorGraft_At)){
    label_At[which(names(label_At) == nm)] <- colorGraft_At[nm]
}

label_Pj <- rep(0, length=n_Pj)
names(label_Pj) <- c(
    rep("Pj_0d", 4),
    rep("Pj_1d", 4), rep("Pj_2d", 4), rep("Pj_3d", 4), rep("Pj_4d", 4),
    rep("Pj_5d", 4), rep("Pj_6d", 4), rep("Pj_7d", 4), rep("Pj_14d", 4))
for(nm in names(colorGraft_Pj)){
    label_Pj[which(names(label_Pj) == nm)] <- colorGraft_Pj[nm]
}

####################################
# DEG conditions
####################################
deg_conditions <- c("0d", "1d", "2d", "3d", "4d", "5d", "6d", "7d", "14d", "graft")

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
    Y_At_time, Y_At_graft, Y_At_all,
    Y_Pj_time, Y_Pj_graft, Y_Pj_all,
    label_At, label_Pj,
    deg_conditions,
    file="data/grafting.RData")
