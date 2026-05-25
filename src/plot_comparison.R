source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outfile <- args[2]

# Load gene set recovery results
infile <- paste0("output/", sample, "/evaluation/geneset.csv")
if(!file.exists(infile)){
    message("No geneset evaluation file found: ", infile)
    quit(status=0)
}

df <- read.csv(infile)

# Filter to F1 scores (non-NA K)
df_f1 <- df[!is.na(df$K), ]

####################################
# Bar plot: F1@K by method and species
####################################
df_f1$method <- factor(df_f1$method,
    levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))
df_f1$K_label <- paste0("K=", df_f1$K)

g1 <- ggplot(df_f1, aes(x=method, y=f1, fill=method)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(species ~ K_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          text = element_text(size=14)) +
    scale_fill_brewer(palette="Set2") +
    ylab("F1 Score (Procambium marker recovery)") +
    xlab("Method") +
    ggtitle(paste0("Gene Set Recovery: ", sample))

ggsave(file=paste0("plot/", sample, "/comparison_f1.png"),
    plot=g1, width=12, height=8)

####################################
# Bar plot: Precision and Recall at K=500
####################################
df_k500 <- df_f1[df_f1$K == 500, ]

df_pr <- rbind(
    data.frame(method=df_k500$method, species=df_k500$species,
        metric="Precision", value=df_k500$precision),
    data.frame(method=df_k500$method, species=df_k500$species,
        metric="Recall", value=df_k500$recall))

g2 <- ggplot(df_pr, aes(x=method, y=value, fill=metric)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~ species) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          text = element_text(size=14)) +
    scale_fill_brewer(palette="Set1") +
    ylab("Score") +
    xlab("Method") +
    ggtitle(paste0("Precision/Recall @K=500: ", sample))

ggsave(file=paste0("plot/", sample, "/comparison_pr.png"),
    plot=g2, width=10, height=6)

####################################
# AUROC comparison
####################################
df_auroc <- df[is.na(df$K) & !is.na(df$auroc), ]
if(nrow(df_auroc) > 0){
    df_auroc$method <- factor(df_auroc$method,
        levels=c("guided_pls", "rbh_pls", "guided_pca", "pca"))

    g3 <- ggplot(df_auroc, aes(x=method, y=auroc, fill=method)) +
        geom_bar(stat="identity") +
        facet_wrap(~ species) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
              text = element_text(size=14)) +
        scale_fill_brewer(palette="Set2") +
        ylab("AUROC") +
        xlab("Method") +
        ggtitle(paste0("AUROC (Procambium markers): ", sample))

    ggsave(file=paste0("plot/", sample, "/comparison_auroc.png"),
        plot=g3, width=8, height=6)
}

file.create(outfile)
