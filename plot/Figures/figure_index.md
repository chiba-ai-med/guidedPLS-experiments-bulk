# Figure Index — guidedPLS-experiments-bulk (Fig.4 in paper)

## Main Figures (plot/Figures/main/)

### Fig4A — Task overview (TO BE CREATED)
- **File**: Fig4A_task_overview.pdf/png
- **Candidate panel**: Fig.4a
- **Dataset**: grafting + parasitism1 (schematic)
- **Guide Z**: experimental condition labels (time-point, treatment)
- **Evaluation label**: DEG sets per condition
- **Method comparison**: N/A (diagram)
- **What it shows**: Host RNA-seq × parasite RNA-seq diagonal integration task. Two species with no shared gene space; guide Z provides condition metadata; evaluation by factor specificity and cross-species concordance.
- **Main/Supplementary**: Main
- **Notes**: Needs to be created. Should illustrate X1 (host), X2 (parasite), Y (shared condition labels), and the diagonal integration concept.

### Fig4B — Latent space comparison
- **File**: Fig4B_latent_{method}_grafting.png (×4 methods)
- **Candidate panel**: Fig.4b
- **Dataset**: grafting
- **Guide Z**: time-point + graft condition
- **Evaluation label**: visual inspection of species/condition clustering
- **Method comparison**: guided-PLS vs RBH-PLS vs guided-PCA vs PCA
- **What it shows**: Score scatter plots (Dim1 vs Dim2) for each method. guided-PLS aligns At and Pj samples by condition in a shared latent space; other methods show worse alignment or no cross-species structure.
- **Main/Supplementary**: Main
- **Notes**: May select 1-2 representative methods for main, rest to supplementary.

### Fig4C — Quantitative performance (factor specificity + concordance)
- **File**: Fig4C_quantitative_performance.pdf/png
- **Candidate panel**: Fig.4c
- **Dataset**: grafting + parasitism1
- **Guide Z**: condition labels
- **Evaluation label**: DEG-derived factor assignment per dimension
- **Method comparison**: guided-PLS vs RBH-PLS vs guided-PCA vs PCA
- **What it shows**: (a) Factor specificity — guided-PLS achieves highest mean specificity (0.77–0.85) across both datasets. (b) Cross-species concordance — guided-PLS achieves 100% (4/4 dimensions concordant) in both datasets; others range 0–55%.
- **Main/Supplementary**: Main
- **Notes**: Already composed as a two-panel bar chart. Key result of this experiment.

### Fig4C alt — Dimension-factor correspondence heatmap
- **File**: Fig4C_dimension_factor_heatmap.pdf/png
- **Candidate panel**: Fig.4c (alternative or additional panel)
- **Dataset**: parasitism1 + grafting
- **Guide Z**: condition labels
- **Evaluation label**: -log10(p) enrichment score per dimension × condition
- **Method comparison**: guided-PLS only
- **What it shows**: Each guided-PLS dimension corresponds to exactly one experimental condition, and this correspondence is consistent across At and Pj. Diagonal-like pattern in the heatmap.
- **Main/Supplementary**: Main
- **Notes**: Strong visual evidence for factor specificity and cross-species concordance.

### Fig4D — Label guidance effect (NOT AVAILABLE)
- **File**: Fig4D_label_guidance_effect.pdf/png
- **Candidate panel**: Fig.4d
- **Dataset**: —
- **Guide Z**: real labels vs shuffled labels vs no labels
- **Evaluation label**: —
- **Method comparison**: guided-PLS (with/without/shuffled labels)
- **What it shows**: Would demonstrate that guidance labels are essential, not just a trivial artifact.
- **Main/Supplementary**: Main
- **Notes**: No shuffled-label experiment exists in this repo. Consider adding, or omit this panel.

## Supplementary Figures (plot/Figures/supplementary/)

### Factor specificity per dataset
- **File**: factor_specificity_{dataset}.png
- **Dataset**: grafting, parasitism1 (separate)
- **What it shows**: Per-species (At/Pj) factor specificity bar charts for all 4 methods.

### Factor concordance per dataset
- **File**: factor_concordance_{dataset}.png
- **Dataset**: grafting, parasitism1 (separate)
- **What it shows**: Cross-species dimension concordance rate for all 4 methods.

### DEG × Loading heatmaps
- **File**: heatmap_{dataset}_{method}.png
- **Dataset**: grafting, parasitism1
- **Method**: guided_pls, rbh_pls, guided_pca, pca
- **What it shows**: Jaccard overlap between DEG sets and loading-based gene sets per dimension. Reveals which dimensions capture which biological signals.

### Score scatter plots (parasitism1)
- **File**: scatter_parasitism1_{method}_1_2.png
- **Dataset**: parasitism1
- **Method**: guided_pls, rbh_pls, guided_pca, pca
- **What it shows**: Latent score scatter (Dim1 vs Dim2) for the second dataset.

### Workflow DAG
- **File**: workflow_dag.png
- **What it shows**: Snakemake pipeline DAG.
