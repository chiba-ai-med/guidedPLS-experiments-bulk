# guidedPLS-experiments-bulk

Benchmarking of guidedPLS and comparison methods (RBH-PLS, guided-PCA, PCA) on bulk RNA-seq data from parasite-host cross-species experiments.

This repository is one of three sibling repositories that together support a single manuscript on guided-PLS:

- [`chiba-ai-med/guidedPLS-experiments-bulk`](https://github.com/chiba-ai-med/guidedPLS-experiments-bulk) — **this repository** (bulk RNA-seq, parasite-host)
- [`chiba-ai-med/guidedPLS-experiments-sc`](https://github.com/chiba-ai-med/guidedPLS-experiments-sc) — single-cell multi-omics experiments
- [`chiba-ai-med/ImageRegistration-experiments3`](https://github.com/chiba-ai-med/ImageRegistration-experiments3) — image-registration experiments

## Overview

This pipeline evaluates cross-species factor analysis methods using two datasets:

- **Parasitism**: *Phelipanche aegyptiaca* on *Arabidopsis thaliana* (NAIST)
- **Grafting**: *Phtheirospermum japonicum* on *Arabidopsis thaliana* (Kurotani et al., 2020)

The following four methods are compared:

| Method | Cross-species integration | Factor guidance |
|--------|--------------------------|-----------------|
| guided-PLS | Yes | Yes |
| RBH-PLS | Yes (1:1 orthologs only) | No |
| guided-PCA | No (within-species) | Yes |
| PCA | No (within-species) | No |

## Workflow

This workflow consists of 4 workflow modules as follows:

- **workflow/preprocess.smk**: Data preprocessing

- **workflow/analysis.smk**: Factor analysis (guided-PLS, RBH-PLS, guided-PCA, PCA) and DEG detection

- **workflow/evaluation.smk**: Gene set evaluation, factor specificity, and GO enrichment analysis

- **workflow/plot.smk**: Score scatter plots, DEG x Loading heatmaps, Sankey plots, and comparison plots

![](https://github.com/chiba-ai-med/guidedPLS-experiments-bulk/blob/main/plot/dag.png?raw=true)

## Requirements

- Bash: GNU bash
- Snakemake: >= 6.5.3
- Singularity: >= 3.8.0

Container image: `docker://koki/guidedpls-experiments:20230502`

## How to reproduce this workflow

### In Local Machine

```bash
snakemake -j 4 --use-singularity
```

### In Open Grid Engine

```bash
snakemake -j 32 --cluster qsub --latency-wait 600 --use-singularity
```

### In Slurm

```bash
snakemake -j 32 --cluster sbatch --latency-wait 600 --use-singularity
```

## License

Copyright (c) 2026 Koki Tsuyuzaki and Artificial Intelligence Medicine, released under the [MIT License](https://opensource.org/licenses/MIT).

## Authors

- Koki Tsuyuzaki
