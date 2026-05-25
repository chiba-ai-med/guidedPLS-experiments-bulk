# CLAUDE.md

## プロジェクト概要

guided-PLSのバルクRNA-seqベンチマーク。寄生植物(Pj)×ホスト植物(At)のRNA-seqデータで、guided-PLSの種間統合における定量的優位性を示す。

このレポジトリは論文の一部で、以下の3レポジトリで構成される:
- `chiba-ai-med/guidedPLS-experiments-bulk` (本レポジトリ)
- `chiba-ai-med/guidedPLS-experiments-sc`
- `chiba-ai-med/ImageRegistration-experiments3`

## データセット

- **parasitism1**: NAIST提供、At(wt/wol) × Pj、時系列(1d/3d/7d)、寄生±
- **grafting**: Kurotani2020、At × Pj接木、時系列(0d-14d)

## 比較手法

| 手法 | 種間統合 | 因子ガイド |
|------|---------|-----------|
| guided-PLS | ○ | ○ |
| RBH-PLS | ○ (1:1オーソログのみ) | × |
| guided-PCA | × (種内) | ○ |
| PCA | × (種内) | × |

## パイプライン

Snakemakeパイプライン（77ジョブ）。Singularityコンテナで実行。

```
Snakefile          # メインエントリポイント（workflow/*.smkをinclude）
workflow/
  preprocess.smk   # 前処理
  analysis.smk     # guided-PLS, RBH-PLS, guided-PCA, PCA, DEG
  evaluation.smk   # geneset評価, factor specificity, GO enrichment
  plot.smk         # scatter, heatmap, sankey, comparison
  dag.sh           # DAG画像生成
  report.sh        # Snakemakeレポート生成
src/               # R/shellスクリプト（各ルールの実体）
```

Dockerイメージ: `docker://koki/guidedpls-experiments:20230502`

## 主な結果

- **因子特異性**: guided-PLS (0.77-0.85) >> 他手法 (0.50-0.81)
- **種間concordance**: guided-PLS 100% (全次元でAtとPjが同じ因子に対応)
- DEG回収率はguided-PLSが低いが、母集団サイズの違いによるバイアスが原因

## Figure構成

- **Figure 2**: guided-PLSの次元×条件 enrichmentヒートマップ (`src/plot_figure3.R`)
- **Figure 3**: 手法間定量比較 — factor specificity & concordance (`src/plot_figure4.R`)
- 出力先: `plot/Figure/{番号}/`

## ディレクトリ構成

```
data/              # 前処理済みRData
output/            # 解析結果（RData, CSV）
benchmarks/        # Snakemakeベンチマーク
logs/              # ログ
plot/              # 全プロット
  Figure/          # 論文用Figure
```

## 実行方法

```bash
# パイプライン全体
snakemake --use-singularity --singularity-args "-B /home/koki/Dev/guidedPLS-experiments/data/naist:/home/koki/Dev/guidedPLS-experiments/data/naist" -j 4

# DAG画像生成
bash workflow/dag.sh

# レポート生成
bash workflow/report.sh
```
