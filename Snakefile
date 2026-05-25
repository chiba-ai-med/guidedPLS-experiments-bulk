import re
from snakemake.utils import min_version

#
# Setting
#
min_version("6.5.3")
container: 'docker://koki/guidedpls-experiments:20230502'

SAMPLES = ['parasitism1', 'grafting']
METHODS = ['guided_pls', 'rbh_pls', 'guided_pca', 'pca']
SPECIES = ['at', 'pj']

# DEG conditions per sample
DEGS_parasitism1 = ['1d', '3d', '7d', 'wol', 'parasm']
DEGS_grafting = ['0d', '1d', '2d', '3d', '4d', '5d', '6d', '7d', '14d', 'graft']
DEGS_MAP = {
    'parasitism1': DEGS_parasitism1,
    'grafting': DEGS_grafting
}

#
# Include workflow modules
#
include: 'workflow/preprocess.smk'
include: 'workflow/analysis.smk'
include: 'workflow/evaluation.smk'
include: 'workflow/plot.smk'

#
# Target rule
#
def all_deg_files(sample):
    return expand('output/{sample}/deg/{s}_{d}.RData',
        sample=sample, s=SPECIES, d=DEGS_MAP[sample])

rule all:
    input:
        # Preprocess
        expand('data/{sample}.RData', sample=SAMPLES),
        # Analysis methods
        expand('output/{sample}/{m}.RData', sample=SAMPLES, m=METHODS),
        # DEG
        [f for s in SAMPLES for f in all_deg_files(s)],
        # Score scatter plots
        expand('plot/{sample}/scatter/{m}/finish', sample=SAMPLES, m=METHODS),
        # DEG x Loading heatmaps
        expand('plot/{sample}/heatmap/{m}.png', sample=SAMPLES, m=METHODS),
        # Sankey plots
        expand('plot/{sample}/sankey/{m}/finish', sample=SAMPLES, m=METHODS),
        # Gene set evaluation
        expand('output/{sample}/evaluation/geneset.csv', sample=SAMPLES),
        # GO enrichment
        expand('output/{sample}/enrichment/finish_{m}', sample=SAMPLES, m=METHODS),
        # Factor specificity
        expand('output/{sample}/evaluation/factor_specificity.csv', sample=SAMPLES),
        # Comparison plots
        expand('plot/{sample}/comparison/finish', sample=SAMPLES)
