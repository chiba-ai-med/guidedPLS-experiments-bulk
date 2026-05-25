rule plot_scores:
    input:
        'data/{sample}.RData',
        'output/{sample}/{method}.RData'
    output:
        'plot/{sample}/scatter/{method}/finish'
    wildcard_constraints:
        method='|'.join([re.escape(x) for x in METHODS])
    benchmark:
        'benchmarks/{sample}_plot_scores_{method}.txt'
    log:
        'logs/{sample}_plot_scores_{method}.log'
    shell:
        'src/plot_scores.sh {wildcards.sample} {wildcards.method} {output} >& {log}'

rule plot_heatmap:
    input:
        'data/{sample}.RData',
        'output/{sample}/{method}.RData',
        lambda wc: expand('output/{sample}/deg/{s}_{d}.RData',
            sample=wc.sample, s=SPECIES,
            d=DEGS_MAP[wc.sample])
    output:
        'plot/{sample}/heatmap/{method}.png'
    wildcard_constraints:
        method='|'.join([re.escape(x) for x in METHODS])
    benchmark:
        'benchmarks/{sample}_plot_heatmap_{method}.txt'
    log:
        'logs/{sample}_plot_heatmap_{method}.log'
    shell:
        'src/plot_heatmap.sh {wildcards.sample} {wildcards.method} {output} >& {log}'

rule plot_sankey:
    input:
        'data/{sample}.RData',
        'output/{sample}/{method}.RData',
        lambda wc: expand('output/{sample}/deg/{s}_{d}.RData',
            sample=wc.sample, s=SPECIES,
            d=DEGS_MAP[wc.sample])
    output:
        'plot/{sample}/sankey/{method}/finish'
    wildcard_constraints:
        method='|'.join([re.escape(x) for x in METHODS])
    benchmark:
        'benchmarks/{sample}_plot_sankey_{method}.txt'
    log:
        'logs/{sample}_plot_sankey_{method}.log'
    shell:
        'src/plot_sankey.sh {wildcards.sample} {wildcards.method} {output} >& {log}'

rule plot_comparison:
    input:
        'output/{sample}/evaluation/geneset.csv'
    output:
        'plot/{sample}/comparison/finish'
    benchmark:
        'benchmarks/{sample}_plot_comparison.txt'
    log:
        'logs/{sample}_plot_comparison.log'
    shell:
        'src/plot_comparison.sh {wildcards.sample} {output} >& {log}'
