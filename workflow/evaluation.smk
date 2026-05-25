rule evaluate_geneset:
    input:
        'data/{sample}.RData',
        expand('output/{{sample}}/{m}.RData', m=METHODS)
    output:
        'output/{sample}/evaluation/geneset.csv'
    benchmark:
        'benchmarks/{sample}_evaluate_geneset.txt'
    log:
        'logs/{sample}_evaluate_geneset.log'
    shell:
        'src/evaluate_geneset.sh {wildcards.sample} {output} >& {log}'

rule evaluate_factor_specificity:
    input:
        'data/{sample}.RData',
        expand('output/{{sample}}/{m}.RData', m=METHODS),
        lambda wc: expand('output/{sample}/deg/{s}_{d}.RData',
            sample=wc.sample, s=SPECIES,
            d=DEGS_MAP[wc.sample])
    output:
        'output/{sample}/evaluation/factor_specificity.csv'
    benchmark:
        'benchmarks/{sample}_evaluate_factor_specificity.txt'
    log:
        'logs/{sample}_evaluate_factor_specificity.log'
    shell:
        'src/evaluate_factor_specificity.sh {wildcards.sample} {output} >& {log}'

rule evaluate_enrichment:
    input:
        'data/{sample}.RData',
        'output/{sample}/{method}.RData'
    output:
        'output/{sample}/enrichment/finish_{method}'
    wildcard_constraints:
        method='|'.join([re.escape(x) for x in METHODS])
    benchmark:
        'benchmarks/{sample}_enrichment_{method}.txt'
    log:
        'logs/{sample}_enrichment_{method}.log'
    shell:
        'src/evaluate_enrichment.sh {wildcards.sample} {wildcards.method} {output} >& {log}'
