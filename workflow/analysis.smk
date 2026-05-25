rule guided_pls:
    input:
        'data/{sample}.RData'
    output:
        'output/{sample}/guided_pls.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/{sample}_guided_pls.txt'
    log:
        'logs/{sample}_guided_pls.log'
    shell:
        'src/guided_pls.sh {wildcards.sample} {output} >& {log}'

rule rbh_pls:
    input:
        'data/{sample}.RData'
    output:
        'output/{sample}/rbh_pls.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/{sample}_rbh_pls.txt'
    log:
        'logs/{sample}_rbh_pls.log'
    shell:
        'src/rbh_pls.sh {wildcards.sample} {output} >& {log}'

rule guided_pca:
    input:
        'data/{sample}.RData'
    output:
        'output/{sample}/guided_pca.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/{sample}_guided_pca.txt'
    log:
        'logs/{sample}_guided_pca.log'
    shell:
        'src/guided_pca.sh {wildcards.sample} {output} >& {log}'

rule pca:
    input:
        'data/{sample}.RData'
    output:
        'output/{sample}/pca.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/{sample}_pca.txt'
    log:
        'logs/{sample}_pca.log'
    shell:
        'src/pca.sh {wildcards.sample} {output} >& {log}'

rule deg:
    input:
        'data/{sample}.RData'
    output:
        'output/{sample}/deg/{species}_{degs}.RData'
    benchmark:
        'benchmarks/{sample}_deg_{species}_{degs}.txt'
    log:
        'logs/{sample}_deg_{species}_{degs}.log'
    shell:
        'src/deg.sh {wildcards.sample} {wildcards.species} {wildcards.degs} {output} >& {log}'
