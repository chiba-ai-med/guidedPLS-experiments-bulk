rule preprocess:
    output:
        'data/{sample}.RData'
    resources:
        mem_gb=50
    benchmark:
        'benchmarks/preprocess_{sample}.txt'
    log:
        'logs/preprocess_{sample}.log'
    shell:
        'src/preprocess_{wildcards.sample}.sh >& {log}'
