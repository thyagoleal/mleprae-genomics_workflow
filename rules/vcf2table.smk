rule bed_snps:
    input:
        "results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz"
    output:
        "results/variant_calls/variants_bed/{sample}.snps.bed.gz",
    log:
        "logs/vcf2bed/{sample}.log",
    threads:
        config['main_config']['threads']
    params:
        config['params']['vcf2bed']['memory']
    conda:
        "../envs/vcf2bed.yml"
    shell:
        "vcf2bed --snvs < <(zcat {input}) --max-mem {params} 2>> {log} | "
        " bgzip -@ {threads} > {output}"

rule bed_insertions:
    input:
        "results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz"
    output:
        "results/variant_calls/variants_bed/{sample}.insertions.bed.gz",
    log:
        "logs/vcf2bed/{sample}.log",
    threads:
        config['main_config']['threads']
    params:
        config['params']['vcf2bed']['memory']
    conda:
        "../envs/vcf2bed.yml"    
    shell:
        "vcf2bed --insertions < <(zcat {input}) --max-mem {params} 2>> {log} | "
        " bgzip -@ {threads} > {output}"

rule bed_deletions:
    input:
        "results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz"
    output:
        "results/variant_calls/variants_bed/{sample}.deletions.bed.gz",
    log:
        "logs/vcf2bed/{sample}.log",
    threads:
        config['main_config']['threads']
    params:
        config['params']['vcf2bed']['memory']
    conda:
        "../envs/vcf2bed.yml"    
    shell:
        "vcf2bed --deletions < <(zcat {input}) --max-mem {params} 2>> {log} | "
        " bgzip -@ {threads} > {output}"
