rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        bam="results/dedup/{sample}.deduped.sorted.bam",
        metrics="results/qc/dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
    conda:
        "../envs/picard.yml"
    shell:
        "picard MarkDuplicates -I {input} --METRICS_FILE {output.metrics} "
        "--OUTPUT {output.bam} 2> {log}"

rule samtools_index:
    input:
        "results/dedup/{sample}.deduped.sorted.bam",
    output:
        "results/dedup/{sample}.deduped.sorted.bai",
    log:
        "logs/samtools_index/{sample}.log",
    threads:
        config['main_config']['threads']
    conda:
        "../envs/mapping_tools.yml"    
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"
