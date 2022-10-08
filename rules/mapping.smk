ruleorder: map_reads_pe > map_reads_se

rule map_reads_se:
    input:
        "results/trimmed/{sample}.trimmed.fastq.gz",
    output:
        temp("results/mapped/{sample}.bam"),
    log:
        log1="logs/bowtie2/{sample}.log",
        log2="logs/samtools_view/{sample}.log",
        log3="logs/fb_leftalign/{sample}.log",
    params:
        idx=config["ref"]["index_prefix"],
        rg=get_read_group,
        ref_fa=config["ref"]["ref-file_path"],
        max_mem_thread=config["params"]["samtools"]["max_mem_thread"],
        platform=config['experiment_opts']['platform']
    threads:
        config["main_config"]["threads"]
    conda:
        "../envs/mapping_tools.yml"
    shell:
        "bowtie2 -x {params.idx} -p {threads} -U {input} "
        "--rg-id {wildcards.sample} "
        "--rg 'SM:{wildcards.sample}\tPL:{params.platform}' 2>> {log.log1} | "
        "samtools view -Su -@ {threads} -F 4 - 2>> {log.log2} | "
        "bamleftalign -c -f {params.ref_fa} 2>> {log.log3} > {output}"

rule map_reads_pe:
    input:
        in1=["results/merged_reads/{sample}_R1.fastq.gz",
        "results/merged_reads/{sample}_R1.NotMerged.fastq.gz"],
        in2=["results/merged_reads/{sample}_R2.fastq.gz",
        "results/merged_reads/{sample}_R2.NotMerged.fastq.gz"],
        U=["results/merged_reads/{sample}.merged.fastq.gz",
        "results/trimmed/{sample}_R1.unpaired.fastq.gz",
        "results/trimmed/{sample}_R2.unpaired.fastq.gz"]
    output:
        temp("results/mapped/{sample}.bam"),
    log:
        log1="logs/bowtie2/{sample}.log",
        log2="logs/samtools_view/{sample}.log",
        log3="logs/bamleftalign/{sample}.log",
    params:
        idx=config["ref"]["index_prefix"],
        rg=get_read_group,
        ref_fa=config["ref"]["ref-file_path"],
        max_mem_thread=config["params"]["samtools"]["max_mem_thread"],
        platform=config['experiment_opts']['platform'],
    threads: config["main_config"]["threads"]
    conda:
        "../envs/mapping_tools.yml"
    shell:
        "bowtie2 -x {params.idx} -p {threads} -U {input.U} -1 {input.in1} "
        "-2 {input.in2} --rg-id {wildcards.sample} "
        "--rg 'SM:{wildcards.sample}\tPL:{params.platform}' 2>> {log.log1} | "
        "samtools view - -@ {threads} -Su -F 4 2>> {log.log2} | "
        "bamleftalign -c -f {params.ref_fa} 2>> {log.log3} > {output}"

rule sort_bam:
    input:
        "results/mapped/{sample}.bam"
    output:
        temp("results/mapped/{sample}.sorted.bam")
    threads:
        config["main_config"]["threads"]
    params:
        max_mem_thread=config["params"]["samtools"]["max_mem_thread"],
    log:
        "logs/samtools_sort/{sample}.log",
    conda:
        "../envs/mapping_tools.yml"    
    shell:
        "samtools sort {input} -@ {threads} -m {params.max_mem_thread} -O BAM -o {output} 2>> {log}"

rule mapping_stats:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        "results/qc/samtools-stats/{sample}.txt",
    log:
        "logs/samtools-stats/{sample}.log",
    threads:
        config["main_config"]["threads"]
    conda:
        "../envs/mapping_tools.yml"    
    shell:
        "samtools stats {input} -@ {threads} > {output} 2> {log}"

rule flag_stats:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        "results/qc/samtools-flagstat/{sample}.txt",
    log:
        "logs/samtools-flagstat/{sample}.log",
    threads:
        config["main_config"]["threads"]
    conda:
        "../envs/mapping_tools.yml"    
    shell:
        "samtools flagstat {input} -@ {threads} > {output} 2> {log}"
