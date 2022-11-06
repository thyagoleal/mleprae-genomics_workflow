ruleorder: trim_reads_PE > trim_reads_SE

rule trim_reads_PE:
    input:
        r1="data/samples/prepared-fastqs/{sample}_R1.fastq.gz",
        r2="data/samples/prepared-fastqs/{sample}_R2.fastq.gz",
    output:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        r1_unpaired="results/trimmed/{sample}_R1.unpaired.fastq.gz",
        r2_unpaired="results/trimmed/{sample}_R2.unpaired.fastq.gz",
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra="-phred33",
    log:
        "logs/trimmomatic/{sample}.log",
    threads: config["main_config"]["threads"]
    # conda: 
    #     "../envs/qc_trimming.yml"
    # wrapper:
    #     "v1.7.0/bio/trimmomatic/pe"
    # script:
    #     "../scripts/trimmomatic_pe_wrapper.py"
    shell:
        "trimmomatic PE -threads {threads} {params.extra} {input.r1} {input.r2} "
        "{output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} "
        " {params} 2> {log}"
        

rule trim_reads_SE:
    input:
        "data/samples/prepared-fastqs/{sample}.fastq.gz",
    output:
        "results/trimmed/{sample}.trimmed.fastq.gz",
    params:
        **config["params"]["trimmomatic"]["se"],
    log:
        "logs/trimmomatic/{sample}.log"
    threads:
        config["main_config"]["threads"]
    # conda: 
    #     "../envs/qc_trimming.yml"
    # wrapper:
    #     "v1.7.0/bio/trimmomatic/se"
    # script:
    #     "../scripts/trimmomatic_se_wrapper.py"
    shell:
        "trimmomatic SE -phred33 -threads {threads} {input} {output} {params} 2> {log}"
