ruleorder: fastqc_se_raw > fastqc_raw_pe_r1 > fastqc_raw_pe_r2 > fastqc_trimmomatic_se > fastqc_trimmomatic_pe_r1 > fastqc_trimmomatic_pe_r2

rule fastqc_raw_pe_r1:
    input:
        "data/samples/prepared-fastqs/{sample}_R1.fastq.gz",
    output:
        html="results/qc/fastqc_raw/{sample}_R1.html",
        zip="results/qc/fastqc_raw/{sample}_R1_fastqc.zip",
    params: "--quiet"
    log:
        "logs/fastqc_raw/{sample}_R1.log",
    threads: config['main_config']['threads']
    # conda: 
    #     "../envs/qc_trimming.yml"
    wrapper:
        "v1.8.0/bio/fastqc"
    # script:
    #     "../scripts/fastqc_wrapper.py"
    # shell:
        # "fastqc {input} --outdir results/qc/fastqc_raw/ {params} > {log} 2>&1"
               

rule fastqc_raw_pe_r2:
    input:
        "data/samples/prepared-fastqs/{sample}_R2.fastq.gz",
    output:
        html="results/qc/fastqc_raw/{sample}_R2.html",
        zip="results/qc/fastqc_raw/{sample}_R2_fastqc.zip",
    params: "--quiet"
    log:
        "logs/fastqc_raw/{sample}_R2.log",
    threads: config['main_config']['threads']
    # conda: 
    #     "../envs/qc_trimming.yml"
    wrapper:
        "v1.8.0/bio/fastqc"
    # script:
    #     "../scripts/fastqc_wrapper.py"
    # shell:
        # "fastqc {input} --outdir results/qc/fastqc_raw/ {params} > {log} 2>&1"

rule fastqc_se_raw:
    input:
        "data/samples/prepared-fastqs/{sample}.fastq.gz",
    output:
        html="results/qc/fastqc_raw/{sample}.html",
        zip="results/qc/fastqc_raw/{sample}_fastqc.zip"
    params: "--quiet"
    # conda: 
    #     "../envs/qc_trimming.yml"
    log:
        "logs/fastqc_raw/{sample}.log",
    threads: config['main_config']['threads']
    wrapper:
        "v1.8.0/bio/fastqc"   
    # script:
    #     "../scripts/fastqc_wrapper.py"
    # shell:
        # "fastqc {input} --outdir results/qc/fastqc_raw/ {params} > {log} 2>&1"

rule fastqc_trimmomatic_se:
    input:
        "results/trimmed/{sample}.trimmed.fastq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}.html",
        zip="results/qc/fastqc_trimmed/{sample}_fastqc.zip",
    params:
        "--quiet"
    threads:
        config['main_config']['threads']
    # conda: 
    #     "../envs/qc_trimming.yml"        
    wrapper:
        "v1.8.0/bio/fastqc"
    # script:
    #     "../scripts/fastqc_wrapper.py"
    # shell:
    #     "fastqc {input} --outdir results/qc/fastqc_trimmed/ {params} > {log} 2>&1"

rule fastqc_trimmomatic_pe_r1:
    input:
        "results/trimmed/{sample}_R1.trimmed.fastq.gz"
    output:
        html="results/qc/fastqc_trimmed/{sample}_R1.html",
        zip="results/qc/fastqc_trimmed/{sample}_R1_fastqc.zip",
    params:
        "--quiet"
    log:
        "logs/fastqc_trimmed/{sample}_R1.log",
    threads:
        config['main_config']['threads']
    # conda: 
    #     "../envs/qc_trimming.yml"        
    wrapper:
        "v1.8.0/bio/fastqc"
    # script:
    #     "../scripts/fastqc_wrapper.py"
    # shell:
    #     "fastqc {input} --outdir results/qc/fastqc_trimmed/ {params} > {log} 2>&1"

rule fastqc_trimmomatic_pe_r2:
    input:
        "results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        html="results/qc/fastqc_trimmed/{sample}_R2.html",
        zip="results/qc/fastqc_trimmed/{sample}_R2_fastqc.zip",
    params:
        "--quiet"
    log:
        "logs/fastqc_trimmed/{sample}_R2.log",
    threads:
        config['main_config']['threads']
    # conda: 
    #     "../envs/qc_trimming.yml"
    wrapper:
        "v1.8.0/bio/fastqc"
    # script:
    #     "../scripts/fastqc_wrapper.py"
    # shell:
    #     "fastqc {input} --outdir results/qc/fastqc_trimmed/ {params} > {log} 2>&1"