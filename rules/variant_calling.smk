rule variant_calling:
    input:
        bam="results/dedup/{sample}.deduped.sorted.bam",
        bai="results/dedup/{sample}.deduped.sorted.bai",
    output:
        temp("results/variant_calls/{sample}.vcf")
    log:
        "logs/samtools_mpileup/{sample}.log"
    params:
        main=config['params']['samtools-mpileup'],
        ref=config['ref']['ref-file_path'],
    threads:
        config['main_config']['threads']
    # conda:
    #     "../envs/mapping_tools.yml"    
    shell:
        "samtools mpileup {params.main} -f {params.ref} {input.bam} 2> {log} "
        "-o {output}"

rule varscan:
    input:
        "results/variant_calls/{sample}.vcf"
    output:
        "results/variant_calls/{sample}.varscan.vcf.gz"
    log:
        log1="logs/varscan/{sample}.log",
        log2="logs/bgzip/{sample}.varscan.log",
    threads:
        config['main_config']['threads']
    params:
        config['params']['varscan']
    # conda:
    #     "../envs/varscan.yml"    
    shell:
        "varscan mpileup2cns {input} {params} --output-vcf 1 2>> {log.log1} | "
        "sed -r -e 's/Sample1/{wildcards.sample}/' -e "
        " 's/^[^#\t]*\t/Chromosome\t/g' | bgzip -@ {threads} > {output} "
        "2> {log.log2}"

rule snpeff:
    input:
        "results/variant_calls/{sample}.varscan.vcf.gz",
    output:
        html="results/variant_calls/snpeff/{sample}_SNPeff-Summary.html",
        vcf="results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz",
        csv="results/variant_calls/snpeff/{sample}_SNPeff.csv",
    log:
        log1="logs/snpEff/{sample}.log",
        log2="logs/bgzip/{sample}.snpeff.log",
    threads:
        config['main_config']['threads']
    params:
        misc=config['params']['snpEff']['misc'],
        db=config['params']['snpEff']['genome_database'],
    # conda:
    #     "../envs/snpeff.yml"    
    shell:
        "snpEff eff {params.db} {input} {params.misc} -htmlStats {output.html} "
        "-csvStats {output.csv} 2> {log.log1} | bgzip -@ {threads} > {output.vcf} 2> {log.log2}"

rule tabix:
    input:
        in1="results/variant_calls/{sample}.varscan.vcf.gz",
        in2="results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz",
    output:
        out1="results/variant_calls/{sample}.varscan.vcf.gz.tbi",
        out2="results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz.tbi",
    threads:
        config['main_config']['threads']
    log:
        "logs/tabix/{sample}.txt"
    # conda:
    #     "../envs/varscan.yml"
    shell:
         "tabix -f -p vcf --verbosity 9 {input.in1} 2>> {log} && "
         "tabix -f -p vcf --verbosity 9 {input.in2} 2>> {log}"
