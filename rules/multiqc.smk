rule multiqc:
    input:
        expand('results/qc/{dir}/{sample}_{read}_fastqc.zip', sample=samples_pe,
        dir=['fastqc_raw', 'fastqc_trimmed'], read=['R1', 'R2']),

        expand("logs/SeqPrep/{sample}.log.clean", sample=samples_pe),

        expand('results/qc/{dir}/{sample}_fastqc.zip', sample=samples_se,
        dir=['fastqc_raw', 'fastqc_trimmed']),

        expand(
        [
            "logs/trimmomatic/{sample}.log",
            "results/qc/dedup/{sample}.metrics.txt",
            "results/qc/samtools-stats/{sample}.txt",
            "results/qc/samtools-flagstat/{sample}.txt",
            "results/variant_calls/snpeff/{sample}_SNPeff.csv",
            "results/qc/qualimap/deduped/{sample}/qualimapReport.html",
            "results/qc/qualimap/all_reads/{sample}/qualimapReport.html", 
            "results/variant_calls/{sample}.varscan.vcf.gz.tbi",
            "results/variant_calls/snpeff/{sample}.varscan_snpeff.vcf.gz.tbi",
            "results/variant_calls/variants_bed/{sample}.{type}.bed.gz",
        ],
        sample=all_samples, type=['snps','insertions', 'deletions']),
    output:
        "results/qc/multiqc/multiqc.html"
    params:
        extra="--dirs"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v1.7.0/bio/multiqc"
    # script:
    #     "../scripts/multiqc_wrapper.py"
    # shell:
    #     "multiqc {params.extra} {input} -o results/qc/multiqc/ -f 2> {log}"