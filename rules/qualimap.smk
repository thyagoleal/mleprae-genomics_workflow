rule qualimap:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        outdir=directory("results/qc/qualimap/all_reads/{sample}"),
        html="results/qc/qualimap/all_reads/{sample}/qualimapReport.html",
    log:
        "logs/qualimap/{sample}.all-reads.log"
    threads:
        config['main_config']['threads']
    params:
        seq_protocol=config['params']['qualimap']['strand']
    # conda:
    #     "../envs/qualimap.yml"    
    shell:
        "qualimap bamqc -bam {input} -nt {threads} -outfile {output.html} "
        "-outdir {output.outdir} -p {params.seq_protocol} -outformat HTML 2> {log}"

rule qualimap_dedup:
    input:
        "results/dedup/{sample}.deduped.sorted.bam",
    output:
        outdir=directory("results/qc/qualimap/deduped/{sample}"),
        html="results/qc/qualimap/deduped/{sample}/qualimapReport.html",
    log:
        "logs/qualimap/{sample}.deduped.log"
    threads:
        config['main_config']['threads']
    params:
        seq_protocol=config['params']['qualimap']['strand']
    # conda:
    #     "../envs/qualimap.yml"    
    shell:
        "qualimap bamqc -bam {input} -nt {threads} -outfile {output.html} "
        "-outdir {output.outdir} -p {params.seq_protocol} -outformat HTML 2> {log}"
