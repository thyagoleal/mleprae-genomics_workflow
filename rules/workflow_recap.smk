rule workflow_recap_metrics:
    input:
        expand(["results/qc/qualimap/all_reads/{sample}/",
        "results/qc/qualimap/deduped/{sample}/"], 
        sample=all_samples),
    output:
        "results/pipeline_output_recap.tsv"
    params:
        samples=all_samples,
    # conda:
    #     "../envs/basic.yml"
    script:
        "../scripts/collect_metrics.py"
    
rule collect_seqprep_metrics:
    input:
        expand("logs/SeqPrep/{sample}.log", sample=samples_pe),
    output:
        expand("logs/SeqPrep/{sample}.log.clean", sample=samples_pe),
    log:
        "logs/collect_seqprep/collect_seqprep_metrics.log", 
    params:
        sample=samples_pe,    
    # conda:
    #     "../envs/basic.yml"
    script:
        "../scripts/collect_seqPrep-metrics.py"