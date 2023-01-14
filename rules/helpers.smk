import pandas as pd
import snakemake.io
import os.path
from platform import system as system_name

# from snakemake.utils import validate
# from snakemake.utils import min_version

configfile: "config/config.yaml"

# validate(config, schema="..schemas/config.schema.yaml")
# validate(samples, schema="schemas/samples.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample"], drop=False)

# units.index = units.index.set_levels(
#     [i.astype(str) for i in units.index.levels]
#     ) # enforce str in index

# validate(units, schema="schemas/units.schema.yaml")

#### wildcard constraints ####

wildcard_constraints:
    sample="|".join(samples.index)

# Helper Functions =============================================================

def get_phyloprep_input():
    if config['merge_prepare_philogeny_script']:
        return ["results/phyloprep_output/output_MERGED_only-informative-sites_plus_LPM.fasta"]
    if not config['merge_prepare_philogeny_script']:
        return ["results/phyloprep_output/output_only-informative-sites_plus_LPM.fasta"]
    else: 
        raise ValueError("Could not determine phyloprep output")

def get_platform():
    """Returns the platform name as one string."""
    sys = system_name()
    if sys:
        return sys
    else:
        raise ValueError("Could not determine system type")

def good_coverage_samples():
    import csv
    from re import sub
    out = []
    
    if config['only_phylogeny_good_coverage'] and config['coverage_threshold_phylogeny'] is not None:
        with open('results/pipeline_output_recap.tsv') as file:
            tsv_file = csv.reader(file, delimiter="\t")
            for line in tsv_file:
                if float(sub('[A-Za-z\=]', '', line[4])) >= float(config['coverage_threshold_phylogeny']):
                    out.append(sub('\-vs\-.*', '', line[0]))
        return out
    if not config['only_phylogeny_good_coverage']:
        return units['sample'].to_list()

    else:
        raise ValueError("Unknown error in good_coverage_samples()")

def is_single_end(sample):
    """Return True if sample is single-end."""
    return pd.isnull(units.loc[sample, "fq2"])

def get_multiqc_input(wildcards):
    multiqc_input = []
    multiqc_input.extend(
        expand(
            ["results/qc/dedup/{sample}.metrics.txt",
            "results/qc/samtools-stats/{sample}.txt",
            "results/qc/samtools-flagstat/{sample}.txt",
            "logs/trimmomatic/{sample}.log",
            "logs/trimmomatic/{sample}.log",
            # 'results/qc/fastqc_raw/{sample}_fastqc.zip',
            # 'results/qc/fastqc_raw/{sample}_R1_fastqc.zip',
            # 'results/qc/fastqc_raw/{sample}_R2_fastqc.zip',
            'results/qc/fastqc_trimmed/{sample}_fastqc.zip',
            'results/qc/fastqc_trimmed/{sample}_R1_fastqc.zip',
            'results/qc/fastqc_trimmed/{sample}_R2_fastqc.zip',
            "results/variant_calls/snpeff/{sample}_SNPeff.csv",
            "results/qc/qualimap/{sample}/qualimapReport.html"],
            sample=units['sample'].to_list()
        )
    )

    return multiqc_input

def get_fastq(wildcards):
    """Get fastq files of given sample"""
    fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_fastq_se(wildcards):
    """Get fastq files for SE sample-unit pair."""
    fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) < 2:
        return {"r1": fastqs.fq1}
    return {}

def get_fastq_pe(wildcards):
    """Get fastq files for PE sample pairs."""
    if not is_single_end(**wildcards):
        fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
        if len(fastqs) == 2:
            return {"r1": fastqs.fq1, "r2": fastqs.fq2}

def get_trimming_input_se(wildcards):
    input = []
    tmp = units.loc[units['design'] == 'single-end', 'sample'].to_list()
    input.extend(
        expand(
            "data/samples/prepared-fastqs/{sample}.fastq.gz",
            sample=tmp)
            )
    return input

def get_trimming_input_pe(wildcards):
    tmp = units.loc[units['design'] == 'paired-end', 'sample']
    input = {
        'r1': expand("data/samples/prepared-fastqs/{sample}_R1.fastq.gz", sample=tmp),
        'r2': expand("data/samples/prepared-fastqs/{sample}_R2.fastq.gz", sample=tmp)
            }
    return input


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"'@RG\tID:{sample}\tSM:{sample}'".format(
        sample=wildcards.sample)

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}_R{read}.trimmed.fastq.gz",
            read=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}.trimmed.fastq.gz".format(**wildcards)

def get_trimmed_reads_2(wildcards):
    """Get fastq files for sample-unit pair."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": 'results/trimmed/{sample}_R1.trimmed.fastq.gz',
                "r2": 'results/trimmed/{sample}_R2.trimmed.fastq.gz'}
    return {"r1": 'results/trimmed/{sample}.trimmed.fastq.gz'}

def get_fastqc_results(wildcards):
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        # paired-end sample
        return expand(
            "results/qc/fastqc_raw/{sample}_R{read}.fastq.gz",
            read=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/qc/fastqc_raw/{sample}.fastq.gz".format(**wildcards)

def get_single_end_raw(wildcards):
    if is_single_end(**wildcards):
        return expand("data/samples/{sample}.fastq.gz",
            **wildcards
        )

def get_paired_end_raw(wildcards):
    if not is_single_end(**wildcards):
        return expand("data/samples/{sample}_{group}.fastq.gz",
            group=['R1', 'R2'],
            **wildcards
        )

def get_single_end_trimmed(wildcards):
    if is_single_end(**wildcards):
        return "results/trimmed/{sample}.trimmed.fastq.gz".format(**wildcards)

def get_all_bams(wildcards):
    if is_single_end(**wildcards):
        return "results/mapped/{sample}.SE.sorted.bam".format(**wildcards)
    return "results/mapped/{sample}.PE.sorted.bam".format(**wildcards)

def get_fastqc_output(type='fastqc_raw'):
    """Get fastqs before fastqc and return their future paths"""
    sample = snakemake.io.glob_wildcards("data/samples/"+"{sample}.fastq.gz")[0]
    out1 = []
    for i,j in enumerate(sample):
        ext = ".fastq.gz"
        if os.path.isfile("data/samples/" + j + ext):
            out1.append("results/qc/" + type + "/" + j + "_fastqc.zip")
    return out1
