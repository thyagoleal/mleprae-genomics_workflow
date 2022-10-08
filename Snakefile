# END USER, do not change this file.

include: "rules/helpers.smk"
configfile: "config/config.yaml"
singularity: "docker://condaforge/mambaforge:4.14.0-0"

# Python3 code -----------------------------------------------------------------

import pandas as pd
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

class bcolors:
    OKGREEN = "\033[92m"
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"

# get all sample names for PE and SE layouts
all_samples = units['sample'].to_list()
samples_pe = units.loc[units['design'] == 'paired-end', 'sample'].to_list()
samples_se = units.loc[units['design'] == 'single-end', 'sample'].to_list()

print("\n=======================================================================")

print(bcolors.HEADER +
    '>>> Processing {} samples in SE format:\n'.format(len(samples_se)) +
    bcolors.ENDC)

print(bcolors.OKGREEN + "\n".join(samples_se) + '\n' + bcolors.ENDC)

print(bcolors.OKBLUE +
    ">>> Processing {} samples in PE format:\n".format(len(samples_pe)) +
    bcolors.ENDC)

print(bcolors.OKGREEN + "\n".join(samples_pe) + bcolors.ENDC)

print("=======================================================================")

# Target rules -----------------------------------------------------------------

phyloprep_output = [
        "results/phyloprep_output/output_MERGED_gaps-per-genome.txt",
        "results/phyloprep_output/output_MERGED_all-positions.txt.gz",
        "results/phyloprep_output/output_MERGED_only-informative-sites.txt",
        "results/phyloprep_output/output_MERGED_only-informative-sites_plus_LPM.fasta",
        "results/phyloprep_output/output_MERGED_only-informative-sites_snpEff.txt"
        ]

if not config['merge_prepare_philogeny_script']:
    phyloprep_output = [i.replace("MERGED", "") for i in phyloprep_output]
          

rule all:
    input:
        ["results/qc/multiqc/multiqc.html",
        "results/pipeline_output_recap.tsv"] + phyloprep_output
        

include: "rules/fastqc.smk"
include: "rules/trimming.smk"
include: "rules/merge_overlapping.smk"
include: "rules/mapping.smk"
include: "rules/markduplicates.smk"
include: "rules/variant_calling.smk"
include: "rules/qualimap.smk"
include: "rules/vcf2table.smk"
include: "rules/multiqc.smk"
include: "rules/prepare_phylogeny.smk"
include: "rules/workflow_recap.smk"