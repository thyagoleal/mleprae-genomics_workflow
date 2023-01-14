if config['merge_prepare_philogeny_script']:
    rule prepare_phylogeny:
        input:
            expand("results/variant_calls/{sample}.varscan.vcf.gz", sample=good_coverage_samples())
        output:
                "results/phyloprep_output/output_MERGED_only-informative-sites_plus_LPM.fasta",
                "results/phyloprep_output/output_MERGED_only-informative-sites.txt",
                "results/phyloprep_output/output_MERGED_gaps-per-genome.txt",
                "results/phyloprep_output/output_MERGED_all-positions.txt.gz",
        params:
            ref_fa=config['ref']['ref-file_path'],
            mlepromatosis_vcf=config['ref']['mlepromatosis_vcf'],
            filter_bad=config['filters']['remove_bad'],
            filter_out=config['filters']['bad-sites'],
            merge_toggle=config['merge_prepare_philogeny_script'], # If False it will ignore previous_genomes
            previous_genomes=config['ref']['previous_genomes_phylogeny'],
            sample_names=good_coverage_samples(), # do not change this
            coverage_toggle=config['only_phylogeny_good_coverage'], 
            coverage_cutoff=config['coverage_threshold_phylogeny'],
        log:
            "logs/prepare_phylogeny/phylogeny-script.log"
        # conda:
        #     "../envs/basic.yml"
        script:
            "../scripts/prepare4phylogeny.py"

if not config['merge_prepare_philogeny_script']:
    rule prepare_phylogeny:
        input:
            expand("results/variant_calls/{sample}.varscan.vcf.gz", sample=good_coverage_samples())
        output:
                "results/phyloprep_output/output_only-informative-sites_plus_LPM.fasta",
                "results/phyloprep_output/output_only-informative-sites.txt",
                "results/phyloprep_output/output_gaps-per-genome.txt",
                "results/phyloprep_output/output_all-positions.txt.gz",
        params:
            ref_fa=config['ref']['ref-file_path'],
            mlepromatosis_vcf=config['ref']['mlepromatosis_vcf'],
            filter_out=config['filters']['bad-sites'],
            merge_toggle=config['merge_prepare_philogeny_script'], # If False it will ignore previous_genomes
            previous_genomes=config['ref']['previous_genomes_phylogeny'],
            sample_names=good_coverage_samples(), # do not change this
            coverage_toggle=config['only_phylogeny_good_coverage'], 
            coverage_cutoff=config['coverage_threshold_phylogeny'],
        log:
            "logs/prepare_phylogeny/phylogeny-script.log"
        # conda:
        #     "../envs/basic.yml"
        script:
            "../scripts/prepare4phylogeny.py"            

if config['merge_prepare_philogeny_script']:
    rule snpEff_phylogeny:
        input:
            infile="results/phyloprep_output/output_MERGED_only-informative-sites.txt",
            script="scripts/SNP-table_to_SNPeff-table_v3.sh",
        output:
            "results/phyloprep_output/output_MERGED_only-informative-sites_snpEff.txt"
        log:
            "logs/prepare_phylogeny/snpeff.log"
        # conda:
        #     "../envs/snpeff.yml"
        shell:
            "{input.script} {input.infile} > {output} 2>> {log} && " 
            "rm snpEff_summary.html snpEff_genes.txt"

if not config['merge_prepare_philogeny_script']:
    rule snpEff_phylogeny:
        input:
            infile="results/phyloprep_output/output_only-informative-sites.txt",
            script="scripts/SNP-table_to_SNPeff-table_v3.sh",
        output:
            "results/phyloprep_output/output_only-informative-sites_snpEff.txt"
        log:
            "logs/prepare_phylogeny/snpeff.log"
        # conda:
        #     "../envs/snpeff.yml"
        shell:
            "{input.script} {input.infile} > {output} 2>> {log} && " 
            "rm snpEff_summary.html snpEff_genes.txt"
