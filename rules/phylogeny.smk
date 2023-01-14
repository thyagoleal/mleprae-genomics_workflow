rule max_parsimony_tree:
    input: 
        fasta=get_phyloprep_input(),
    output:
        "results/max-parsimony-tree/max-parsimony-tree.nwk", 
    params:
        mao_script="scripts/max-parsimony-tree.mao",
        outdir="results/max-parsimony-tree/max-parsimony-tree",
    log:
        "logs/mega11/max-parsimony-tree.log"
    shell:
        "megacc -a {params.mao_script} -d {input.fasta} -o {params.outdir} > {log} 2>&1"
