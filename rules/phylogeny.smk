rule max_parsimony_tree:
    input: 
        fasta=get_phyloprep_output(),
        mao_script="scripts/max-parsimony-tree.mao"
    output:
        directory("results/max-parsimony-tree")
    shell:
        "megacc -a {input.mao_script} -d {input.fasta} -o {output}"
