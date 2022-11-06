rule merge_overlapping_reads_pe:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
    output:
        out1="results/merged_reads/{sample}_R1.fastq.gz",
        out2="results/merged_reads/{sample}_R2.fastq.gz",
        out3="results/merged_reads/{sample}_R1.NotMerged.fastq.gz",
        out4="results/merged_reads/{sample}_R2.NotMerged.fastq.gz",
        out5="results/merged_reads/{sample}.merged.fastq.gz",
        out6="results/merged_reads/{sample}.alignments.txt.gz",

    log: "logs/SeqPrep/{sample}.log"
    params:
        forward_primer="AGATCGGAAGAGCACACGTCT",
        reverse_primer="AGATCGGAAGAGCGTCGTGTA",
    # conda:
    #     "../envs/seqprep.yml"    
    shell:
        "SeqPrep -f {input.r1} -r {input.r2} -1 {output.out1} -2 {output.out2} "
        "-3 {output.out3} -4 {output.out4} -A {params.forward_primer} "
        "-B {params.reverse_primer} -L 40 -o 40 -y I -s {output.out5} -E "
        "{output.out6} 2> {log}"
