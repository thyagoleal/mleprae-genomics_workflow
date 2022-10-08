# Get statistics from run

from datetime import date

samples = set(snakemake.params.samples)

for sample in samples:

    # bowtie
    ifile_mapping = open('logs/bowtie2/' + sample + '.log', 'r')
    table_mapping = []
    for i in ifile_mapping:
        table_mapping += [i]
    mapping_percentage = table_mapping[-1].rstrip()

    # trimmomatic
    ifile_trimming = open('logs/trimmomatic/' + sample + '.log', 'r')
    table_trimming = []
    for i in ifile_trimming:
        table_trimming += [i]
    trimmomatic_info = table_trimming[-2].rstrip()

    # bamqc deduped data
    ifile_bamqc = open('results/qc/qualimap/deduped/' + sample + '/genome_results.txt', 'r')
    table_bamqc = []
    for i in ifile_bamqc:
        table_bamqc += [i]
    mean_cov = str()
    stdev_cov = str()
    for i in table_bamqc:
        if 'mean coverageData' in i:
            mean_cov = i.rstrip().strip('\t')
        if 'std coverageData' in i:
            stdev_cov = i.rstrip()

    # bamqc all reads
    ifile_bamqc_nodup = open('results/qc/qualimap/all_reads/' + sample + '/genome_results.txt', 'r')
    table_bamqc_nodup = []
    for j in ifile_bamqc_nodup:
        table_bamqc_nodup += [j]
    mean_cov_nodup = str()
    stdev_cov_nodup = str()
    for i in table_bamqc_nodup:
        if 'mean coverageData' in i:
            mean_cov_nodup = i.rstrip().strip('\t')
        if 'std coverageData' in i:
            stdev_cov_nodup = i.rstrip()

    ofile = open('results/pipeline_output_recap.tsv', 'a') # append to file

    ofile.write(sample + '-vs-' + snakemake.config["ref"]["name"] +
     '\t' + str(date.today()) + '\t' + mapping_percentage + '\t' + 
     trimmomatic_info + '\t' + mean_cov + '\t' + stdev_cov + '\t'+ 
     'without duplicated reads:' + mean_cov_nodup + '\t'+ 
     'without duplicated reads:' + stdev_cov_nodup+ '\n')

    ofile.close()