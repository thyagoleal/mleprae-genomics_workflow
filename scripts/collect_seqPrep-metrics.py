# Modify SeqPrep logs and create new log file

import re
from pathlib import Path
import os

file = ['logs/SeqPrep/{sample}.log'.format(sample = x) for x in snakemake.params.sample]
# logs = list(Path('.').glob(file))
logs = file

for log in logs:

    table = []
    logfile = open(log, 'r+')
    
    for line in logfile:
        table += [line]
    
    logfile.close()
    
    table = table[1:]
    pairs_processed = table[0]
    pairs_merged = table[1]

    regex = r'[0-9]{1,}'

    match_merged = re.search(regex, pairs_merged)
    results_merged = float(match_merged.group())

    match_processed = re.search(regex, pairs_processed)
    results_processed = float(match_processed.group())

    if results_processed == 0 or results_merged == 0:
        total_merged = 'Error! No data found.'
    else:
        total_merged = round(results_merged/results_processed * 100, 2)

    with open(log + '.clean', "w+") as output:
        for i in table:
            output.write(i)
        output.write('Merging is: ' + str(total_merged) + '% \n')

    logfile.close()

exit(0)

