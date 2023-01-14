# adapted from Chloeloiiseau's 2_Phylogeny_v7.property
# Converted from python2 to python3
# gzip now reads rt instead of rb

import fileinput
import gzip
import os, sys, subprocess, getopt, re, glob, math, shutil
from subprocess import call
from pathlib import Path
from pandas import read_table


class bcolors:
    OKGREEN = "\033[92m"
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    FAIL = '\033[91m'

# import from snakemake rule 

print(bcolors.HEADER +
    '\n>>> Running Philogeny preparation script...\n' + bcolors.ENDC)

MERGE_TOGGLE = snakemake.params.merge_toggle
REF = snakemake.params.ref_fa
BASENAMES = snakemake.params.sample_names # list of sample names
INPUT_VCFS = snakemake.input # space delimited list of vcf files
LPM = snakemake.params.mlepromatosis_vcf
FILTER_TOGGLE = snakemake.params.filter_bad # if False, next is ignored
FILTER = snakemake.params.filter_out
MERGE = snakemake.params.previous_genomes
OUTDIR = Path(snakemake.output[0]).parent.absolute()
OUTPUT = "output"
OUTPUT = os.path.join(OUTDIR, OUTPUT)

# convert list to space separated string
# BASENAMES_str = ' '.join([str(name) for name in BASENAMES])

# separate string to list
INPUT_VCFS_str = INPUT_VCFS
INPUT_VCFS = ' '.join(INPUT_VCFS).split()


if len(INPUT_VCFS) != len(BASENAMES):
    print(bcolors.FAIL + '\nERROR: number of entries in basenames and VCFs are not equal!\n' + bcolors.ENDC)

# the script will run not MERGED, default is to merge
if not MERGE_TOGGLE:
    ''' First step: Make a table of all positions '''
    FINAL_LIST = [list(range(0,3268204))] # FINAL_LIST is a list which will be appended gradually with each genome inputed. First we add the coordinates for each nucleotide. Ex: [[coordinates 0 to 3268204],[genome1],[genome2],etc]
    # add to TN to FINAL_LIST
    ifile  = open(REF,'r') # open TN reference fasta file
    next(ifile) # skip the first line of the file which is the header
    reference_list = [] # create an empty list which will contain sublists and where each sublist contains one line from the fasta
    for i in ifile: # loop into the TN fasta file
        reference_list += i.rstrip() # rstrip removes all whitespace characters contained at the end of the string 'i'. We do that because the fasta goes to the line every 70 char
    TN = [] # Make a new empty list
    for i in reference_list:
        TN+=i
    TN.insert(0,'TN')
    FINAL_LIST += [TN]

    for basename,vcf in zip(list(range(len(BASENAMES))),list(range(len(INPUT_VCFS)))):
        print(BASENAMES[basename], INPUT_VCFS[vcf])
        vcf_file = gzip.open(INPUT_VCFS[vcf],'rt') # open vcf file
        table_vcf = []
        for i in vcf_file:
            if i[0] != '#': # do not take headers of vcf
                table_vcf += [i.split('\t')] # this list contains sublists where each sublists is a row in the vcf file


        for cat in range(len(table_vcf)):
            table_vcf[cat][9] = table_vcf[cat][9].split(":") # split the 10th column of the vcf according to the ':'. i.e 1/1:3:1:1:0:1:14,29%:5E-1:0:37:0:0:0:1
            table_vcf[cat][7] = table_vcf[cat][7].split(";")

            if len(table_vcf[cat][9])>=7: # if there is more than 7 elements in this column:
                table_vcf[cat][9][6]=table_vcf[cat][9][6].replace(",",".") # replace the frequencies. Example: 14,29 --> 14.29 (accepted by python)
                table_vcf[cat][9][7] = table_vcf[cat][9][7].replace(",",".") # replace the pvalue
                table_vcf[cat][9][7] = float(table_vcf[cat][9][7])




        fasta_list = []
        pos = 1
        for i in range(1,len(table_vcf)): #i corresponds to each row of the vcf file
            a = int(table_vcf[i][1]) # the variable 'a' will be the position in the vcf file (so the 2nd column of the vcf).
            b = pos
            index_equalsign_ADP = table_vcf[i][7][0].index('=') # index of the equal sign in ADP=n
            index_equalsign_NC = table_vcf[i][7][4].index('=')  # index of the equal sign in NC=n

            if a == b:

                if table_vcf[i][4] == ".": # REF
                    # if coverage above 5:
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if len(table_vcf[i][9])>=7 and 20.0<=float(table_vcf[i][9][6][:-1])<=80.0: # MIX SITUATION: reference frequency is between 20% and 80%
                            fasta_list += ['-']
                        elif len(table_vcf[i][9])>=7 and float(table_vcf[i][9][6][:-1])<=20.0: # clear REF base: variant frequency below 20%
                            if len(table_vcf[i][3])>1: # deletion
                                fasta_list +=['D']
                            elif len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                            else:
                                fasta_list += ['-']
                        else:
                            fasta_list += ['-']
                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']



                elif table_vcf[i][4] != ".": # ALT
                    # if true variant (i.e coverage at least 5, variant frequency at least 80%, min supporting reads 3, quality of base above 15 and p value below 0.05)
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if table_vcf[i][9][7]<=float(0.05) and float(table_vcf[i][9][6][:-1])>=float(80) and int(table_vcf[i][9][5])>=3 and int(table_vcf[i][9][9])>15:
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1:
                                fasta_list += [table_vcf[i][4]]
                        elif float(table_vcf[i][9][6][:-1])<float(20): # if variant frequency below 20% (i.e reference)
                            if len(table_vcf[i][4])>1 and len(table_vcf[i][3])==1: # insertion
                                fasta_list += [table_vcf[i][3]] # write the reference
                            elif len(table_vcf[i][3])>1 and len(table_vcf[i][4])==1: # deletion
                               fasta_list += [table_vcf[i][4]] # write the alt as it corresponds to the reference
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]] # write the reference
                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])>=3: # if variant frequency between 20% and 80% (i.e mix situation) and min supporting reads 3.
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += ['-']
                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])<3:
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += [table_vcf[i][3]]
                        else:
                            fasta_list += ['-']


                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']

            elif a!=b:
                pos = int(a)
                fasta_list += list('-'*(pos-b))

                if table_vcf[i][4] == ".": # REF
                    # if coverage above 5:
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if len(table_vcf[i][9])>=7 and 20.0<=float(table_vcf[i][9][6][:-1])<=80.0: # MIX SITUATION: reference frequency is between 20% and 80%
                            fasta_list += ['-']
                        elif len(table_vcf[i][9])>=7 and float(table_vcf[i][9][6][:-1])<=20.0: # clear REF base: variant frequency below 20%
                            if len(table_vcf[i][3])>1: # Deletion
                                fasta_list +=['D']
                            elif len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                        else:
                            fasta_list += ['-']
                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']



                elif table_vcf[i][4] != ".": # ALT
                    # if true variant (i.e coverage at least 5, variant frequency at least 80%, min supporting reads 3, quality of base above 15 and p value below 0.05)
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if table_vcf[i][9][7]<=float(0.05) and float(table_vcf[i][9][6][:-1])>=float(80) and int(table_vcf[i][9][5])>=3 and int(table_vcf[i][9][9])>15:
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1:  # Deletion
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1:
                                fasta_list += [table_vcf[i][4]]
                        elif float(table_vcf[i][9][6][:-1])<float(20): # if variant frequency below 20% (i.e reference)
                            if len(table_vcf[i][4])>1 and len(table_vcf[i][3])==1: # insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1 and len(table_vcf[i][4])==1: # deletion
                               fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])>=3: # if variant frequency between 20% and 80% (i.e mix situation) and min supporting reads 3.
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += ['-']
                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])<3:
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += [table_vcf[i][3]]
                        else:
                            fasta_list += ['-']

                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']

            pos += 1

        extension_list = []
        if len(fasta_list) < 3268203:
            extension_list += list('-'*(3268203-len(fasta_list)))
        fasta_list.extend(extension_list)
        fasta_list.insert(0,BASENAMES[basename])
        FINAL_LIST += [fasta_list]



    # write in columns
    with open(OUTPUT+'_all-positions.txt','w') as f:
        for x in zip(*FINAL_LIST):
            strs = '\t'.join(str(value) for value in x)
            f.write(strs+'\n')
        f.close()

    if FILTER_TOGGLE:
        
        ''' Second Step: Edit the file created in the First step, to remove  the repetitive regions, rRNAs and uninformative sites -> create a new file with only informative positions!'''
        input_file_filter = open(FILTER,'r')
        # filter table is a list of positions to remove (rRNAS and repeats positions)

        filter_table = []
        for posi in input_file_filter:
            filter_table += [posi.strip()]

        # pos_dic is the same than filter_table except in dictionary as much faster to parse.

        pos_dic = {}
        for j in filter_table:
            pos_dic[j] = ''


        for line in fileinput.input(OUTPUT + '_all-positions.txt',inplace=True, backup='.bak'): # modify the file created in the first step
            index = line.index('\t')
            # Here I remove from the file the positions that are in repetitive regions of the genome or that are rRNAs
            if line[0] == '0':
                print(line, end=' ')
            else:
                if line[:index] not in pos_dic:
                #Here I remove uninformative sites : only A, only G, only T, only C, only - or mix of - and A, - and T, - and C, - and G
                    if not (('-' in line[index+1:] and 'A' not in line[index+1:] and 'C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:]) or ('A' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('C' in line[index+1:] and ('A' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('T' in line[index+1:] and ('C' not in line[index+1:] and 'A' not in line[index+1:] and 'G' not in line[index+1:])) or ('G' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'A' not in line[index+1:])) or ('D' in line[index+1:]) or ('I' in line[index+1:])):
                        print(line, end='')

        subprocess.call('mv '+OUTPUT+'_all-positions.txt '+ OUTPUT+'_only-informative-sites.txt',shell=True)
        subprocess.call('mv '+OUTPUT+'_all-positions.txt.bak '+ OUTPUT+'_all-positions.txt',shell=True)

        fileinput.close()

    if not FILTER_TOGGLE:
        
        for line in fileinput.input(OUTPUT + '_all-positions.txt',inplace=True, backup='.bak'): # modify the file created in the first step
            index = line.index('\t')
            # Here I remove from the file the positions that are in repetitive regions of the genome or that are rRNAs
            if line[0] == '0':
                print(line, end=' ')
            else:
                #Here I remove uninformative sites : only A, only G, only T, only C, only - or mix of - and A, - and T, - and C, - and G
                if not (('-' in line[index+1:] and 'A' not in line[index+1:] and 'C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:]) or ('A' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('C' in line[index+1:] and ('A' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('T' in line[index+1:] and ('C' not in line[index+1:] and 'A' not in line[index+1:] and 'G' not in line[index+1:])) or ('G' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'A' not in line[index+1:])) or ('D' in line[index+1:]) or ('I' in line[index+1:])):
                    print(line, end='')

        subprocess.call('mv '+OUTPUT+'_all-positions.txt '+ OUTPUT+'_only-informative-sites.txt',shell=True)
        subprocess.call('mv '+OUTPUT+'_all-positions.txt.bak '+ OUTPUT+'_all-positions.txt',shell=True)

        fileinput.close()


    ''' Third Step: Make a fasta from the new txt file created in the Second step --> fasta file with only informative positions'''
    # Ported to python3 

    # only could make it work using pandas (Thyago)
    df = read_table(OUTPUT+'_only-informative-sites.txt', header=None)
    informative_positions = list(tuple([tuple(df[col]) for col in df]))

    # print the length of informative positions to see if all have the same length
    print(bcolors.OKBLUE + '\nindex ;','name ;','number of informative positions\n' + bcolors.ENDC)

    for i in range(len(informative_positions)):
        print(i,';',informative_positions[i][0],';',len(informative_positions[i]))

    #create a fasta file
    ofile2 = open(OUTPUT+'_only-informative-sites_plus_LPM.fasta','w+')

    for j in range(1,len(informative_positions)): # for each column of the informative position tabular file (except the positions)
        ofile2.write('\n>'+informative_positions[j][0]+'\n') # write the first line of the fasta i.e >TN
        for nucleotide in range(1,len(informative_positions[j])):
            ofile2.write(informative_positions[j][nucleotide]) # write the base
    ofile2.close()
    ifile.close()

    ''' Fourth Step: select the LPM bases corresponding to each informative positions of M. leprae dataset --> append the fasta with LPM sequence'''
    print(bcolors.OKBLUE + '\nAdding LPM to fasta...\n' + bcolors.ENDC)

    informative_sites_only = list(informative_positions[0][1:]) # this is a list of the informative positions (ex: 73, 8637, 837272, 33992928)

    # lpm_file = open(LPM,'r') # open the LPM VCF file
    lpm_file = gzip.open(LPM,'rt') # Thyago

    # extract information of vcf in lists
    table = []
    for i in lpm_file:
        table += [i]
    for j in range(len(table)):
        table[j]=table[j].split('\t')
    newtable2 = []
    for k in table:
        if len(k)>2:
            newtable2 += [k]


    fasta_list = []
    LPM_pos_dic = {}

    for i in newtable2:
        LPM_pos_dic[i[1]] = ''

    for l in informative_sites_only:
            if str(l) in LPM_pos_dic:
                for i in newtable2:
                    if i[1] == str(l):
                        if i[4] == '.': # REF
                                fasta_list += [i[3]]
                        elif i[4] != '.': # ALT
                                fasta_list += [i[4]]
            elif str(l) not in LPM_pos_dic:
                fasta_list += ['-']

    ifile2 = open(OUTPUT+'_only-informative-sites_plus_LPM.fasta','a') # append (a) the fasta file with the LPM sequence corresponding to the M. leprae information sites
    ifile2.write('\n>LPM\n')
    for base in fasta_list:
      ifile2.write(base)
    ifile2.close()
    lpm_file.close()
    # cut every 70 characters using unix code.. i diddn't manage in python :-(
    # Andrej: I commented this out because I prefer to have 1 line/sample

    #subprocess.call("grep -oE '.{1,70}' "+ OUTPUT+"_only-informative-sites_plus_LPM.fasta > "+OUTPUT+'_2.fasta',shell=True)
    #subprocess.call("rm "+OUTPUT+"_only-informative-sites_plus_LPM.fasta", shell=True)
    #subprocess.call("mv "+OUTPUT+"_2.fasta "+OUTPUT+"_only-informative-sites_plus_LPM.fasta",shell=True )


    ''' Fifth Step : Get the number of gaps per genome '''
    print(bcolors.OKBLUE + '\nCounting gaps...\n' + bcolors.ENDC)
    # The code below breaks when fasta is not wrapped. I will include a bash command that will also calculate the % of gaps.

    #if args.gaps:
        #fasta=open(OUTPUT+"_only-informative-sites_plus_LPM.fasta",'r')
        #file=fasta.readlines()
        #newtable3=[]

        #for lines in file:
                #if lines[0]!='>':
                    #newtable3[-1]+=lines
                #else:
                    #newtable3+=[lines]


        #for i in range(len(newtable3)):
                #newtable3[i]=newtable3[i].split('\n')
                #newtable3[i][1:]=[''.join(newtable3[i][1:])]


        #outfile2 = open(OUTPUT+'_gaps-per-genome.txt','w')
        #n = 0
        #for i in range(len(newtable3)):
            #for j in newtable3[i][1]:
                #if j == '-':
                    #n += 1

            #outfile2.write(newtable3[i][0][1:] + '\t'+ str(n)+'\n')
            #n = 0


        #for base in range(len(newtable3)):
            #print newtable3[base][0][1:], len(newtable3[base][1])

        #fasta.close()
        #outfile2.close()

    # Need to use bash becasue the default sh will make Syntax error: "(" unexpected
    subprocess.call("cat <(echo -e 'Sample\\tTotal_gaps\\tGaps_percent') <(paste <(grep '>' "+OUTPUT+"_only-informative-sites_plus_LPM.fasta | tr -d '>') <(grep -A1 '>' "+OUTPUT+"_only-informative-sites_plus_LPM.fasta | sed '/>/d' | while read p; do paste <(echo $p | grep -o '-' | wc -l) <(echo \"scale=4; $(echo $p | grep -o '-' | wc -l)/$(echo $p | wc -m)\" | bc -l | sed 's/^\./0./');done)) > "+OUTPUT+"_gaps-per-genome.txt", shell=True, executable='/bin/bash')

    # Finally, since I removed fasta wrapping, the file starts with an empty line, which makes MEGA complain. Remove the 1st line from the fasta:
    # sed -i would be good, but sed behaves differently in Mac, so we will do it the safe way:
    print(bcolors.OKBLUE + '\nPolish fasta...\n' + bcolors.ENDC)
    subprocess.call("tail -n+2 "+OUTPUT+"_only-informative-sites_plus_LPM.fasta > temp.fasta &&mv temp.fasta "+OUTPUT+"_only-informative-sites_plus_LPM.fasta", shell=True)

    # Finally, gzip the *all-positions.txt
    print(bcolors.OKBLUE + '\nCompressing '+OUTPUT+'_all-positions.txt...\n' + bcolors.ENDC)
    subprocess.call('gzip '+OUTPUT+'_all-positions.txt',shell=True)
    
    print(bcolors.OKGREEN + "ALL DONE!\n" + bcolors.ENDC)
    exit(0)

if MERGE_TOGGLE:

    # subprocess.call("cp " + MERGE + './' + MERGE, shell=True)

    # this file will be updated at each step
    # MERGE = Path("MERGE").resolve()

    ''' First step: Make a table of all positions for current dataset'''
    FINAL_LIST = [list(range(0,3268204))] # FINAL_LIST is a list which will be appended gradually with each genome inputed. First we add the coordinates for each nucleotide. Ex: [[coordinates 0 to 3268204],[genome1],[genome2],etc]
    # add to TN to FINAL_LIST
    ifile  = open(REF,'r') # open TN reference fasta file
    next(ifile) # skip the first line of the file which is the header
    reference_list = [] # create an empty list which will contain sublists and where each sublist contains one line from the fasta
    for i in ifile: # loop into the TN fasta file
        reference_list += i.rstrip() # rstrip removes all whitespace characters contained at the end of the string 'i'. We do that because the fasta goes to the line every 70 char
    TN = [] # Make a new empty list
    for i in reference_list:
        TN+=i
    TN.insert(0,'TN')
    FINAL_LIST += [TN]

    ifile.close()

    for basename,vcf in zip(list(range(len(BASENAMES))),list(range(len(INPUT_VCFS)))):
        vcf_file = gzip.open(INPUT_VCFS[vcf],'rt') # open vcf file
        table_vcf = []
        for i in vcf_file:
            if i[0] != '#': # do not take headers of vcf
                table_vcf += [i.split('\t')] # this list contains sublists where each sublists is a row in the vcf file

        vcf_file.close()

        for cat in range(len(table_vcf)):
            table_vcf[cat][9] = table_vcf[cat][9].split(":") # split the 10th column of the vcf according to the ':'. i.e 1/1:3:1:1:0:1:14,29%:5E-1:0:37:0:0:0:1
            table_vcf[cat][7] = table_vcf[cat][7].split(";")

            if len(table_vcf[cat][9])>=7: # if there is more than 7 elements in this column:
                table_vcf[cat][9][6]=table_vcf[cat][9][6].replace(",",".") # replace the frequencies. Example: 14,29 --> 14.29 (accepted by python)
                table_vcf[cat][9][7] = table_vcf[cat][9][7].replace(",",".") # replace the pvalue
                table_vcf[cat][9][7] = float(table_vcf[cat][9][7])

        fasta_list = []
        pos = 1
        for i in range(1,len(table_vcf)): #i corresponds to each row of the vcf file
            a = int(table_vcf[i][1]) # the variable 'a' will be the position in the vcf file (so the 2nd column of the vcf).
            b = pos
            index_equalsign_ADP = table_vcf[i][7][0].index('=')
            index_equalsign_NC = table_vcf[i][7][4].index('=')

            if a == b:

                if table_vcf[i][4] == ".": # REF
                    # if coverage above 5:
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if len(table_vcf[i][9])>=7 and 20.0<=float(table_vcf[i][9][6][:-1])<=80.0: # MIX SITUATION: reference frequency is between 20% and 80%
                            fasta_list += ['-']
                        elif len(table_vcf[i][9])>=7 and float(table_vcf[i][9][6][:-1])<=20.0: # clear REF base: variant frequency below 20%
                            if len(table_vcf[i][3])>1: # del
                                fasta_list +=['D']
                            elif len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                        else:
                            fasta_list += ['-']
                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']

                elif table_vcf[i][4] != ".": # ALT
                    # if true variant (i.e coverage at least 5, variant frequency at least 80%, min supporting reads 3, quality of base above 15 and p value below 0.05)
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if table_vcf[i][9][7]<=float(0.05) and float(table_vcf[i][9][6][:-1])>=float(80) and int(table_vcf[i][9][5])>=3 and int(table_vcf[i][9][9])>15:
                            if len(table_vcf[i][4])>1:
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1:
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1:
                                fasta_list += [table_vcf[i][4]]
                        elif float(table_vcf[i][9][6][:-1])<float(20): # if variant frequency below 20% (i.e reference)
                            if len(table_vcf[i][4])>1 and len(table_vcf[i][3])==1: # insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1 and len(table_vcf[i][4])==1: # deletion
                               fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])>=3: # if variant frequency between 20% and 80% (i.e mix situation) and min supporting reads 3.
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += ['-']

                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])<3:
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += [table_vcf[i][3]]

                        else:
                            fasta_list += ['-']

                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']

            elif a!=b:
                pos = int(a)
                fasta_list += list('-'*(pos-b))

                if table_vcf[i][4] == ".": # REF
                    # if coverage above 5:
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if len(table_vcf[i][9])>=7 and 20.0<=float(table_vcf[i][9][6][:-1])<=80.0: # MIX SITUATION: reference frequency is between 20% and 80%
                            fasta_list += ['-']
                        elif len(table_vcf[i][9])>=7 and float(table_vcf[i][9][6][:-1])<=20.0: # clear REF base: variant frequency below 20%
                            if len(table_vcf[i][3])>1: # indel
                                fasta_list +=['D']
                            elif len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                        else:
                            fasta_list += ['-']
                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']

                elif table_vcf[i][4] != ".": # ALT
                    # if true variant (i.e coverage at least 5, variant frequency at least 80%, min supporting reads 3, quality of base above 15 and p value below 0.05)
                    if int(table_vcf[i][7][0][index_equalsign_ADP+1:]) >= 5:
                        if table_vcf[i][9][7]<=float(0.05) and float(table_vcf[i][9][6][:-1])>=float(80) and int(table_vcf[i][9][5])>=3 and int(table_vcf[i][9][9])>15:
                            if len(table_vcf[i][4])>1:
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1:
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1:
                                fasta_list += [table_vcf[i][4]]
                        elif float(table_vcf[i][9][6][:-1])<float(20): # if variant frequency below 20% (i.e reference)
                            if len(table_vcf[i][4])>1 and len(table_vcf[i][3])==1: # insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1 and len(table_vcf[i][4])==1: # deletion
                               fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1: # SNP
                                fasta_list += [table_vcf[i][3]]
                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])>=3: # if variant frequency between 20% and 80% (i.e mix situation) and min supporting reads 3.
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += ['I']
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += ['D']
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += ['-']

                        elif 20.0<=float(table_vcf[i][9][6][:-1])<=80.0 and int(table_vcf[i][9][5])<3:
                            if len(table_vcf[i][4])>1: # Insertion
                                fasta_list += [table_vcf[i][3]]
                            elif len(table_vcf[i][3])>1: # Deletion
                                fasta_list += [table_vcf[i][4]]
                            elif len(table_vcf[i][4])==1 and len(table_vcf[i][3])==1:
                                fasta_list += [table_vcf[i][3]]
                        else:
                            fasta_list += ['-']


                    elif int(table_vcf[i][7][0][index_equalsign_ADP+1:]) < 5:
                        fasta_list += ['-']

            pos += 1
        extension_list = []
        if len(fasta_list) < 3268203:
            extension_list += list('-'*(3268203-len(fasta_list)))
        fasta_list.extend(extension_list)
        fasta_list.insert(0,BASENAMES[basename])
        FINAL_LIST += [fasta_list]

    # write in columns
    with open(OUTPUT+'_all-positions.txt','w') as f:
        for x in zip(*FINAL_LIST):
            strs = '\t'.join(str(value) for value in x)
            f.write(strs+'\n')

    f.close()

    # if you give an existing tabular file with all positions and you wish to merge it to the new tabular file with all positions (in order to find informative positions)

    # subprocess.call('cut -f 3- '+OUTPUT+'_all-positions.txt > '+OUTPUT+'2_all-positions.txt',shell=True)
    # subprocess.call('paste -d "\t" '+MERGE+' '+OUTPUT+'2_all-positions.txt > ' + OUTPUT+'_MERGED_all-positions.txt',shell=True)
    # subprocess.call('rm '+OUTPUT+'2_all-positions.txt',shell=True)
    subprocess.call('paste <(gzip -dc '+MERGE+') <(cut -f 3- '+OUTPUT+'_all-positions.txt) > '+OUTPUT+'_MERGED_all-positions.txt',shell=True, executable='/bin/bash')

    if FILTER_TOGGLE:

        input_file_filter = open(FILTER,'r')
        filter_table = []
        for posi in input_file_filter:
            filter_table += [posi.strip()]
        # pos_dic is the same than filter_table except in dictionary as much faster to parse.

        input_file_filter.close()

        pos_dic = {}
        for j in filter_table:
            pos_dic[j] = ''

        for line in fileinput.input(OUTPUT+'_MERGED_all-positions.txt',inplace=True, backup='.bak'): # modify the file created in the first step
            index = line.index('\t')
                # Here I remove from the file the positions that are in repetitive regions of the genome or that are rRNAs
            if line[0] == '0':
                print(line, end=' ')
            else:
                if line[:index] not in pos_dic:
            #Here I remove uninformative sites : only A, only G, only T, only C, only - or mix of - and A, - and T, - and C, - and G
                    if not (('-' in line[index+1:] and 'A' not in line[index+1:] and 'C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:]) or ('A' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('C' in line[index+1:] and ('A' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('T' in line[index+1:] and ('C' not in line[index+1:] and 'A' not in line[index+1:] and 'G' not in line[index+1:])) or ('G' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'A' not in line[index+1:])) or ('D' in line[index+1:]) or ('I' in line[index+1:])):
                        print(line, end=' ')
        # rename files: _MERGED_all-positions.txt contains informative sites and _MERGED_all-positions.txt.bak contains all positions
        subprocess.call('mv '+OUTPUT+'_MERGED_all-positions.txt '+ OUTPUT+'_MERGED_only-informative-sites.txt',shell=True)
        subprocess.call('mv '+OUTPUT+'_MERGED_all-positions.txt.bak '+ OUTPUT+'_MERGED_all-positions.txt',shell=True)

        fileinput.close()
    
    if not FILTER_TOGGLE:

        for line in fileinput.input(OUTPUT+'_MERGED_all-positions.txt',inplace=True, backup='.bak'): # modify the file created in the first step
            index = line.index('\t')
                # Here I remove from the file the positions that are in repetitive regions of the genome or that are rRNAs
            if line[0] == '0':
                print(line, end=' ')
            else:
            #Here I remove uninformative sites : only A, only G, only T, only C, only - or mix of - and A, - and T, - and C, - and G
                if not (('-' in line[index+1:] and 'A' not in line[index+1:] and 'C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:]) or ('A' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('C' in line[index+1:] and ('A' not in line[index+1:] and 'T' not in line[index+1:] and 'G' not in line[index+1:])) or ('T' in line[index+1:] and ('C' not in line[index+1:] and 'A' not in line[index+1:] and 'G' not in line[index+1:])) or ('G' in line[index+1:] and ('C' not in line[index+1:] and 'T' not in line[index+1:] and 'A' not in line[index+1:])) or ('D' in line[index+1:]) or ('I' in line[index+1:])):
                    print(line, end=' ')
                    
        # rename files: _MERGED_all-positions.txt contains informative sites and _MERGED_all-positions.txt.bak contains all positions
        subprocess.call('mv '+OUTPUT+'_MERGED_all-positions.txt '+ OUTPUT+'_MERGED_only-informative-sites.txt',shell=True)
        subprocess.call('mv '+OUTPUT+'_MERGED_all-positions.txt.bak '+ OUTPUT+'_MERGED_all-positions.txt',shell=True)

        fileinput.close()

    '''  Make a fasta file with only informative positions'''
    # Ported to python3 

    # only could make it work using pandas (Thyago)
    df = read_table(OUTPUT+'_MERGED_only-informative-sites.txt', header=None)
    informative_positions = list(tuple([tuple(df[col]) for col in df]))

    # print the length of informative positions to see if all have the same length
    print(bcolors.OKBLUE + '\nindex ;','name ;','number of informative positions\n' + bcolors.ENDC)

    for i in range(len(informative_positions)):
        print(i,';',informative_positions[i][0],';',len(informative_positions[i]))

    #create a fasta file
    ofile2 = open(OUTPUT + '_MERGED_only-informative-sites_plus_LPM.fasta',mode='w+')

    for j in range(1,len(informative_positions)): # for each column of the informative position tabular file (except the positions)
            ofile2.write('\n>'+informative_positions[j][0]+'\n') # write the first line of the fasta i.e >TN
            for nucleotide in range(1,len(informative_positions[j])):
                ofile2.write(informative_positions[j][nucleotide]) # write the base

    ofile2.close()
    ifile.close()
    

    ''' Fourth Step: select the LPM bases corresponding to each informative positions of M. leprae dataset --> append the fasta with LPM sequence'''
    print(bcolors.OKBLUE + '\nAdding LPM to fasta...\n' + bcolors.ENDC)

    informative_sites_only = list(informative_positions[0][1:]) # this is a list of the informative positions (ex: 73, 8637, 837272, 33992928)

    # lpm_file = open(LPM,'r') # open the LPM VCF file
    lpm_file = gzip.open(LPM,'rt') # Thyago

    # extract information of vcf in lists
    table = []
    for i in lpm_file:
        table += [i]
    for j in range(len(table)):
        table[j]=table[j].split('\t')
    newtable2 = []
    for k in table:
        if len(k)>2: # exclude header
            newtable2 += [k]


    fasta_list = []
    LPM_pos_dic = {}

    for i in newtable2:
        LPM_pos_dic[i[1]] = ''

    for l in informative_sites_only:
            if str(l) in LPM_pos_dic:
                for i in newtable2:
                    if i[1] == str(l):
                        if i[4] == '.': # REF
                                fasta_list += [i[3]]
                        elif i[4] != '.': # ALT
                                fasta_list += [i[4]]
            elif str(l) not in LPM_pos_dic:
                fasta_list += ['-']

    ifile2 = open(OUTPUT+'_MERGED_only-informative-sites_plus_LPM.fasta','a') # append (a) the fasta file with the LPM sequence corresponding to the M. leprae information sites
    ifile2.write('\n>LPM\n')
    for base in fasta_list:
        ifile2.write(base)

    ifile2.close()
    lpm_file.close()

    # cut every 70 characters using unix code.. i diddn't manage in python :-(
    # Andrej: I commented this out because I prefer to have 1 line/sample

    #subprocess.call("grep -oE '.{1,70}' "+ OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta > "+OUTPUT+'_MERGED_2.fasta',shell=True)
    #subprocess.call("rm "+OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta", shell=True)
    #subprocess.call("mv "+OUTPUT+"_MERGED_2.fasta "+OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta",shell=True )

    ''' Get the number of gaps per genome '''
    print(bcolors.OKBLUE + '\nCounting gaps...\n' + bcolors.ENDC)
    # The code below breaks when fasta is not wrapped. I will include a bash command that will also calculate the % of gaps.

    #if args.gaps:
        #fasta=open(OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta",'r')
        #file=fasta.readlines()
        #newtable3=[]

        #for lines in file:
                #if lines[0]!='>':
                    #newtable3[-1]+=lines
                #else:
                    #newtable3+=[lines]


        #for i in range(len(newtable3)):
                #newtable3[i]=newtable3[i].split('\n')
                #newtable3[i][1:]=[''.join(newtable3[i][1:])]


        #outfile2 = open(OUTPUT+'_gaps-per-genome.txt','w')
        #n = 0
        #for i in range(len(newtable3)):
            #for j in newtable3[i][1]:
                #if j == '-':
                    #n += 1

            #outfile2.write(newtable3[i][0][1:] + '\t'+ str(n)+'\n')
            #n = 0


        #for base in range(len(newtable3)):
            #print newtable3[base][0][1:], len(newtable3[base][1])

        #fasta.close()
        #outfile2.close()

    # Need to use bash becasue the default sh will make Syntax error: "(" unexpected
    subprocess.call("cat <(echo -e 'Sample\\tTotal_gaps\\tGaps_percent') <(paste <(grep '>' "+OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta | tr -d '>') <(grep -A1 '>' "+OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta | sed '/>/d' | while read p; do paste <(echo $p | grep -o '-' | wc -l) <(echo \"scale=4; $(echo $p | grep -o '-' | wc -l)/$(echo $p | wc -m)\" | bc -l | sed 's/^\./0./');done)) > "+OUTPUT+"_MERGED_gaps-per-genome.txt", shell=True, executable='/bin/bash')

    # Finally, since I removed fasta wrapping, the file starts with an empty line, which makes MEGA complain. Remove the 1st line from the fasta:
    # sed -i would be good, but sed behaves differently in Mac, so we will do it the safe way:
    subprocess.call("tail -n+2 "+OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta > temp.fasta &&mv temp.fasta "+OUTPUT+"_MERGED_only-informative-sites_plus_LPM.fasta", shell=True)

    fileinput.close()

    # Finally, gzip the *all-positions.txt
    print(bcolors.OKBLUE + '\nCompressing '+OUTPUT+'_MERGED_all-positions.txt...\n' + bcolors.ENDC)
    subprocess.call('gzip '+OUTPUT+'_MERGED_all-positions.txt',shell=True)

    # And remove OUTPUT+'_all-positions.txt' because it's redundant with '+OUTPUT+'_MERGED_all-positions.txt'
    subprocess.call('rm '+OUTPUT+'_all-positions.txt',shell=True)

print(bcolors.OKGREEN + "ALL DONE!\n" + bcolors.ENDC)
exit(0)