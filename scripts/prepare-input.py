import os
import re
import glob
import argparse
from collections import defaultdict
from shutil import copyfile
from shutil import move
import logging
import pandas as pd
logging.basicConfig(filename="prepare-input_script.log", level=logging.DEBUG, \
format=' %(asctime)s - %(levelname)s - %(message)s')

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Parsing arguments ------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=bcolors.BOLD + bcolors.HEADER + '***_***_***_*** Prepare FASTQs \
     for the pipeline. ***_***_***_***\n\n' + bcolors.ENDC)
parser.add_argument(
    '-i', '--input', help='Absolute path to directory storing FASTQs.\n \
    Ex.: /home/usuario/projeto1/fastqs', required=True)
parser.add_argument('-o', '--output',
                    help='Output directory if different than the default: \
                    ./prepared-fastqs/', default='./prepared-fastqs/')

args = parser.parse_args()

# defining variables -----------------------------------------------------------
FASTQ_END = ('*.fastq.gz', '*.fq.gz', '*.fq', '*.fastq')
INPUT = os.path.abspath(args.input)
OUTPUT = os.path.abspath(os.path.join(INPUT, args.output))

print(bcolors.HEADER + '\n=================================\n'
                              'Running...\n' +
                         '=================================\n' +
    bcolors.ENDC + '\n')

# check dir
assert os.path.isdir(INPUT), bcolors.FAIL + \
    'ERROR: {} directory does not exist. Please check.'.format(INPUT) \
    + bcolors.ENDC

# change to input dir
os.chdir(os.path.abspath(INPUT))

print(bcolors.OKGREEN + "Changed working directory to: {}\n\n".format(os.getcwd()) +
bcolors.ENDC)

# create output dir
if not os.path.isdir(OUTPUT):
    os.mkdir(os.path.join(os.getcwd(), OUTPUT))

assert os.path.isdir(OUTPUT), bcolors.FAIL + \
    '{} directory does not exist. Please check.'.format(OUTPUT) \
    + bcolors.ENDC

# Regex dictionary containing patterns. This can be expanded as needed

patterns = {
    'sample_code_lane_fixed':
        '(^[A-Za-z0-9\.\-]+)\_(S[\d]+)\_(L[00\d]+)\_(001)(\.[a-z\.]+)$',
    'sample_code_lane_read_fixed':
        '(^[A-Za-z0-9\.\-]+)\_(S[\d]+)\_(L[00\d]+)\_(R[12])\_(001)(\.[a-z\.]+)$',
    'sample_lane_fixed':
        '(^[A-Za-z0-9\-\.]+)\_(L[00\d]+)\_(001)(\.[a-z\.]+)$',
    'sample_lane_read':
        '(^[A-Z-a-z0-9\-\.]+)\_(L00\d+)\_(R[12])(\.[a-z\.]+)$',
    'sample_code_lane_read':
        '(^[A-Z-a-z0-9\-\.]+)\_(S\d{1,2})\_(L00\d+)\_(R[12])(\.[a-z\.]+)$',
    'sample':
        '(^[a-zA-Z0-9\-\.]+)(\.[a-z\.]+)$',
    'sample_else':
        '^([a-zA-Z0-9\-\.]+)\_(([a-zA-Z\-]+)(?=.fastq.*|.fq.*))([.a-z]+)$',    
    'sample_lane_read_fixed':
        '(^[a-zA-Z0-9\-\.]+)\_(L00\d)\_(R[12])\_(001)(\.[a-z\.]+)$',
    'sample_code_fixed':
        '(^[A-Za-z0-9\-\.]+)\_(S\d{1,2})\_(001)(\.[a-z\.]+)$',
    'sample_lane':
        '(^[0-9A-Za-z\-\.]+)\_(L00\d)(\.[a-z\.]+)$',
    'sample_read':
        '(^[A-Za-z0-9\-\.]+)\_(R[12])(\.[a-z\.]+)$',
    'sample_fixed':
        '(^[A-Za-z0-9\-\.]+)\_(001)(\.[a-z\.]+)$',
    'sample_code_read_fixed':
        '(^[A-Za-z0-9\-\.]+)\_(S\d{1,2})\_(R[12])\_(001)(\.[a-z\.]+)$',
    'sample_fixed_read':
        '(^[A-Za-z0-9\-\.]+)\_(001)_(R[12])(\.[a-z\.]+)$',
    'sample_read_fixed':
            '(^[A-Za-z0-9\-\.]+)\_(R[12])\_(001)(\.[a-z\.]+)$',
    'sample_code_read':
        '(^[A-Za-z0-9\-\.]+)\_(S\d{1,2})\_(R[12])(\.[a-z\.]+)$',
}

# Use glob to get the files in the specified, i.e., current directory

def get_fastqs_paths(path=os.getcwd()):
    """ List FASTQ files in current directory."""
    fls = []
    for ext in FASTQ_END:
        fls.extend(glob.glob(os.path.join(path, ext)))
    return sorted(fls)

# Get files from workdir -------------------------------------------------------
files = [os.path.basename(i) for i in get_fastqs_paths()]

# check if any files was captured
assert len(files) > 0, bcolors.FAIL + 'No expected files found in the specified directory.\n \
        Make sure the files end with .fastq.gz, .fq.gz, .fq, or .fastq' + bcolors.ENDC

print(bcolors.OKGREEN + 'Found {} files in the specified directory.'.format(len(files)) + "\n" + bcolors.ENDC)

# Fix Read naming --------------------------------------------------------------

def detect_wrong_pair(file_list):
    """
    Detect if paired files uses _1/_2 instead of _R1/_R2 naming convention \
    and renames it in the same dir. Otherwise, the regex wont capture.
    """

    global files
    wrong_pair = []
    for file in file_list:
        if re.search('\_([12]+)\.([a-z\.]+)$', file) is not None:
            wrong_pair.append(file)

    if len(wrong_pair) > 0:
        print(bcolors.BOLD + bcolors.OKCYAN + 'Found {} files with wrong pair naming convention.'.format(
            len(wrong_pair)) + "\n" + bcolors.ENDC)

        for file in wrong_pair:
            src = file
            dest = re.sub('\_([12]+)', '_R\\1', file)
            try:
                copyfile(src, dest)
                print(bcolors.OKCYAN + 'Renamed {} to: {}'.format(file, dest) + bcolors.ENDC)
            except FileExistsError:
                print(bcolors.WARNING + "File already exists. Skipping..." + bcolors.ENDC)
                continue

    print(bcolors.OKGREEN + '\n=========================================\n'
                            'Step1: Done fixing read pair naming scheme \n' +
                            '=========================================\n' +
bcolors.ENDC + '\n')

# call fix read naming
detect_wrong_pair(files)

# get files in directory again to capture newly renamed
files = [os.path.basename(i) for i in get_fastqs_paths()]

# Identify files from patterns -------------------------------------------------

fastqs = defaultdict(list)
tmp = defaultdict(list)
for file in files:
    for name, pattern in patterns.items():
        match = re.search(pattern, file)
        if type(match) == re.Match:
            fastqs[name].append(list(match.groups()))
            tmp[name].append(match.group())


for key,value in fastqs.items():
    print(key, value)

# check if a file was captured by more than one pattern and stop
vals = [file for l in tmp.values() for file in l]
assert len(list(set(vals))) == len(vals), bcolors.FAIL + "Some files were matched by more than one pattern! This is dangerous. Please check!\n" + bcolors.ENDC

# GZIP uncompressed files ------------------------------------------------------

def gzip_uncompressed(input_file):
    """ Test if gile is gzipped and if not compress it in same dir."""

    if isinstance(input_file, list):
        length = len(input_file)
    else:
        length = 1
        input_file = [input_file]

    index = 0
    while index < length:
        if input_file[index].endswith('.fastq') or input_file[index].endswith('.fq'):
            assert os.path.isfile(input_file[index]), bcolors.FAIL + 'File not found: {}'.format(input_file[index]) + bcolors.ENDC
            print(bcolors.OKCYAN + 'Compressing {}'.format(input_file[index]) + bcolors.ENDC)
            os.system("gzip -v {}".format(input_file[index]))
        index += 1
        continue

    print(bcolors.OKGREEN + '\n=========================================\n'
                                  'Step2: Done compressing files\n' +
                            '=========================================\n' +
bcolors.ENDC + '\n')

gzip_uncompressed(files)

# Capture files again to account for newly gzipped
files = [os.path.basename(i) for i in get_fastqs_paths()]

# get files again to newly gzipped files
fastqs = defaultdict(list)
for file in files:
    for name, pattern in patterns.items():
        match = re.search(pattern, file)
        if type(match) == re.Match:
            fastqs[name].append(list(match.groups()))

# Resolve files-----------------------------------------------------------------

def resolve_files(matches):
    """
    This function detects fastq file in a directory and merge their lanes.
    The script also tries to fix names and remove unwanted parts of official
    Illumina FASTQ names, like _001 fixed region and sample code (SXX).

    Once the files are listed, the user can choose to proceed or not.
    """

    print("Listing files...\n")

    inst = list()
    global files

    if os.path.isfile('cmds_used.txt'):
        os.unlink('cmds_used.txt')

    # loop over each key in matches
    for key in matches.keys():
        # single-end files

        if key == 'sample_code_lane_fixed':
            for file in matches[key]:
                old = file[0] + '_' + file[1] + \
                    '_' + "L00*_" + file[3] + file[-1]
                new = os.path.join(OUTPUT, file[0] + file[-1])
                cmd_str = 'os.system("cat {} > {}")'.format(old, new)
                inst.append(cmd_str)

        elif key == 'sample_code_fixed':
            for file in matches[key]:
                old = file[0] + '_' + file[1] + '_' + file[2] + file[-1]
                new = os.path.join(OUTPUT, file[0] + file[-1])
                inst.append('copyfile("{}", "{}")'.format(old, new))

        elif key == 'sample_fixed':
            for file in matches[key]:
                old = file[0] + '_' + file[1] + file[-1]
                new = os.path.join(OUTPUT, file[0] + file[-1])
                inst.append('copyfile("{}", "{}")'.format(old, new))

        elif key == 'sample_lane_fixed':
            for file in matches[key]:
                old = file[0] + "_" + "L00*" + "_" + file[2] + file[-1]
                new = os.path.join(OUTPUT, file[0] + file[-1])
                cmd_str = 'os.system("cat {} > {}")'.format(old, new)
                inst.append(cmd_str)

        elif key == "sample_lane":
            for file in matches[key]:
                old = file[0] + "_L00*" + file[-1]
                new = os.path.join(OUTPUT, file[0] + file[-1])
                cmd_str = 'os.system("cat {} > {}")'.format(old, new)
                inst.append(cmd_str)

        elif key == "sample":
            for file in matches[key]:
                orig = file[0] + file[-1]
                inst.append('copyfile("{}", "{}")'.format(orig, os.path.join(OUTPUT, orig)))

        elif key == "sample_else":
            for file in matches[key]:
                old = file[0] + '_' + file[1] + file[-1]
                new = file[0] + file[-1]
                inst.append('copyfile("{}", "{}")'.format(old, os.path.join(OUTPUT, new)))

        # paired-end files
        elif key == 'sample_code_lane_read_fixed':
            for file in matches[key]:
                old_r1 = file[0] + '_' + file[1] + '_' + file[2] + \
                    '_R1_' + file[4] + file[-1]
                old_r2 = file[0] + '_' + file[1] + '_' + file[2] + \
                    '_R2_' + file[4] + file[-1]
                old_r1_sub = re.sub('L00\d', 'L00*', old_r1)
                old_r2_sub = re.sub('L00\d', 'L00*', old_r2)
                new1 = os.path.join(OUTPUT, file[0] + '_R1' + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + '_R2' + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str)
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str)
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str)
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str)

        elif key == 'sample_code_lane_read':
            for file in matches[key]:
                old_r1 = file[0] + '_' + file[1] + '_' + file[2] + \
                    '_R1' + file[-1]
                old_r2 = file[0] + '_' + file[1] + '_' + file[2] + \
                    '_R2' + file[-1]
                old_r1_sub = re.sub('L00\d', 'L00*', old_r1)
                old_r2_sub = re.sub('L00\d', 'L00*', old_r2)
                new1 = os.path.join(OUTPUT, file[0] + '_R1' + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + '_R2' + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str)
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str)
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str)
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str)

        elif key == 'sample_code_read':
            for file in matches[key]:
                old_r1 = file[0] + "_" + file[1] + "_R1" + file[-1]
                old_r2 = file[0] + "_" + file[1] + "_R2" + file[-1]
                new1 = os.path.join(OUTPUT, file[0] + "_R1" + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + "_R2" + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))

        elif key == 'sample_fixed_read':
            for file in matches[key]:
                old_r1 = file[0] + "_" + file[1] + "_R1" + file[-1]
                old_r2 = file[0] + "_" + file[1] + "_R2" + file[-1]
                new1 = os.path.join(OUTPUT, file[0] + "_R1" + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + "_R2" + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))

        elif key == 'sample_code_read_fixed':
            for file in matches[key]:
                old_r1 = file[0] + '_' + file[1] + '_R1_' + file[3] + file[-1]
                old_r2 = file[0] + '_' + file[1] + '_R2_' + file[3] + file[-1]
                new1 = os.path.join(OUTPUT, file[0] + '_R1' + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + '_R2' + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    cmd_str1 = 'copyfile("{}", "{}")'.format(old_r1, new1)
                    cmd_str2 = 'copyfile("{}", "{}")'.format(old_r2, new2)
                    inst.append(cmd_str1)
                    inst.append(cmd_str2)
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    cmd_str1 = 'copyfile("{}", "{}")'.format(old_r1, new1)
                    inst.append(cmd_str1)
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    cmd_str2 = 'copyfile("{}", "{}")'.format(old_r2, new2)
                    inst.append(cmd_str2)

        elif key == "sample_lane_read":
            for file in matches[key]:
                old_r1 = file[0] + "_" + file[1] + "_R1" + file[-1]
                old_r2 = file[0] + "_" + file[1] + "_R2" + file[-1]
                old_r1_sub = re.sub('L00\d', 'L00*', old_r1)
                old_r2_sub = re.sub('L00\d', 'L00*', old_r2)
                new1 = os.path.join(OUTPUT, file[0] + "_R1" + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + "_R2" + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    cmd_str = 'os.system("cat {}  > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str)
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str)
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str)
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read2
                    cmd_str = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str)

        elif key == "sample_lane_read_fixed":
            for file in matches[key]:
                old_r1 = file[0] + "_" + file[1] + "_" + file[2] + "_" + file[3] + file[-1]
                old_r2 = file[0] + "_" + file[1] + "_" + file[2] + "_" + file[3] + file[-1]
                old_r1_sub = re.sub('L00\d', 'L00*', old_r1)
                old_r2_sub = re.sub('L00\d', 'L00*', old_r2)
                new1 = os.path.join(OUTPUT, file[0] + "_R1" + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + "_R2" + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    cmd_str1 = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    cmd_str2 = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str1)
                    inst.append(cmd_str2)
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    cmd_str1 = 'os.system("cat {} > {}")'.format(old_r1_sub, new1)
                    inst.append(cmd_str1)
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    cmd_str2 = 'os.system("cat {} > {}")'.format(old_r2_sub, new2)
                    inst.append(cmd_str2)

        elif key == "sample_read_fixed":
            for file in matches[key]:
                old_r1 = file[0] + "_R1" + '_' + file[2] + file[-1]
                old_r2 = file[0] + "_R2" + '_' + file[2] + file[-1]
                new1 = os.path.join(OUTPUT, file[0] + "_R1" + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + "_R2" + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))

        elif key == "sample_read":
            for file in matches[key]:
                old_r1 = file[0] + "_R1" + file[-1]
                old_r2 = file[0] + "_R2" + file[-1]
                new1 = os.path.join(OUTPUT, file[0] + "_R1" + file[-1])
                new2 = os.path.join(OUTPUT, file[0] + "_R2" + file[-1])
                if os.path.isfile(old_r1) and os.path.isfile(old_r2):
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))
                elif os.path.isfile(old_r1):
                    new1 = re.sub('\_R[1,2]', '', new1) # remove orphan 'paired' read2
                    inst.append('copyfile("{}", "{}")'.format(old_r1, new1))
                elif os.path.isfile(old_r2):
                    new2 = re.sub('\_R[1,2]', '', new2) # remove orphan 'paired' read1
                    inst.append('copyfile("{}", "{}")'.format(old_r2, new2))

    print(bcolors.OKBLUE + "The following code will be executed:\n" + bcolors.ENDC)
    cmds = open("cmds_used.txt", 'a')

    for j,cmd in enumerate(sorted(list(set(inst)))):
        print(j+1, "\t", cmd)
        cmds.write(cmd+'\n')
        try:
            exec(cmd)
        except Exception as e:
            logging.error(traceback.format_exc())
    cmds.close()

    print(bcolors.OKGREEN + '\n=========================================\n'
                                 'Step 3: Done resolving filenames\n' +
                              '=========================================\n' +
    bcolors.ENDC + '\n')

resolve_files(fastqs)

print(bcolors.OKCYAN + "Normalizing extension of files..." + bcolors.ENDC)

prep_files = glob.glob(os.path.join(OUTPUT, "*"))

def normalized_fastqs(prep_files):
    """ Make sure all FASTQs ends with fastq.gz and not fq.gz """

    if isinstance(prep_files, list):
        length = len(prep_files)
    else:
        length = 1
        prep_files = [prep_files]

    index = 0
    while index < length:
        if prep_files[index].endswith('fq.gz'):
            assert os.path.isfile(prep_files[index]), bcolors.FAIL + \
            'File not found: {}'.format(prep_files[index]) + bcolors.ENDC
            print(bcolors.OKCYAN + 'Normalizing FASTQ extension {}'.format(
                os.path.basename(prep_files[index])) + "  ->  " +
                os.path.basename(prep_files[index].replace('fq.gz', 'fastq.gz')) +
                bcolors.ENDC)
            move(prep_files[index], prep_files[index].replace('fq.gz', 'fastq.gz'))
            assert os.path.isfile(prep_files[index].replace('fq.gz', 'fastq.gz'))
        index += 1
        continue

    print(bcolors.OKGREEN + '\n=========================================\n'
                              'Step4: Done normalizing extensions\n' +
                            '=========================================\n' +
    bcolors.ENDC + '\n')

normalized_fastqs(prep_files)

print(bcolors.HEADER + "Prepared inputs are in {}".format(OUTPUT) + bcolors.ENDC + '\n')

print(bcolors.OKGREEN + "Creating samples.tsv and units.tsv files..." + bcolors.ENDC + '\n')

output_files = glob.glob(os.path.join(OUTPUT, "*.fastq.gz"))

samples, fq1, fq2, design = [], [], [], []

for file in output_files:
    name = os.path.basename(file)
    sample = re.search('(^[A-Za-z0-9\-]+)', name).group(1)

    if bool(re.search('_R\d', name)) and (file not in fq1 and file not in fq2):
        samples.append(sample)
        fq1.append(glob.glob(os.path.join(OUTPUT, sample) + "*R1*")[0])
        fq2.append(glob.glob(os.path.join(OUTPUT, sample) + "*R2*")[0])
        design.append("paired-end")
    elif file not in fq1 and file not in fq2:
        samples.append(sample)
        fq1.append(file)
        fq2.append('')
        design.append("single-end")

samples_df = pd.DataFrame(columns=['sample'], data=samples)
units_df = pd.DataFrame({'sample': samples, 'fq1': fq1, 'fq2': fq2, 'design': design})

samples_df.to_csv(os.path.join(OUTPUT, 'samples.tsv'), sep='\t', index=False)
units_df.to_csv(os.path.join(OUTPUT, 'units.tsv'), sep='\t', index=False)

print(bcolors.OKCYAN + "samples.tsv and units.tsv files are in {}".format(OUTPUT) + bcolors.ENDC + '\n')

print(bcolors.OKGREEN + '\n=========================================\n'
                              'Step5: Done creating sample.tsv and units.tsv \n' +
                            '=========================================\n' +
                            bcolors.ENDC + '\n')

assert all([(x in units_df['fq1'].to_list() or x in units_df['fq2'].to_list()) for x in output_files]), bcolors.FAIL + \
    "Some files are missing from samples.tsv and units.tsv. Please check.\n" + bcolors.ENDC


exit(0)
