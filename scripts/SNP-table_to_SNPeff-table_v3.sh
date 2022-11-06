#!/bin/bash

# Author Chloe
# Modified by Thyago in v3

# v3: 
#   1) removed arpgarsing to be included within the pipeline directly
#   2) tested on snpEff 5.1 

# v2:
#  1) Added the location of snpEff as new obligatory argument, so that other users can use this script.
#  2) snpEff code changed to fit v4.3 (GCA_000195855.1.29 -> Mycobacterium_leprae_tn)
#  3) MaxOS comes with old coreutils, and some programs work differently than in GNU. Sed has different arguments, cannot understand tabs etc. The easiest is to install GNU-sed (and other GNU coreutils for convenience) in Mac.
#     Check https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/
#     Basically you just execute brew install coreutils, and then brew install gnu-sed (because sed does not come by default wiht coreutils). Like this, all GNU tools start with 'g', so GNU-sed is now gsed.
#  3) Check wheter on Linux or MacOS, and use sed or gsed accordingly.
#  4) Don't refer to my custom snpEff wrapper, rather execute the program through java

INFILE=$1

if [ "$(uname)" == "Linux" ] #if on Linux use sed
	then
		#stitch snpEff restuls (almost the entire command) + rest of the table
		paste <(cat <(echo -e 'POS\tREF\tALT\tEFF\tImpact\tGene\tCodon_change\tA.acid_change') `#include a header for the snpEff part of the table, extract the 1st two columns (coordinate and TN) and infer the alternate alleles:` \
		<(tail -n+2 $INFILE | while read p #read every line, store it in variable $p, and do this:
			do #echo converts tabs to spaces, therefore paste each space-separated column
				paste <(echo $p | cut -d ' ' -f1 ) \
				<(echo $p | cut -d ' ' -f2) \
				<(echo $p | cut -d ' ' -f2- | #grab from 2nd column onward
				tr $(echo $p | cut -d ' ' -f2) '-'| `#translate the letter from the 2nd column (TN) to -, so that all such letters don't go into ALT` \
				tr ' ' '\n' | sort -u | `#get a list of unique letters (ALTs)` \
				tr '\n' '\t') | `#after sorting, the first character after the coordinate will be '-', we don't need it` \
				cut -f1,2,4- | `#after sorting, the first character after the coordinate will be '-', we don't need it` \
				tr '\t' ',' | sed -r 's/^([0-9]+),([ATGC]),(.+),$/Chromosome\t\1\t.\t\2\t\3/g' #convert to format for snpEff (Chromosome	POS	.	REF	ALT[s])
			done |
		snpEff eff -ud 0 Mycobacterium_leprae_tn | #run snpEff
		grep -v '#' | #omit header
		cut -f 2,4,5,8 | #grab only columns that we want
		sed -r -e 's/WARNING_TRANSCRIPT_NO_START_CODON//g' -e 's/ANN=//g' -e 's/\|,/\|\t/g' | #remove 'ANN=' so that ALT is first character of that column, and in case there were more ALTs, separate them with a TAB
		while read x
			do
				if [ $(echo $x | wc -w) = 4 ] #if there was only one ALT
				then paste <(echo $x | cut -d ' ' -f1-3) <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 2,3,4,10,11)
				elif [ $(echo $x | wc -w) = 5 ] #if there were 2 ALTs
				then paste <(echo $x | cut -d ' ' -f1-3) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 2) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 2)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 3) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 3)) <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 4) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 10) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 10)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 11) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 11))
				elif [ $(echo $x | wc -w) = 6 ] #if there were 3 ALTs
				then paste <(echo $x | cut -d ' ' -f1-3) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 2) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 2) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 2)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 3) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 3) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 3)) <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 4) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 10) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 10) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 10)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 11) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 11) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 11))
				else
					echo >&2 "##############################################
There were more than 3 ALTs at position $(echo $x | cut -d ' ' -f1)!
This is not normal for SNPs.
Check what is going on or modifiy this script to handle it!
##############################################
"
					exit 0
				fi
		done | # Clean up the output:
		sed -r -e 's/[\| ]/\t/g' -e 's/\t[cn]\./\t/g' -e 's/\tp\./\t/g' -e 's/;[cn]\./;/g' -e 's/;p\./;/g')) `#finally append the rest of the table` \
		<(cut -f2- $INFILE)

# ----------------------------- MAC OS VERSION -----------------------------

elif [ "$(uname)" == "Darwin" ] #if on MacOS use gsed
	then
		#stitch snpEff restuls (almost the entire command) + rest of the table
		paste <(cat <(echo -e 'POS\tREF\tALT\tEFF\tImpact\tGene\tCodon_change\tA.acid_change') `#include a header for the snpEff part of the table, extract the 1st two columns (coordinate and TN) and infer the alternate alleles:` \
		<(tail -n+2 $INFILE | while read p #read every line, store it in variable $p, and do this:
			do #echo converts tabs to spaces, therefore paste each space-separated column
				paste <(echo $p | cut -d ' ' -f1 ) \
				<(echo $p | cut -d ' ' -f2) \
				<(echo $p | cut -d ' ' -f2- | #grab from 2nd column onward
				tr $(echo $p | cut -d ' ' -f2) '-'| `#translate the letter from the 2nd column (TN) to -, so that all such letters don't go into ALT` \
				tr ' ' '\n' | sort -u | `#get a list of unique letters (ALTs)` \
				tr '\n' '\t') | `#after sorting, the first character after the coordinate will be '-', we don't need it` \
				cut -f1,2,4- | `#after sorting, the first character after the coordinate will be '-', we don't need it` \
				tr '\t' ',' | gsed -r 's/^([0-9]+),([ATGC]),(.+),$/Chromosome\t\1\t.\t\2\t\3/g' #convert to format for snpEff (Chromosome	POS	.	REF	ALT[s])
			done |
		snpEff eff -ud 0 Mycobacterium_leprae_tn | #run snpEff
		grep -v '#' | #omit header
		cut -f 2,4,5,8 | #grab only columns that we want
		gsed -r -e 's/WARNING_TRANSCRIPT_NO_START_CODON//g' -e 's/ANN=//g' -e 's/\|,/\|\t/g' | #remove 'ANN=' so that ALT is first character of that column, and in case there were more ALTs, separate them with a TAB
		while read x
			do
				if [ $(echo $x | wc -w) = 4 ] #if there was only one ALT
				then paste <(echo $x | cut -d ' ' -f1-3) <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 2,3,4,10,11)
				elif [ $(echo $x | wc -w) = 5 ] #if there were 2 ALTs
				then paste <(echo $x | cut -d ' ' -f1-3) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 2) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 2)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 3) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 3)) <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 4) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 10) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 10)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 11) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 11))
				elif [ $(echo $x | wc -w) = 6 ] #if there were 3 ALTs
				then paste <(echo $x | cut -d ' ' -f1-3) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 2) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 2) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 2)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 3) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 3) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 3)) <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 4) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 10) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 10) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 10)) <(paste -d ';' <(echo $x | cut -d ' ' -f4 | cut -d '|' -f 11) <(echo $x | cut -d ' ' -f5 | cut -d '|' -f 11) <(echo $x | cut -d ' ' -f6 | cut -d '|' -f 11))
				else
					echo >&2 "##############################################
There were more than 3 ALTs at position $(echo $x | cut -d ' ' -f1)!
This is not normal for SNPs.
Check what is going on or modifiy this script to handle it!
##############################################
"
					exit 0
				fi
		done | # Clean up the output:
		gsed -r -e 's/[\| ]/\t/g' -e 's/\t[cn]\./\t/g' -e 's/\tp\./\t/g' -e 's/;[cn]\./;/g' -e 's/;p\./;/g')) `#finally append the rest of the table` \
		<(cut -f2- $INFILE)

	else
		echo>&2 "##############################################
Cannot determine whether on Linux or MacOS :(
Check what is going on or modifiy this script to handle it!
##############################################
"
		exit 0

	fi

exit 0