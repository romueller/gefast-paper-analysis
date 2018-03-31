#!/bin/bash

# Adapted from:
# https://github.com/torognes/vsearch-eval
#
# File: https://github.com/torognes/vsearch-eval/blob/master/cluster/scripts/tax.sh
#
# Determines the ground truth by matching against the reference sequences



USEARCH=$1
THREADS=$2
FASTA=$3
REFS=$4
ASSIGNMENTS=$5
THRESHOLD=$6

${USEARCH} --usearch_global ${FASTA} \
        --db ${REFS} \
        --id ${THRESHOLD} \
        --blast6out ${ASSIGNMENTS} \
        --strand plus \
        --threads ${THREADS}
