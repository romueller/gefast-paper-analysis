#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
#
# Builds confusion table and computes metrics



# scripts
BUILD_CONFUSION_TABLE=scripts/confusion_table_reduced.py
COMPUTE_CLUSTER_METRICS=scripts/CValidate_reduced.pl

# inputs
CLUSTERING_RESULTS=$1
METHOD=$2
THRESHOLD=$3
TAXONOMIC_ASSIGNMENTS=$4
CONFUSION_TABLE=$5
METRICS_FILE=$6
REPETITION=$7

# build confusion table; derive and store metric values
python ${BUILD_CONFUSION_TABLE} -t ${TAXONOMIC_ASSIGNMENTS} -s ${CLUSTERING_RESULTS} > ${CONFUSION_TABLE}
METRICS_DATA=$(perl ${COMPUTE_CLUSTER_METRICS} --cfile=${CONFUSION_TABLE})
echo "${METHOD},${THRESHOLD},${REPETITION},${METRICS_DATA}" >> ${METRICS_FILE}

