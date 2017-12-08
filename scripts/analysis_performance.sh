#!/bin/bash

# Performs a performance analysis on the given (dereplicated) data set
# in terms of memory and runtime (non-fastidious and fastidious clustering).



SWARM=tools/Swarm-2.1.13
GEFAST=tools/GeFaST


DATA_SET=$1
MIN_T=$2
MAX_T=$3

INFILE=data/${DATA_SET}/${DATA_SET}_derep.fasta
OUTDIR=results
LOG_FILE=${OUTDIR}/${DATA_SET}-performance-measurements.csv


# run tools
echo "time,memory,cmd" > ${LOG_FILE}

for ((T=${MIN_T}; T<=${MAX_T}; T++))
do

	# Swarm
	bash scripts/swarm.sh ${SWARM} ${T} ${INFILE} ${OUTDIR}/swarm ${LOG_FILE} 1 1

	# GeFaST
	bash scripts/gefast.sh ${GEFAST} ${T} ${INFILE} ${OUTDIR}/gefast ${LOG_FILE} 0 1 1 1 0

done

# evaluate performance
Rscript --vanilla scripts/eval_performance.R ${LOG_FILE} ${OUTDIR}/${DATA_SET}-performance ${MIN_T} ${MAX_T}