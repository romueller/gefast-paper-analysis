#!/bin/bash

# Performs a performance analysis on the given (dereplicated) data set
# in terms of memory and runtime (non-fastidious and fastidious clustering).



SWARM_V1=tools/Swarm-1.2.3
SWARM_V2=tools/Swarm-2.1.13
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

	# Swarm v1
#	bash scripts/swarm-v1.sh ${SWARM_V1} ${T} ${INFILE} ${OUTDIR}/swarm-v1 ${LOG_FILE}

	# Swarm v2
	bash scripts/swarm-v2.sh ${SWARM_V2} ${T} ${INFILE} ${OUTDIR}/swarm-v2 ${LOG_FILE} 1 1

	# GeFaST
	bash scripts/gefast.sh ${GEFAST} ${T} ${INFILE} ${OUTDIR}/gefast ${LOG_FILE} 0 1 1 1 0 0 1

	# others
#	bash scripts/others.sh ${T} ${INFILE} ${OUTDIR} ${LOG_FILE}

done

# evaluate performance
Rscript --vanilla scripts/eval_performance.R ${LOG_FILE} ${OUTDIR}/${DATA_SET}-performance ${MIN_T} ${MAX_T}