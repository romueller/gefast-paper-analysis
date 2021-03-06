#!/bin/bash

# Performs a performance analysis on the given (dereplicated) data set
# and subsets of it in order to investigate the development of memory 
# consumption and runtime with increasing data set size.



SWARM_V1=tools/Swarm-1.2.3
SWARM_V2=tools/Swarm-2.1.13
GEFAST=tools/GeFaST


DATA_SET=$1
MIN_T=$2
MAX_T=$3
REPETITIONS=$4

FASTA_FILE=data/${DATA_SET}/${DATA_SET}_derep.fasta
PERCENTAGES=5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100
SUBSAMPLE_STEM=data/${DATA_SET}/${DATA_SET}_derep_sub
REMOVE_SUBSAMPLES=1
OUTDIR=results
LOG_FILE=${OUTDIR}/${DATA_SET}-subsampling-measurements.csv


# subsample data
./tools/FastaSampler -r ${FASTA_FILE} ${PERCENTAGES} ${REPETITIONS} ${SUBSAMPLE_STEM}


# run tools on subsamples
echo "time,memory,cmd" > ${LOG_FILE}

for ((T=${MIN_T}; T<=${MAX_T}; T++))
do
	for ((P=5; P<=100; P=P+5))
	do
		for ((R=0; R<${REPETITIONS}; R++))
		do

			INFILE=${SUBSAMPLE_STEM}_${P}_${R}.fasta

			# Swarm v1
#			bash scripts/swarm.sh ${SWARM} ${T} ${INFILE} ${OUTDIR}/swarm-v1 ${LOG_FILE}

			# Swarm v2
			bash scripts/swarm.sh ${SWARM} ${T} ${INFILE} ${OUTDIR}/swarm-v2 ${LOG_FILE} 1 1

			# GeFaST
			bash scripts/gefast.sh ${GEFAST} ${T} ${INFILE} ${OUTDIR}/gefast ${LOG_FILE} 0 1 1 1 0 0 1

			# others
#			bash scripts/others.sh ${T} ${INFILE} ${OUTDIR} ${LOG_FILE}

		done
	done
done

# evaluate performance
Rscript --vanilla scripts/eval_subsampling.R ${LOG_FILE} ${OUTDIR}/${DATA_SET}-subsampling

# clean up
if [ "${REMOVE_SUBSAMPLES}" == "1" ]; then
	rm ${SUBSAMPLE_STEM}*
fi
