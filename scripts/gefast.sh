#!/bin/bash

# Organises the measured (runtime & memory) execution of GeFaST with 
# different parameters / options.



### configuration
GEFAST=$1
THRESHOLD=$2
INFILE=$3
RESULTS=$4

LOG_FILE=$5
LOG_CMD="/usr/bin/time -f %e,%M,%C -a -o ${LOG_FILE}"

CONF=gefast.conf

PREPROCESSING_ONLY=$6
NON_FASTIDIOUS=$7
FASTIDIOUS_INC=$8
FASTIDIOUS_TWICE=$9
FASTIDIOUS_CONST=${10}



### clustering

# preprocessing only
if [ "${PREPROCESSING_ONLY}" == "1" ]; then
	${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --use-score -c ${CONF} --preprocessing-only
fi

# scoring function, non-fastidious
if [ "${NON_FASTIDIOUS}" == "1" ]; then
	${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --use-score -c ${CONF} -so ${RESULTS}_o_${THRESHOLD}.csv
	rm ${RESULTS}_o_${THRESHOLD}.csv
fi

# scoring function, fastidious (t + 1)
if [ "${FASTIDIOUS_INC}" == "1" ]; then
	${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --use-score -c ${CONF} -so ${RESULTS}_o_${THRESHOLD}_f1.csv \
		-sf --swarm-fastidious-threshold $((${THRESHOLD} + 1))
	rm ${RESULTS}_o_${THRESHOLD}_f1.csv
fi

# scoring function, fastidious (2 * t)
if [ "${FASTIDIOUS_TWICE}" == "1" ]; then
	${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --use-score -c ${CONF} -so ${RESULTS}_o_${THRESHOLD}_2f.csv \
		-sf --swarm-fastidious-threshold $((2 * ${THRESHOLD}))
	rm ${RESULTS}_o_${THRESHOLD}_2f.csv
fi


# scoring function, constant fastidious (independent of t)
if [ "${FASTIDIOUS_CONST}" != "0" ]; then
	${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --use-score -c ${CONF} -so ${RESULTS}_o_${THRESHOLD}_cf${FASTIDIOUS_CONST}.csv \
		-sf --swarm-fastidious-threshold ${FASTIDIOUS_CONST}
	rm ${RESULTS}_o_${THRESHOLD}_cf${FASTIDIOUS_CONST}.csv
fi

