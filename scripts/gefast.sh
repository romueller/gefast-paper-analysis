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
EDIT=${11}
SCORING=${12}



### clustering

# preprocessing only
if [ "${PREPROCESSING_ONLY}" == "1" ]; then
	${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} --use-score --preprocessing-only
fi

# edit-distance mode
if [ "${EDIT}" == "1" ]; then

	# non-fastidious
	if [ "${NON_FASTIDIOUS}" == "1" ]; then
		RES=${RESULTS}-e_o_${THRESHOLD}.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES}
		rm ${RES}
	fi

	# fastidious (t + 1)
	if [ "${FASTIDIOUS_INC}" == "1" ]; then
		RES=${RESULTS}-e_o_${THRESHOLD}_f1.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} \
			-sf --swarm-fastidious-threshold $((${THRESHOLD} + 1))
		rm ${RES}
	fi

	# fastidious (2 * t)
	if [ "${FASTIDIOUS_TWICE}" == "1" ]; then
		RES=${RESULTS}-e_o_${THRESHOLD}_2f.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} \
			-sf --swarm-fastidious-threshold $((2 * ${THRESHOLD}))
		rm ${RES}
	fi

	# constant fastidious (independent of t)
	if [ "${FASTIDIOUS_CONST}" != "0" ]; then
		RES=${RESULTS}-e_o_${THRESHOLD}_cf${FASTIDIOUS_CONST}.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} \
			-sf --swarm-fastidious-threshold ${FASTIDIOUS_CONST}
		rm ${RES}
	fi

fi

# scoring-function mode
if [ "${SCORING}" == "1" ]; then

	# non-fastidious
	if [ "${NON_FASTIDIOUS}" == "1" ]; then
		RES=${RESULTS}-s_o_${THRESHOLD}.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} --use-score
		rm ${RES}
	fi

	# fastidious (t + 1)
	if [ "${FASTIDIOUS_INC}" == "1" ]; then
		RES=${RESULTS}-s_o_${THRESHOLD}_f1.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} --use-score \
			-sf --swarm-fastidious-threshold $((${THRESHOLD} + 1))
		rm ${RES}
	fi

	# fastidious (2 * t)
	if [ "${FASTIDIOUS_TWICE}" == "1" ]; then
		RES=${RESULTS}-s_o_${THRESHOLD}_2f.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} --use-score \
			-sf --swarm-fastidious-threshold $((2 * ${THRESHOLD}))
		rm ${RES}
	fi

	# constant fastidious (independent of t)
	if [ "${FASTIDIOUS_CONST}" != "0" ]; then
		RES=${RESULTS}-s_o_${THRESHOLD}_cf${FASTIDIOUS_CONST}.csv
		${LOG_CMD} ${GEFAST} ${INFILE} -t ${THRESHOLD} --config ${CONF} -so ${RES} --use-score \
			-sf --swarm-fastidious-threshold ${FASTIDIOUS_CONST}
		rm ${RES}
	fi

fi

