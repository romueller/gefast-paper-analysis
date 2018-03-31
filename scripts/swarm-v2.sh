#!/bin/bash

# Organises the measured (time & memory) execution of Swarm (v2) with 
# different parameters / options.



### configuration
SWARM=$1
THRESHOLD=$2
INFILE=$3
RESULTS=$4

LOG_FILE=$5
LOG_CMD="/usr/bin/time -f %e,%M,%C -a -o ${LOG_FILE}"

NON_FASTIDIOUS=$6
FASTIDIOUS=$7


### clustering

# non-fastidious
if [ "${NON_FASTIDIOUS}" == "1" ]; then
	RES=${RESULTS}_o_${THRESHOLD}.csv
	${LOG_CMD} ${SWARM} -d ${THRESHOLD} -o ${RES} -a 1 ${INFILE}
	rm ${RES}
fi

# fastidious
if [ "${THRESHOLD}" == "1" ] && [ "${FASTIDIOUS}" == "1" ]; then
	RES=${RESULTS}_o_${THRESHOLD}_f.csv
	${LOG_CMD} ${SWARM} -f -d ${THRESHOLD} -o ${RES} -a 1 ${INFILE}
	rm ${RES}
fi

