#!/bin/bash

# Organises the measured (time & memory) execution of Swarm with 
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
	${LOG_CMD} ${SWARM} -d ${THRESHOLD} -a 1 -o ${RESULTS}_o_${THRESHOLD}.csv ${INFILE}
	rm ${RESULTS}_o_${THRESHOLD}.csv
fi

# fastidious
if [ "${THRESHOLD}" == "1" ] && [ "${FASTIDIOUS}" == "1" ]; then
	${LOG_CMD} ${SWARM} -f -d ${THRESHOLD} -a 1 -o ${RESULTS}_o_${THRESHOLD}_f.csv ${INFILE}
	rm ${RESULTS}_o_${THRESHOLD}_f.csv
fi

