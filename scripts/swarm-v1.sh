#!/bin/bash

# Organises the measured (time & memory) execution of Swarm (v1) with 
# different parameters / options.



### configuration
SWARM=$1
THRESHOLD=$2
INFILE=$3
RESULTS=$4

LOG_FILE=$5
LOG_CMD="/usr/bin/time -f %e,%M,%C -a -o ${LOG_FILE} /bin/sh -c "

SWARM_BREAKER=scripts/swarm_breaker.py


### clustering

# non-fastidious
RES=${RESULTS}_o_${THRESHOLD}.csv
${LOG_CMD} "${SWARM} -d ${THRESHOLD} -o ${RES}.tmp ${INFILE}; \
python ${SWARM_BREAKER} -b ${SWARM} -f ${INFILE} -s ${RES}.tmp -d ${THRESHOLD} > ${RES}"
rm ${RES}*
