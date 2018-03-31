#!/bin/bash

# Organises the measured (runtime & memory) execution of a collection 
# clustering tools using a percentage threshold.



### tools
USEARCH=tools/usearch10.0.240_i86linux32 # path to USEARCH binary, adjust to your system
VSEARCH=tools/vsearch
CDHIT=tools/cd-hit
DNACLUST=tools/dnaclust
SUMACLUST=tools/sumaclust

### configuration
THRESHOLD=0.$((100 - $1))
INFILE=$2
OUTDIR=$3

LOG_FILE=$4
LOG_CMD="/usr/bin/time -f %e,%M,%C -a -o ${LOG_FILE}"

INFILE_ALT=${INFILE}_alt
INFILE_ALT_LENGTH=${INFILE}_alt_length
INFILE_ALT_ABUNDANCE=${INFILE}_alt_abundance



### create input files for USEARCH and VSEARCH (with VSEARCH to avoid memory limit of 32-bit USEARCH)
sed 's/_/;size=/' ${INFILE} > ${INFILE_ALT} 
${VSEARCH} -sortbylength ${INFILE_ALT} -output ${INFILE_ALT_LENGTH} -minsize 1
${VSEARCH} -sortbysize ${INFILE_ALT} -output ${INFILE_ALT_ABUNDANCE} -minsize 1


### clustering

## USEARCH

# USEARCH (cluster_fast, sort by length)
RES=${OUTDIR}/usearch-fast-length_${THRESHOLD}.csv
${LOG_CMD} ${USEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id ${THRESHOLD} -uc ${RES} -sort length -fulldp
rm ${RES}

# USEARCH (cluster_fast, sort by abundance)
RES=${OUTDIR}/usearch-fast-abund_${THRESHOLD}.csv
${LOG_CMD} ${USEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id ${THRESHOLD} -uc ${RES} -sort size -fulldp
rm ${RES}

# USEARCH (cluster_smallmem, presorted by length)
RES=${OUTDIR}/usearch-small-length_${THRESHOLD}.csv
${LOG_CMD} ${USEARCH} -cluster_smallmem ${INFILE_ALT_LENGTH} -id ${THRESHOLD} -uc ${RES} -sortedby length -fulldp
rm ${RES}

# USEARCH (cluster_smallmem, presorted by abundance)
RES=${OUTDIR}/usearch-small-abund_${THRESHOLD}.csv
${LOG_CMD} ${USEARCH} -cluster_smallmem ${INFILE_ALT_ABUNDANCE} -id ${THRESHOLD} -uc ${RES} -sortedby size -fulldp
rm ${RES}


## VSEARCH

# VSEARCH (cluster_fast, sort by length)
RES=${OUTDIR}/vsearch-fast-length_${THRESHOLD}.csv
${LOG_CMD} ${VSEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id ${THRESHOLD} -uc ${RES} -minsize 1
rm ${RES}

# VSEARCH (cluster_size, sort by abundance)
RES=${OUTDIR}/vsearch-size-abundance_${THRESHOLD}.csv
${LOG_CMD} ${VSEARCH} -threads 1 -cluster_size ${INFILE_ALT} -id ${THRESHOLD} -uc ${RES} -minsize 1
rm ${RES}

# VSEARCH (cluster_smallmem, presorted by length)
RES=${OUTDIR}/vsearch-small-length_${THRESHOLD}.csv
${LOG_CMD} ${VSEARCH} -threads 1 -cluster_smallmem ${INFILE_ALT_LENGTH} -id ${THRESHOLD} -uc ${RES} -minsize 1
rm ${RES}

# VSEARCH (cluster_smallmem, presorted by abundance)
RES=${OUTDIR}/vsearch-small-abund_${THRESHOLD}.csv
${LOG_CMD} ${VSEARCH} -threads 1 -cluster_smallmem ${INFILE_ALT_ABUNDANCE} -id ${THRESHOLD} -uc ${RES} -minsize 1 -usersort
rm ${RES}


## misc

# CD-HIT (-d 0 = use full sequence name; -T 1 = use one thread; -M 0 = no memory limit)
RES=${OUTDIR}/cdhit_${THRESHOLD}.csv
${LOG_CMD} ${CDHIT} -i ${INFILE} -o ${RES} -c ${THRESHOLD} -d 0 -T 1 -M 0
rm ${RES}*

# DNACLUST (-t 1 = use one thread)
RES=${OUTDIR}/dnaclust_${THRESHOLD}.csv
${LOG_CMD} ${DNACLUST} -s ${THRESHOLD} -i ${INFILE} -t 1 > ${RES}
rm ${RES}

# Sumaclust
RES=${OUTDIR}/sumaclust_${THRESHOLD}.csv
${LOG_CMD} ${SUMACLUST} -t ${THRESHOLD} -O ${RES} ${INFILE} > /dev/null
rm ${RES}


### clean up 
rm ${INFILE_ALT} ${INFILE_ALT_LENGTH} ${INFILE_ALT_ABUNDANCE}
