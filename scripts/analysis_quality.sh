#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
# 
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
#
# Performs a clustering-quality analysis on the (un)even data set
# similar to the one in the paper cited above



SWARM_V1=tools/Swarm-1.2.3
SWARM_V2=tools/Swarm-2.1.13
GEFAST=tools/GeFaST
CONF=gefast.conf

USEARCH=tools/usearch7.0.1090_i86linux32 # path to USEARCH binary, adjust to your system
SWARM_BREAKER=scripts/swarm_breaker.py
COMPUTE_METRICS=scripts/compute_metrics.sh
REFS=data/rrna_reference.fasta


DATA_SET=$1
MIN_T=$2
MAX_T=$3

INFILE=data/${DATA_SET}/${DATA_SET}_derep.fasta
TAXA_FILE=data/${DATA_SET}/${DATA_SET}.spec.ualn
OUTDIR=results
METRICS_FILE=${OUTDIR}/${DATA_SET}-metrics.csv
CONF_TABLE=${OUTDIR}/confusion-table.conf


# ===== Taxonomic assignment =====

# get reference data
wget -nc -P data/ https://raw.githubusercontent.com/torognes/vsearch-eval/master/cluster/data/rrna_reference.fasta

# assign taxa
bash scripts/assign_taxa.sh ${USEARCH} 8 ${INFILE} ${REFS} ${TAXA_FILE}



# ===== Clustering & confusion table =====

# initialise metrics file
echo "method,threshold,recall,precision,nmi,randindex,adjrandindex" > ${METRICS_FILE}

# run tools, compute & evaluate confusion table
for ((T=${MIN_T}; T<=${MAX_T}; T++))
do

	## Swarm	

	# Swarm (version 1, non-fastidious)
	RES=${OUTDIR}/swarm-v1_o_${T}.csv
	${SWARM_V1} -d ${T} -o ${RES}.tmp ${INFILE}
	python ${SWARM_BREAKER} -b ${SWARM_V1} -f ${INFILE} -s ${RES}.tmp -d ${T} > ${RES}
	rm ${RES}.tmp
	bash ${COMPUTE_METRICS} ${RES} "swarm-v1" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

	# Swarm (version 2, non-fastidious)
	RES=${OUTDIR}/swarm-v2_o_${T}.csv
	${SWARM_V2} -d ${T} -a 1 -o ${RES} ${INFILE} 
	bash ${COMPUTE_METRICS} ${RES} "swarm-v2" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

	# Swarm (version 2, fastidious)
	if [ "${T}" == "1" ]; then
		RES=${OUTDIR}/swarm-v2_s_${T}_f.csv
		${SWARM_V2} -f -d ${T} -a 1 -o ${RES} ${INFILE}
		bash ${COMPUTE_METRICS} ${RES} "swarm-2f-v2" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
		rm ${RES} ${CONF_TABLE}
	fi
	


	## GeFaST

	# GeFaST (edit distance, non-fastidious)
	RES=${OUTDIR}/gefast-e_o_${T}_e.csv
	${GEFAST} ${INFILE} -t ${T} -c ${CONF} -so ${RES}
	bash ${COMPUTE_METRICS} ${RES} "gefast-e" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

	# GeFaST (edit distance, fastidious, t + 1)
	RES=${OUTDIR}/gefast-e_o_${T}_ef1.csv
	${GEFAST} ${INFILE} -t ${T} -c ${CONF} -sf --swarm-fastidious-threshold $((${T} + 1)) -so ${RES}
	bash ${COMPUTE_METRICS} ${RES} "gefast-f1-e" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

	# GeFaST (edit distance, fastidious, 2 * t)
	RES=${OUTDIR}/gefast-e_o_${T}_e2f.csv
	${GEFAST} ${INFILE} -t ${T} -c ${CONF} -sf --swarm-fastidious-threshold $((2 * ${T})) -so ${RES}
	bash ${COMPUTE_METRICS} ${RES} "gefast-2f-e" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}


	# GeFaST (scoring function, non-fastidious)
	RES=${OUTDIR}/gefast-s_o_${T}_s.csv
	${GEFAST} ${INFILE} -t ${T} -c ${CONF} --use-score -so ${RES}
	bash ${COMPUTE_METRICS} ${RES} "gefast-s" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

	# GeFaST (scoring function, fastidious, t + 1)
	RES=${OUTDIR}/gefast-s_o_${T}_sf1.csv
	${GEFAST} ${INFILE} -t ${T} -c ${CONF} --use-score -sf --swarm-fastidious-threshold $((${T} + 1)) -so ${RES}
	bash ${COMPUTE_METRICS} ${RES} "gefast-f1-s" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

	# GeFaST (scoring function, fastidious, 2 * t)
	RES=${OUTDIR}/gefast-s_o_${T}_s2f.csv
	${GEFAST} ${INFILE} -t ${T} -c ${CONF} --use-score -sf --swarm-fastidious-threshold $((2 * ${T})) -so ${RES}
	bash ${COMPUTE_METRICS} ${RES} "gefast-2f-s" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE}
	rm ${RES} ${CONF_TABLE}

done

# evaluate quality
Rscript --vanilla scripts/eval_quality.R ${METRICS_FILE} ${OUTDIR}/${DATA_SET}-metrics ${MIN_T} ${MAX_T}






