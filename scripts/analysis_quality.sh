#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
# 
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
#
# Performs a clustering-quality analysis on the given data set
# similar to the analysis in the paper cited above



SWARM_V1=tools/Swarm-1.2.3
SWARM_V2=tools/Swarm-2.1.13
GEFAST=tools/GeFaST
CONF=gefast.conf
USEARCH=tools/usearch10.0.240_i86linux32 # path to USEARCH binary, adjust to your system
VSEARCH=tools/vsearch
CDHIT=tools/cd-hit
DNACLUST=tools/dnaclust
SUMACLUST=tools/sumaclust


SWARM_BREAKER=scripts/swarm_breaker.py
COMPUTE_METRICS=scripts/compute_metrics.sh


DATA_SET=$1
REFS=$2
MIN_T=$3
MAX_T=$4


GROUND_TRUTH_THRESHOLDS=(0.95 0.97 0.99)

INFILE=data/${DATA_SET}/${DATA_SET}_derep.fasta
INFILE_ALT=${INFILE}_alt
INFILE_ALT_LENGTH=${INFILE}_alt_length
INFILE_ALT_ABUNDANCE=${INFILE}_alt_abundance
OUTDIR=results
CONF_TABLE=${OUTDIR}/confusion-table.conf

LOG_FILE=${OUTDIR}/${DATA_SET}-quality-log.csv
LOG_CMD="/usr/bin/time -f %e,%M,%C -a -o ${LOG_FILE} /bin/sh -c "
echo "time,memory,cmd" > ${LOG_FILE}



# ===== Clustering & confusion table =====

# create input files for USEARCH and VSEARCH (with VSEARCH to avoid memory limit of 32-bit USEARCH)
sed 's/_/;size=/' ${INFILE} > ${INFILE_ALT} 
${VSEARCH} -sortbylength ${INFILE_ALT} -output ${INFILE_ALT_LENGTH} -minsize 1
${VSEARCH} -sortbysize ${INFILE_ALT} -output ${INFILE_ALT_ABUNDANCE} -minsize 1


# prepare ground truths and metrics files
TAXA_FILES=()
METRICS_FILES=()
for P in "${GROUND_TRUTH_THRESHOLDS[@]}"
do

	# assign taxa (with VSEARCH to avoid memory limit of 32-bit USEARCH)
	TAXA_FILE=data/${DATA_SET}/${DATA_SET}_${P}.spec.ualn
	TAXA_FILES+=(${TAXA_FILE})
	bash scripts/assign_taxa.sh ${VSEARCH} 8 ${INFILE} ${REFS} ${TAXA_FILE} ${P}

	# initialise metrics & performance-log file
	METRICS_FILE=${OUTDIR}/${DATA_SET}_${P}-metrics.csv
	METRICS_FILES+=(${METRICS_FILE})
	echo "method,threshold,rep,recall,precision,nmi,randindex,adjrandindex" > ${METRICS_FILE}

done

# run tools, compute & evaluate confusion table
for ((T=${MIN_T}; T<=${MAX_T}; T++))
do

	## Swarm	

	# Swarm (version 1, non-fastidious)
	RES=${OUTDIR}/swarm-v1_o_${T}.csv
	${LOG_CMD} "${SWARM_V1} -d ${T} -o ${RES}.tmp ${INFILE}; \
	python ${SWARM_BREAKER} -b ${SWARM_V1} -f ${INFILE} -s ${RES}.tmp -d ${T} > ${RES}"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "swarm-v1" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# Swarm (version 2, non-fastidious)
	RES=${OUTDIR}/swarm-v2_o_${T}.csv
	${LOG_CMD} "${SWARM_V2} -d ${T} -o ${RES} -a 1 ${INFILE}" 
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "swarm-v2" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}

	# Swarm (version 2, fastidious)
	if [ "${T}" == "1" ]; then
		RES=${OUTDIR}/swarm-v2_s_${T}_f.csv
		${LOG_CMD} "${SWARM_V2} -f -d ${T} -o ${RES} -a 1 ${INFILE}"
		for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
		do
			bash ${COMPUTE_METRICS} ${RES} "swarm-v2-2f" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
		done
		rm ${RES} ${CONF_TABLE}
	fi
	


	## GeFaST

	# GeFaST (edit distance, non-fastidious)
	RES=${OUTDIR}/gefast-e_o_${T}.csv
	${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES}"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "gefast-e" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}

	# GeFaST (edit distance, fastidious, t + 1)
	RES=${OUTDIR}/gefast-e_o_${T}_f1.csv
	${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} -sf --swarm-fastidious-threshold $((${T} + 1))"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "gefast-e-f1" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}

	# GeFaST (edit distance, fastidious, 2 * t)
	RES=${OUTDIR}/gefast-e_o_${T}_2f.csv
	${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} -sf --swarm-fastidious-threshold $((2 * ${T}))"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "gefast-e-2f" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}


	# GeFaST (scoring function, non-fastidious)
	RES=${OUTDIR}/gefast-s_o_${T}.csv
	${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} --use-score"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "gefast-s" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}

	# GeFaST (scoring function, fastidious, t + 1)
	RES=${OUTDIR}/gefast-s_o_${T}_f1.csv
	${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} --use-score -sf --swarm-fastidious-threshold $((${T} + 1))"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "gefast-s-f1" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}

	# GeFaST (scoring function, fastidious, 2 * t)
	RES=${OUTDIR}/gefast-s_o_${T}_2f.csv
	${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} --use-score -sf --swarm-fastidious-threshold $((2 * ${T}))"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "gefast-s-2f" ${T} ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}



	## USEARCH

	# USEARCH (cluster_fast, sort by length)
	RES=${OUTDIR}/usearch-fast-length_${T}.csv
	${LOG_CMD} "${USEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sort length -fulldp"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "usearch-fast-length" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# USEARCH (cluster_fast, sort by abundance)
	RES=${OUTDIR}/usearch-fast-abund_${T}.csv
	${LOG_CMD} "${USEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sort size -fulldp"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "usearch-fast-abund" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# USEARCH (cluster_smallmem, presorted by length)
	RES=${OUTDIR}/usearch-small-length_${T}.csv
	${LOG_CMD} "${USEARCH} -cluster_smallmem ${INFILE_ALT_LENGTH} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sortedby length -fulldp"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "usearch-small-length" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# USEARCH (cluster_smallmem, presorted by abundance)
	RES=${OUTDIR}/usearch-small-abund_${T}.csv
	${LOG_CMD} "${USEARCH} -cluster_smallmem ${INFILE_ALT_ABUNDANCE} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sortedby size -fulldp"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "usearch-small-abund" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}



	## VSEARCH

	# VSEARCH (cluster_fast, sort by length)
	RES=${OUTDIR}/vsearch-fast-length_${T}.csv
	${LOG_CMD} "${VSEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "vsearch-fast-length" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# VSEARCH (cluster_size, sort by abundance)
	RES=${OUTDIR}/vsearch-size-abundance_${T}.csv
	${LOG_CMD} "${VSEARCH} -threads 1 -cluster_size ${INFILE_ALT} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "vsearch-size-abund" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# VSEARCH (cluster_smallmem, presorted by length)
	RES=${OUTDIR}/vsearch-small-length_${T}.csv
	${LOG_CMD} "${VSEARCH} -threads 1 -cluster_smallmem ${INFILE_ALT_LENGTH} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "vsearch-small-length" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# VSEARCH (cluster_smallmem, presorted by abundance)
	RES=${OUTDIR}/vsearch-small-abund_${T}.csv
	${LOG_CMD} "${VSEARCH} -threads 1 -cluster_smallmem ${INFILE_ALT_ABUNDANCE} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1 -usersort"
	awk 'BEGIN {FS="\t"}{
		if ($1 == "S") clusters[$2] = $9
		if ($1 == "H") clusters[$2] = clusters[$2] " " $9
	}
	END {
		for (cluster in clusters) {
			print clusters[cluster]
		}
	}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "vsearch-small-abund" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}



	## misc

	# CD-HIT (-d 0 = use full sequence name; -T 1 = use one thread; -M 0 = no memory limit)
	RES=${OUTDIR}/cdhit_${T}.csv
	${LOG_CMD} "${CDHIT} -i ${INFILE} -o ${RES}.tmp -c 0.$((100 - ${T})) -d 0 -T 1 -M 0"
	awk ' BEGIN {FS = ">|(\\.\\.\\.)"}
		NR != 1 {printf /^>Cluster/ ? "\n" : $2" "}
	END {printf "\n"}' ${RES}.tmp.clstr | sed -e 's/ $//' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "cd-hit" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

	# DNACLUST (-t 1 = use one thread)
	RES=${OUTDIR}/dnaclust_${T}.csv
	${LOG_CMD} "${DNACLUST} -s 0.$((100 - ${T})) -i ${INFILE} -t 1 | sed 's/\t/ /g; s/ $//g' > ${RES}"
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "dnaclust" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES} ${CONF_TABLE}

	# Sumaclust
	RES=${OUTDIR}/sumaclust_${T}.csv
	${LOG_CMD} "${SUMACLUST} -t 0.$((100 - ${T})) -O ${RES}.tmp ${INFILE} > /dev/null"
	cut -d$'\t' -f2- ${RES}.tmp | sed 's/\t/ /g' > ${RES}
	for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
	do
		bash ${COMPUTE_METRICS} ${RES} "sumaclust" 0.$((100 - ${T})) ${TAXA_FILES[${I}]} ${CONF_TABLE} ${METRICS_FILES[${I}]} 0
	done
	rm ${RES}* ${CONF_TABLE}

done

rm ${INFILE_ALT} ${INFILE_ALT_LENGTH} ${INFILE_ALT_ABUNDANCE}


# evaluate quality
for ((I=0; I<${#GROUND_TRUTH_THRESHOLDS[@]}; I++))
do
	Rscript --vanilla scripts/eval_quality.R ${METRICS_FILES[${I}]} ${OUTDIR}/${DATA_SET}_${GROUND_TRUTH_THRESHOLDS[$[I]]}-metrics ${MIN_T} ${MAX_T}
done
