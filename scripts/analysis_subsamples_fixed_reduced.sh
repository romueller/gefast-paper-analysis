#!/bin/bash

# Adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
# 
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
#
# Performs a performance and clustering-quality analysis on the given data set
# using a ground truth and multiple random subsamples of a specified size



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
COMPUTE_METRICS=scripts/compute_metrics_reduced.sh


DATA_SET=$1
MIN_T=$2
MAX_T=$3
GROUND_TRUTH_THRESHOLD=$4
PROPORTION=$5
REPETITIONS=$6

OUTDIR=results
CONF_TABLE=${OUTDIR}/confusion-table.conf
FASTA_FILE=data/${DATA_SET}/${DATA_SET}_derep.fasta

SUBSAMPLE_STEM=data/${DATA_SET}/${DATA_SET}_derep_sub
REMOVE_SUBSAMPLES=1



# ===== Generate subsamples and ground truth =====

# reference configuration
REP_SET=data/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta
TAXONOMY=data/SILVA_128_QIIME_release/taxonomy/16S_only/97/majority_taxonomy_7_levels.txt
REFORMAT_SCRIPT=scripts/reformat-silva-taxonomy.R

# assign taxa (with VSEARCH to avoid memory limit of 32-bit USEARCH)
# and replace identifiers of reference representatives with the actual taxonomy information
TAXA_FILE=data/${DATA_SET}/${DATA_SET}_derep.spec.ualn
bash scripts/assign_taxa.sh ${VSEARCH} 8 ${FASTA_FILE} ${REP_SET} ${TAXA_FILE}.tmp ${GROUND_TRUTH_THRESHOLD}
python scripts/assign_taxa_rep.py -r ${TAXONOMY} -t ${TAXA_FILE}.tmp > ${TAXA_FILE}
rm ${TAXA_FILE}.tmp

# reduce data to sequences with clean, "complete" taxonomy
RED_FASTA_FILE=data/${DATA_SET}/${DATA_SET}_derep.reduced.fasta
RED_TAXA_FILE=data/${DATA_SET}/${DATA_SET}_derep.reduced.spec.ualn
RED_IDS_FILE=data/${DATA_SET}/${DATA_SET}_derep.reduced.ids
Rscript --vanilla ${REFORMAT_SCRIPT} ${TAXA_FILE} ${RED_TAXA_FILE}
cut -d$'\t' -f1 ${RED_TAXA_FILE} > ${RED_IDS_FILE}
./tools/FastaSampler -l ${FASTA_FILE} ${RED_IDS_FILE} ${RED_FASTA_FILE}

# generate subsamples
SUB_SIZE=$((($(wc -l ${RED_TAXA_FILE} | cut -d' ' -f1) * ${PROPORTION}) / 100))
for ((R=0; R<${REPETITIONS}; R++))
do

	SUB_RED_TAXA_FILE=${SUBSAMPLE_STEM}_${PROPORTION}_${R}.reduced.spec.ualn
	SUB_RED_IDS_FILE=${SUBSAMPLE_STEM}_${PROPORTION}_${R}.reduced.ids
	SUB_RED_FASTA_FILE=${SUBSAMPLE_STEM}_${PROPORTION}_${R}.reduced.fasta

	shuf -n ${SUB_SIZE} ${RED_TAXA_FILE} > ${SUB_RED_TAXA_FILE}
	cut -d$'\t' -f1 ${SUB_RED_TAXA_FILE} > ${SUB_RED_IDS_FILE}
	./tools/FastaSampler -l ${RED_FASTA_FILE} ${SUB_RED_IDS_FILE} ${SUB_RED_FASTA_FILE}

done



# ===== Clustering & confusion table =====

# prepare metrics & performance-log file
METRICS_FILE=${OUTDIR}/${DATA_SET}-sub-fixed-red-metrics.csv
echo "method,threshold,rep,recall,precision,adjrandindex" > ${METRICS_FILE}

LOG_FILE=${OUTDIR}/${DATA_SET}-sub-fixed-red-log.csv
LOG_CMD="/usr/bin/time -f %e,%M,%C -a -o ${LOG_FILE} /bin/sh -c "
echo "time,memory,cmd" > ${LOG_FILE}


# run tools, compute & evaluate confusion table
for ((T=${MIN_T}; T<=${MAX_T}; T++))
do
	for ((R=0; R<${REPETITIONS}; R++))
	do

		INFILE=${SUBSAMPLE_STEM}_${PROPORTION}_${R}.reduced.fasta
		INFILE_ALT=${INFILE}_alt
		INFILE_ALT_LENGTH=${INFILE}_alt_length
		INFILE_ALT_ABUNDANCE=${INFILE}_alt_abundance
		TAXA_FILE=${SUBSAMPLE_STEM}_${PROPORTION}_${R}.reduced.spec.ualn

		# create input files for USEARCH and VSEARCH (with VSEARCH to avoid memory limit of 32-bit USEARCH)
		sed 's/_/;size=/' ${INFILE} > ${INFILE_ALT} 
		${VSEARCH} -sortbylength ${INFILE_ALT} -output ${INFILE_ALT_LENGTH} -minsize 1
		${VSEARCH} -sortbysize ${INFILE_ALT} -output ${INFILE_ALT_ABUNDANCE} -minsize 1



		## Swarm

#		# Swarm (version 1, non-fastidious)
#		RES=${OUTDIR}/swarm-v1_o_${T}.csv
#		${LOG_CMD} "${SWARM_V1} -d ${T} -o ${RES}.tmp ${INFILE}; \
#		python ${SWARM_BREAKER} -b ${SWARM_V1} -f ${INFILE} -s ${RES}.tmp -d ${T} > ${RES}"
#		rm ${RES}.tmp
#		bash ${COMPUTE_METRICS} ${RES} "swarm-v1" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES} ${CONF_TABLE}
#
#		# Swarm (version 2, non-fastidious)
#		RES=${OUTDIR}/swarm-v2_o_${T}.csv
#		${LOG_CMD} "${SWARM_V2} -d ${T} -o ${RES} -a 1 ${INFILE}" 
#		bash ${COMPUTE_METRICS} ${RES} "swarm-v2" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES} ${CONF_TABLE}
#
#		# Swarm (version 2, fastidious)
#		if [ "${T}" == "1" ]; then
#			RES=${OUTDIR}/swarm-v2_s_${T}_f.csv
#			${LOG_CMD} "${SWARM_V2} -f -d ${T} -o ${RES} -a 1 ${INFILE}"
#			bash ${COMPUTE_METRICS} ${RES} "swarm-v2-2f" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#			rm ${RES} ${CONF_TABLE}
#		fi



		## GeFaST

#		# GeFaST (edit distance, non-fastidious)
#		RES=${OUTDIR}/gefast-e_o_${T}_e.csv
#		${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES}"
#		bash ${COMPUTE_METRICS} ${RES} "gefast-e" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES} ${CONF_TABLE}
#
#		# GeFaST (edit distance, fastidious, t + 1)
#		RES=${OUTDIR}/gefast-e_o_${T}_ef1.csv
#		${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} -sf --swarm-fastidious-threshold $((${T} + 1))"
#		bash ${COMPUTE_METRICS} ${RES} "gefast-e-f1" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES} ${CONF_TABLE}
#
#		# GeFaST (edit distance, fastidious, 2 * t)
#		RES=${OUTDIR}/gefast-e_o_${T}_e2f.csv
#		${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} -sf --swarm-fastidious-threshold $((2 * ${T}))"
#		bash ${COMPUTE_METRICS} ${RES} "gefast-e-2f" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES} ${CONF_TABLE}


		# GeFaST (scoring function, non-fastidious)
		RES=${OUTDIR}/gefast-s_o_${T}_s.csv
		${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} --use-score"
		bash ${COMPUTE_METRICS} ${RES} "gefast-s" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES} ${CONF_TABLE}

		# GeFaST (scoring function, fastidious, t + 1)
		RES=${OUTDIR}/gefast-s_o_${T}_sf1.csv
		${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} --use-score -sf --swarm-fastidious-threshold $((${T} + 1))"
		bash ${COMPUTE_METRICS} ${RES} "gefast-s-f1" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES} ${CONF_TABLE}

		# GeFaST (scoring function, fastidious, 2 * t)
		RES=${OUTDIR}/gefast-s_o_${T}_s2f.csv
		${LOG_CMD} "${GEFAST} ${INFILE} -t ${T} --config ${CONF} -so ${RES} --use-score -sf --swarm-fastidious-threshold $((2 * ${T}))"
		bash ${COMPUTE_METRICS} ${RES} "gefast-s-2f" ${T} ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES} ${CONF_TABLE}



		## USEARCH

#		# USEARCH (cluster_fast, sort by length)
#		RES=${OUTDIR}/usearch-fast-length_${T}.csv
#		${LOG_CMD} "${USEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sort length -fulldp"
#		awk 'BEGIN {FS="\t"}{
#			if ($1 == "S") clusters[$2] = $9
#			if ($1 == "H") clusters[$2] = clusters[$2] " " $9
#		}
#		END {
#			for (cluster in clusters) {
#				print clusters[cluster]
#			}
#		}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
#		bash ${COMPUTE_METRICS} ${RES} "usearch-fast-length" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES}* ${CONF_TABLE}

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
		bash ${COMPUTE_METRICS} ${RES} "usearch-fast-abund" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES}* ${CONF_TABLE}

#		# USEARCH (cluster_smallmem, presorted by length)
#		RES=${OUTDIR}/usearch-small-length_${T}.csv
#		${LOG_CMD} "${USEARCH} -cluster_smallmem ${INFILE_ALT_LENGTH} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sortedby length -fulldp"
#		awk 'BEGIN {FS="\t"}{
#			if ($1 == "S") clusters[$2] = $9
#			if ($1 == "H") clusters[$2] = clusters[$2] " " $9
#		}
#		END {
#			for (cluster in clusters) {
#				print clusters[cluster]
#			}
#		}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
#		bash ${COMPUTE_METRICS} ${RES} "usearch-small-length" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES}* ${CONF_TABLE}

#		# USEARCH (cluster_smallmem, presorted by abundance)
#		RES=${OUTDIR}/usearch-small-abund_${T}.csv
#		${LOG_CMD} "${USEARCH} -cluster_smallmem ${INFILE_ALT_ABUNDANCE} -id 0.$((100 - ${T})) -uc ${RES}.tmp -sortedby size -fulldp"
#		awk 'BEGIN {FS="\t"}{
#			if ($1 == "S") clusters[$2] = $9
#			if ($1 == "H") clusters[$2] = clusters[$2] " " $9
#		}
#		END {
#			for (cluster in clusters) {
#				print clusters[cluster]
#			}
#		}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
#		bash ${COMPUTE_METRICS} ${RES} "usearch-small-abund" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES}* ${CONF_TABLE}



		## VSEARCH

#		# VSEARCH (cluster_fast, sort by length)
#		RES=${OUTDIR}/vsearch-fast-length_${T}.csv
#		${LOG_CMD} "${VSEARCH} -threads 1 -cluster_fast ${INFILE_ALT} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1"
#		awk 'BEGIN {FS="\t"}{
#			if ($1 == "S") clusters[$2] = $9
#			if ($1 == "H") clusters[$2] = clusters[$2] " " $9
#		}
#		END {
#			for (cluster in clusters) {
#				print clusters[cluster]
#			}
#		}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
#		bash ${COMPUTE_METRICS} ${RES} "vsearch-fast-length" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES}* ${CONF_TABLE}

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
		bash ${COMPUTE_METRICS} ${RES} "vsearch-size-abund" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES}* ${CONF_TABLE}

#		# VSEARCH (cluster_smallmem, presorted by length)
#		RES=${OUTDIR}/vsearch-small-length_${T}.csv
#		${LOG_CMD} "${VSEARCH} -threads 1 -cluster_smallmem ${INFILE_ALT_LENGTH} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1"
#		awk 'BEGIN {FS="\t"}{
#			if ($1 == "S") clusters[$2] = $9
#			if ($1 == "H") clusters[$2] = clusters[$2] " " $9
#		}
#		END {
#			for (cluster in clusters) {
#				print clusters[cluster]
#			}
#		}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
#		bash ${COMPUTE_METRICS} ${RES} "vsearch-small-length" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES}* ${CONF_TABLE}

#		# VSEARCH (cluster_smallmem, presorted by abundance)
#		RES=${OUTDIR}/vsearch-small-abund_${T}.csv
#		${LOG_CMD} "${VSEARCH} -threads 1 -cluster_smallmem ${INFILE_ALT_ABUNDANCE} -id 0.$((100 - ${T})) -uc ${RES}.tmp -minsize 1 -usersort"
#		awk 'BEGIN {FS="\t"}{
#			if ($1 == "S") clusters[$2] = $9
#			if ($1 == "H") clusters[$2] = clusters[$2] " " $9
#		}
#		END {
#			for (cluster in clusters) {
#				print clusters[cluster]
#			}
#		}' ${RES}.tmp | sed 's/;size=/_/g' > ${RES}
#		bash ${COMPUTE_METRICS} ${RES} "vsearch-small-abund" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
#		rm ${RES}* ${CONF_TABLE}



		## misc

		# CD-HIT (-d 0 = use full sequence name; -T 1 = use one thread; -M 0 = no memory limit)
		RES=${OUTDIR}/cdhit_${T}.csv
		${LOG_CMD} "${CDHIT} -i ${INFILE} -o ${RES}.tmp -c 0.$((100 - ${T})) -d 0 -T 1 -M 0"
		awk ' BEGIN {FS = ">|(\\.\\.\\.)"}
		NR != 1 {printf /^>Cluster/ ? "\n" : $2" "}
		END {printf "\n"}' ${RES}.tmp.clstr | sed -e 's/ $//' > ${RES}
		bash ${COMPUTE_METRICS} ${RES} "cd-hit" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES}* ${CONF_TABLE}

		# DNACLUST (-t 1 = use one thread)
		RES=${OUTDIR}/dnaclust_${T}.csv
		${LOG_CMD} "${DNACLUST} -s 0.$((100 - ${T})) -i ${INFILE} -t 1 | sed 's/\t/ /g; s/ $//g' > ${RES}"
		bash ${COMPUTE_METRICS} ${RES} "dnaclust" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES} ${CONF_TABLE}

		# Sumaclust
		RES=${OUTDIR}/sumaclust_${T}.csv
		${LOG_CMD} "${SUMACLUST} -t 0.$((100 - ${T})) -O ${RES}.tmp ${INFILE} > /dev/null"
		cut -d$'\t' -f2- ${RES}.tmp | sed 's/\t/ /g' > ${RES}
		bash ${COMPUTE_METRICS} ${RES} "sumaclust" 0.$((100 - ${T})) ${TAXA_FILE} ${CONF_TABLE} ${METRICS_FILE} ${R}
		rm ${RES}* ${CONF_TABLE}



		rm ${INFILE_ALT} ${INFILE_ALT_LENGTH} ${INFILE_ALT_ABUNDANCE}

	done
done

# evaluate quality and performance
Rscript --vanilla scripts/eval_subsamples_fixed_quality.R ${METRICS_FILE} ${OUTDIR}/${DATA_SET}-sub-fixed-red-quality method,threshold,rep,recall,precision,adjrandindex
Rscript --vanilla scripts/eval_subsamples_fixed_performance.R ${LOG_FILE} ${OUTDIR}/${DATA_SET}-sub-fixed-red-performance

# clean up
if [ "${REMOVE_SUBSAMPLES}" == "1" ]; then
	rm ${SUBSAMPLE_STEM}*
fi
