#!/bin/bash

# Runs the complete analysis for the paper: 
# downloading and preparing the data,
# running the tools, and visualising the results.



# ===== Preparations =====

# get & prepare sequence data sets
bash scripts/prepare.sh

# get reference data for (un)even
wget -nc -P data/ https://raw.githubusercontent.com/torognes/vsearch-eval/master/cluster/data/rrna_reference.fasta

# get reference data for eldermet
wget -nc -P data/ https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz
tar -zxvf data/Silva_128_release.tgz
rm data/Silva_128_release.tgz


# ===== Clustering quality =====

# even data set, threshold 1 - 10
bash scripts/analysis_quality.sh even data/rrna_reference.fasta 1 10

# uneven data set, threshold 1 - 10
bash scripts/analysis_quality.sh uneven data/rrna_reference.fasta 1 10


# ===== Performance (runtime & memory) =====

# eldermet data set, threshold 1 - 10
bash scripts/analysis_performance.sh eldermet 1 10

# subsets of eldermet data set, threshold 1 - 2, 3 repetitions
bash scripts/analysis_subsampling.sh eldermet 1 2 3


# ===== Clustering quality & performance (runtime & memory) =====

# uneven data set, 80 % subsets, 97% ground truth, threshold 1 - 10, 5 repetitions
bash scripts/analysis_subsamples_fixed.sh uneven data/rrna_reference.fasta 1 10 0.97 80 5

# even data set, 80 % subsets, 97% ground truth, threshold 1 - 10, 5 repetitions
bash scripts/analysis_subsamples_fixed.sh even data/rrna_reference.fasta 1 10 0.97 80 5

# eldermet data set / genus level, 80 % subsets, 95% ground truth, threshold 1 - 10, 5 repetitions
bash scripts/analysis_subsamples_fixed_reduced.sh eldermet 1 10 0.95 80 5



