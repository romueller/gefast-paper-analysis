#!/bin/bash

# Redraws the visualisations from the existing metrics / log files.



# ===== Clustering quality =====

GROUND_TRUTH_THRESHOLDS=(0.95 0.97 0.99)

for P in "${GROUND_TRUTH_THRESHOLDS[@]}"
do

	# even data set, threshold 1 - 10
	Rscript --vanilla scripts/eval_quality.R results/even_${P}-metrics.csv results/even_${P}-metrics 1 10

	# uneven data set, threshold 1 - 10
	Rscript --vanilla scripts/eval_quality.R results/uneven_${P}-metrics.csv results/uneven_${P}-metrics 1 10

done



# ===== Performance (runtime & memory) =====

# eldermet data set, threshold 1 - 10
Rscript --vanilla scripts/eval_performance.R results/eldermet-performance-measurements.csv results/eldermet-performance 1 10

# subsets of eldermet data set, threshold 1 - 2, 3 repetitions
Rscript --vanilla scripts/eval_subsampling.R results/eldermet-subsampling-measurements.csv results/eldermet-subsampling



# ===== Clustering quality & performance (runtime & memory) =====

# uneven data set, 80 % subsets, 97% ground truth, threshold 1 - 10, 5 repetitions
Rscript --vanilla scripts/eval_subsamples_fixed_quality.R results/uneven-sub-fixed-metrics.csv results/uneven-sub-fixed-quality
Rscript --vanilla scripts/eval_subsamples_fixed_performance.R results/uneven-sub-fixed-log.csv results/uneven-sub-fixed-performance

# even data set, 80 % subsets, 97% ground truth, threshold 1 - 10, 5 repetitions
Rscript --vanilla scripts/eval_subsamples_fixed_quality.R results/even-sub-fixed-metrics.csv results/even-sub-fixed-quality
Rscript --vanilla scripts/eval_subsamples_fixed_performance.R results/even-sub-fixed-log.csv results/even-sub-fixed-performance

# eldermet data set / genus level, 80 % subsets, 95% ground truth, threshold 1 - 10, 5 repetitions
Rscript --vanilla scripts/eval_subsamples_fixed_quality.R results/eldermet-sub-fixed-red-metrics.csv results/eldermet-sub-fixed-red-quality method,threshold,rep,recall,precision,adjrandindex
Rscript --vanilla scripts/eval_subsamples_fixed_performance.R results/eldermet-sub-fixed-red-log.csv results/eldermet-sub-fixed-red-performance

