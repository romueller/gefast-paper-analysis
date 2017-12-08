#!/bin/bash

# Redraws the visualisations from the existing metrics / log files.



# ===== Clustering quality =====

# even data set, threshold 1 - 10
Rscript --vanilla scripts/eval_quality.R results/even-metrics.csv results/even-metrics 1 10

# uneven data set, threshold 1 - 10
Rscript --vanilla scripts/eval_quality.R results/uneven-metrics.csv results/uneven-metrics 1 10


# ===== Performance (runtime & memory) =====

# eldermet data set, threshold 1 - 10
Rscript --vanilla scripts/eval_performance.R results/eldermet-performance-measurements.csv results/eldermet-performance 1 10

# subsets of eldermet data set, threshold 1 - 2
Rscript --vanilla scripts/eval_subsampling.R results/eldermet-subsampling-measurements.csv results/eldermet-subsampling

