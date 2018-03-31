#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(reshape2)

# Extracts the parameters from the executed GeFaST command
# and creates a new, extended (and more readable) data frame
prepare_gefast_data <- function(data) {
  
  method <- rep("gefast-e", nrow(data))
  method[grepl("gefast-s", data$method)] <- "gefast-s"
  threshold <- data$threshold
  
  submethod_type <- rep("nf", nrow(data))
  submethod_type[grepl("-[es]-f1", data$method)] <- "f1"
  submethod_type[grepl("-[es]-2f", data$method)] <- "2f"
  
  fastidious <- submethod_type != "-nf"
  fastidious_threshold <- rep(NA, nrow(data))
  fastidious_threshold[submethod_type == "-f1"] <- data$threshold[submethod_type == "-f1"] + 1
  fastidious_threshold[submethod_type == "-2f"] <- 2 * data$threshold[submethod_type == "-2f"]
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Swarm (v1) command
# and creates a new, extended (and more readable) data frame
prepare_swarm1_data <- function(data) {
  
  method <- rep("swarm-v1", nrow(data))
  threshold <- data$threshold
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("nf", nrow(data))
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Swarm (v2) command
# and creates a new, extended (and more readable) data frame
prepare_swarm2_data <- function(data) {
  
  method <- rep("swarm-v2", nrow(data))
  threshold <- data$threshold
  fastidious <- grepl("-2f", data$method)
  fastidious_threshold <- rep(NA, nrow(data))
  fastidious_threshold[fastidious] <- 2
  submethod_type <- rep("nf", nrow(data))
  submethod_type[fastidious] <- "2f"
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed USEARCH command
# and creates a new, extended (and more readable) data frame
prepare_usearch_data <- function(data) {
  
  method <- rep("usearch", nrow(data))
  threshold <- round(100 * (1 - data$threshold))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("nf", nrow(data))
  submethod_type[grepl("-fast-length", data$method)] <- "fast-length"
  submethod_type[grepl("-fast-abund", data$method)] <- "fast-abund"
  submethod_type[grepl("-small-length", data$method)] <- "small-length"
  submethod_type[grepl("-small-abund", data$method)] <- "small-abund"
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed VSEARCH command
# and creates a new, extended (and more readable) data frame
prepare_vsearch_data <- function(data) {
  
  method <- rep("vsearch", nrow(data))
  threshold <- round(100 * (1 - data$threshold))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("nf", nrow(data))
  submethod_type[grepl("-fast-length", data$method)] <- "fast-length"
  submethod_type[grepl("-size-abund", data$method)] <- "size-abund"
  submethod_type[grepl("-small-length", data$method)] <- "small-length"
  submethod_type[grepl("-small-abund", data$method)] <- "small-abund"
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed CD-HIT command
# and creates a new, extended (and more readable) data frame
prepare_cdhit_data <- function(data) {
  
  method <- rep("cd-hit", nrow(data))
  threshold <- round(100 * (1 - data$threshold))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("nf", nrow(data))
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed DNACLUST command
# and creates a new, extended (and more readable) data frame
prepare_dnaclust_data <- function(data) {
  
  method <- rep("dnaclust", nrow(data))
  threshold <- round(100 * (1 - data$threshold))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("nf", nrow(data))
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Sumaclust command
# and creates a new, extended (and more readable) data frame
prepare_sumaclust_data <- function(data) {
  
  method <- rep("sumaclust", nrow(data))
  threshold <- round(100 * (1 - data$threshold))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("nf", nrow(data))
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Reads and prepares the data for subsequent significance tests.
read_data <- function(metrics_file, col_names = c("method", "threshold", "rep", "recall", "precision", 
                                                    "nmi", "randindex", "adjrandindex"), 
                      dev = "pdf") {
  
  data <- read.csv(metrics_file, header = T, sep = ",", colClasses = c("character", rep("numeric", length(col_names) - 1)))
  colnames(data) <- col_names
  
  gefast_data <- prepare_gefast_data(data[grepl("gefast", data$method),])
  swarm1_data <- prepare_swarm1_data(data[grepl("swarm-v1", data$method),])
  swarm2_data <- prepare_swarm2_data(data[grepl("swarm-v2", data$method),])
  usearch_data <- prepare_usearch_data(data[grepl("usearch", data$method),])
  vsearch_data <- prepare_vsearch_data(data[grepl("vsearch", data$method),])
  cdhit_data <- prepare_cdhit_data(data[grepl("cd-hit", data$method),])
  dnaclust_data <- prepare_dnaclust_data(data[grepl("dnaclust", data$method),])
  sumaclust_data <- prepare_sumaclust_data(data[grepl("sumaclust", data$method),])
  data <- rbind(gefast_data, swarm1_data, swarm2_data, usearch_data, vsearch_data, cdhit_data, dnaclust_data, sumaclust_data)
  colnames(data) <- c("method", "threshold", "fastidious", "fastidious_threshold", "submethod_type",
                      "repetition", "recall", "precision", "adjrandindex")

  
  data
  
}

# Helper method to easily subset the data
select_data <- function(data, sel_method, sel_threshold, sel_submethod) {
  subset(data, method == sel_method & threshold == sel_threshold & submethod_type == sel_submethod)
}

# Selects specified parts of the data and performs a test of statistical significance
compare <- function(data, method1, threshold1, submethod1, method2, threshold2, submethod2) {
  
  data_method1 <- select_data(data, method1, threshold1, submethod1)
  data_method2 <- select_data(data, method2, threshold2, submethod2)
  
  metrics <- c("recall", "precision", "adjrandindex")
  cols <- c("method1", "threshold1", "submethod1", "method2", "threshold2", "submethod2", "metric", 
            "mean_diff", "sd_diff", "p_value", "power1", "power2")
  
  
  res <- matrix(nrow = length(metrics), ncol = length(cols))
  
  for (i in 1:length(metrics)) {

    v1 <- data_method1[, metrics[i]]
    v2 <- data_method2[, metrics[i]]
    test_res <- t.test(v1, v2, var.equal = F, paired = T)
    res[i, ] <- c(method1, threshold1, submethod1, method2, threshold2, submethod2, metrics[i],
                  mean(v1 - v2), sd(v1 - v2), test_res$p.value, 
                  abs(mean(v1 - v2) / sd(v1 - v2)), abs(mean(v1 - v2) / ((mean(v1) + mean(v2)) / 2)))
    
  }
  
  res <- data.frame(res)
  colnames(res) <- cols
  
  res
  
}

# Performs tests of statistical significance between the methods specified in comps and targets
compute_significances <- function(metrics_file, ouput_file, thresholds, targets, comps, 
                                  col_names = c("method", "threshold", "rep", "recall", "precision", "nmi", "randindex", "adjrandindex")) {
  
  data <- read_data(metrics_file, col_names)
  
  res <- data.frame(method1 = as.factor(character()), threshold1 = integer(), submethod1 = as.factor(character()),
                    method2 = as.factor(character()), threshold2 = integer(), submethod2 = as.factor(character()),
                    metric = as.factor(character()), mean_diff = numeric(), sd_diff = numeric(),
                    p_value = numeric(), power1 = numeric(), power2 = numeric()) 
  
  for (i in 1:nrow(comps)) {
    for (j in 1:nrow(targets)) {
      for (t in thresholds) {
        
        if (!(comps$method[i] == "swarm-v2" & comps$submethod[i] == "2f") 
            & !(targets$method[j] == "swarm-v2" & targets$submethod[j] == "2f") 
            | t == 1) {
          res <- rbind(res, compare(data, comps$method[i], t, comps$submethod[i], targets$method[j], t, targets$submethod[j]))
        }
        
      }
    }
  }
  
  write.table(res, output_file, sep = ",", row.names = F, quote = F)
  
  res
  
}

## Analyses

# even
metrics_file <- "results/even-sub-fixed-metrics.csv"
output_file <- "results/even-sub-fixed-significances.csv"

thresholds <- 1:10

targets <- data.frame(rbind(
  c("gefast-s", "nf"), 
  c("gefast-s", "f1"), 
  c("gefast-s", "2f")
))
names(targets) <- c("method", "submethod")
targets$method <- as.character(targets$method)
targets$submethod <- as.character(targets$submethod)

comps <- data.frame(rbind( 
               c("swarm-v1", "nf"),
               c("swarm-v2", "nf"),
               c("swarm-v2", "2f"),
               c("gefast-e", "nf"),
               c("gefast-e", "f1"),
               c("gefast-e", "2f"),
               c("gefast-s", "nf"),
               c("gefast-s", "f1"),
               c("gefast-s", "2f"),
               c("usearch", "fast-abund"),
               c("usearch", "fast-length"),
               c("usearch", "small-abund"),
               c("usearch", "small-length"),
               c("vsearch", "size-abund"),
               c("vsearch", "fast-length"),
               c("vsearch", "small-abund"),
               c("vsearch", "small-length"),
               c("cd-hit", "nf"),
               c("dnaclust", "nf"),
               c("sumaclust", "nf")
))
names(comps) <- c("method", "submethod")
comps$method <- as.character(comps$method)
comps$submethod <- as.character(comps$submethod)

compute_significances(metrics_file, output_file, thresholds, targets, comps)



# uneven
metrics_file <- "results/uneven-sub-fixed-metrics.csv"
output_file <- "results/uneven-sub-fixed-significances.csv"

thresholds <- 1:10

targets <- data.frame(rbind(
  c("gefast-s", "nf"), 
  c("gefast-s", "f1"), 
  c("gefast-s", "2f")
))
names(targets) <- c("method", "submethod")
targets$method <- as.character(targets$method)
targets$submethod <- as.character(targets$submethod)

comps <- data.frame(rbind( 
  c("swarm-v1", "nf"),
  c("swarm-v2", "nf"),
  c("swarm-v2", "2f"),
  c("gefast-e", "nf"),
  c("gefast-e", "f1"),
  c("gefast-e", "2f"),
  c("gefast-s", "nf"),
  c("gefast-s", "f1"),
  c("gefast-s", "2f"),
  c("usearch", "fast-abund"),
  c("usearch", "fast-length"),
  c("usearch", "small-abund"),
  c("usearch", "small-length"),
  c("vsearch", "size-abund"),
  c("vsearch", "fast-length"),
  c("vsearch", "small-abund"),
  c("vsearch", "small-length"),
  c("cd-hit", "nf"),
  c("dnaclust", "nf"),
  c("sumaclust", "nf")
))
names(comps) <- c("method", "submethod")
comps$method <- as.character(comps$method)
comps$submethod <- as.character(comps$submethod)

compute_significances(metrics_file, output_file, thresholds, targets, comps)



# reduced eldermet
metrics_file <- "results/eldermet-sub-fixed-red-metrics.csv"
output_file <- "results/eldermet-sub-fixed-red-significances.csv"

thresholds <- 1:10

targets <- data.frame(rbind(
  c("gefast-s", "nf"), 
  c("gefast-s", "f1"), 
  c("gefast-s", "2f")
))
names(targets) <- c("method", "submethod")
targets$method <- as.character(targets$method)
targets$submethod <- as.character(targets$submethod)

comps <- data.frame(rbind( 
  c("gefast-s", "nf"),
  c("gefast-s", "f1"),
  c("gefast-s", "2f"),
  c("usearch", "fast-abund"),
  c("vsearch", "size-abund"),
  c("cd-hit", "nf"),
  c("dnaclust", "nf"),
  c("sumaclust", "nf")
))
names(comps) <- c("method", "submethod")
comps$method <- as.character(comps$method)
comps$submethod <- as.character(comps$submethod)

compute_significances(metrics_file, output_file, thresholds, targets, comps,
                      col_names = c("method", "threshold", "rep", "recall", "precision", "adjrandindex"))

