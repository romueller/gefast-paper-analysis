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
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- rep("gefast-e", nrow(data))
  cmd[grepl("--use-score", data$cmd)] <- "gefast-s"
  threshold <- as.numeric(substr(data$cmd, regexpr("-t ", data$cmd) + 3, regexpr("--config ", data$cmd) - 2))
  fastidious <- grepl("-sf", data$cmd)
  suppressWarnings(fastidious_threshold <- as.numeric(substr(data$cmd, 
                                                             regexpr("--swarm-fastidious-threshold", data$cmd) + 28, 
                                                             nchar(data$cmd))))
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[grepl("_f1", data$cmd)] <- "-f1"
  submethod_type[grepl("_2f", data$cmd)] <- "-2f"

  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Swarm (v1) command
# and creates a new, extended (and more readable) data frame
prepare_swarm1_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- rep("swarm-v1", nrow(data))
  threshold <- as.numeric(substr(data$cmd, regexpr("-d ", data$cmd) + 3, regexpr("-o ", data$cmd) - 2))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("-nf", nrow(data))
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Swarm (v2) command
# and creates a new, extended (and more readable) data frame
prepare_swarm2_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- rep("swarm-v2", nrow(data))
  threshold <- as.numeric(substr(data$cmd, regexpr("-d ", data$cmd) + 3, regexpr("-o ", data$cmd) - 2))
  fastidious <- grepl("-f ", data$cmd)
  fastidious_threshold <- rep(NA, nrow(data))
  fastidious_threshold[fastidious] <- 2
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[fastidious] <- "-2f"
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed USEARCH command
# and creates a new, extended (and more readable) data frame
prepare_usearch_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- 100 * (1 - as.numeric(substr(data$cmd, regexpr("-id ", data$cmd) + 4, regexpr("-uc ", data$cmd) - 2)))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[grepl("-cluster_fast", data$cmd) & grepl("-sort length", data$cmd)] <- "-fast-length"
  submethod_type[grepl("-cluster_fast", data$cmd) & grepl("-sort size", data$cmd)] <- "-fast-abund"
  submethod_type[grepl("-cluster_smallmem", data$cmd) & grepl("-sortedby length", data$cmd)] <- "-small-length"
  submethod_type[grepl("-cluster_smallmem", data$cmd) & grepl("-sortedby size", data$cmd)] <- "-small-abund"
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed VSEARCH command
# and creates a new, extended (and more readable) data frame
prepare_vsearch_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- 100 * (1 - as.numeric(substr(data$cmd, regexpr("-id ", data$cmd) + 4, regexpr("-uc ", data$cmd) - 2)))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[grepl("-cluster_fast", data$cmd)] <- "-fast-length"
  submethod_type[grepl("-cluster_size", data$cmd)] <- "-size-abund"
  submethod_type[grepl("-cluster_smallmem", data$cmd) & grepl("_alt_length", data$cmd)] <- "-small-length"
  submethod_type[grepl("-cluster_smallmem", data$cmd) & grepl("_alt_abundance", data$cmd)] <- "-small-abund"
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed CD-HIT command
# and creates a new, extended (and more readable) data frame
prepare_cdhit_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- 100 * (1 - as.numeric(substr(data$cmd, regexpr("-c ", data$cmd) + 3, regexpr("-d ", data$cmd) - 2)))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("-nf", nrow(data))
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed DNACLUST command
# and creates a new, extended (and more readable) data frame
prepare_dnaclust_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- 100 * (1 - as.numeric(substr(data$cmd, regexpr("-s ", data$cmd) + 3, regexpr("-i ", data$cmd) - 2)))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("-nf", nrow(data))
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Sumaclust command
# and creates a new, extended (and more readable) data frame
prepare_sumaclust_data <- function(data) {
  
  data$cmd <- gsub("/bin/sh -c ", "", data$cmd)
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- 100 * (1 - as.numeric(substr(data$cmd, regexpr("-t ", data$cmd) + 3, regexpr("-O ", data$cmd) - 2)))
  fastidious <- rep(F, nrow(data))
  fastidious_threshold <- rep(NA, nrow(data))
  submethod_type <- rep("-nf", nrow(data))
  
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, submethod_type, percentage, repetition, stringsAsFactors = F)
  
}

# Loads and prepares data using above functions.
# Creates visualisations of the development of runtime and 
# memory consumption with increasing data set size.
facet_point_plots <- function(log_file, plot_file, with_legend = T, with_error_bars = F, dev = "pdf") {

  # read & prepare data
  data <- read.csv(log_file, header = T, sep = ",", colClasses = c("numeric", "numeric", "character"))
  gefast_data <- prepare_gefast_data(data[grepl("GeFaST", data$cmd),])
  swarm1_data <- prepare_swarm1_data(data[grepl("Swarm-1", data$cmd),])
  swarm2_data <- prepare_swarm2_data(data[grepl("Swarm-2", data$cmd),])
  usearch_data <- prepare_usearch_data(data[grepl("usearch", data$cmd),])
  vsearch_data <- prepare_vsearch_data(data[grepl("vsearch", data$cmd),])
  cdhit_data <- prepare_cdhit_data(data[grepl("cd-hit", data$cmd),])
  dnaclust_data <- prepare_dnaclust_data(data[grepl("dnaclust", data$cmd),])
  sumaclust_data <- prepare_sumaclust_data(data[grepl("sumaclust", data$cmd),])
  data <- rbind(gefast_data, swarm1_data, swarm2_data, usearch_data, vsearch_data, cdhit_data, dnaclust_data, sumaclust_data)

  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)", 
                    "USEARCH", "VSEARCH", "CD-HIT", "DNACLUST", "Sumaclust")
  data$cmd[grepl("swarm-v1", data$cmd)] <- method_levels[1]
  data$cmd[grepl("swarm-v2", data$cmd)] <- method_levels[2]
  data$cmd[grepl("gefast-e", data$cmd)] <- method_levels[3]
  data$cmd[grepl("gefast-s", data$cmd)] <- method_levels[4]
  data$cmd[grepl("usearch", data$cmd)] <- method_levels[5]
  data$cmd[grepl("vsearch", data$cmd)] <- method_levels[6]
  data$cmd[grepl("cd-hit", data$cmd)] <- method_levels[7]
  data$cmd[grepl("dnaclust", data$cmd)] <- method_levels[8]
  data$cmd[grepl("sumaclust", data$cmd)] <- method_levels[9]
  data$cmd <- factor(data$cmd, levels = method_levels)
  
  submethod_levels <- c("non-fastidious", "fastidious (t + 1)", "fastidious (2 * t)", 
                        "*SEARCH (fast, length)", "*SEARCH (smallmem, length)", 
                        "*SEARCH (fast, abundance)", "*SEARCH (smallmem, abundance)", "*SEARCH (size, abundance)")
  data$submethod_type[grepl("-nf", data$submethod_type)] <- submethod_levels[1] 
  data$submethod_type[grepl("-f1", data$submethod_type)] <- submethod_levels[2] 
  data$submethod_type[grepl("-2f", data$submethod_type)] <- submethod_levels[3]
  data$submethod_type[grepl("-fast-length", data$submethod_type)] <- submethod_levels[4]
  data$submethod_type[grepl("-small-length", data$submethod_type)] <- submethod_levels[5]
  data$submethod_type[grepl("-fast-abund", data$submethod_type)] <- submethod_levels[6]
  data$submethod_type[grepl("-small-abund", data$submethod_type)] <- submethod_levels[7]
  data$submethod_type[grepl("-size-abund", data$submethod_type)] <- submethod_levels[8]
  data$submethod_type <- factor(data$submethod_type, levels = submethod_levels)
  
  # average over repetitions
  averaged_data <- data.frame(group_by(data, cmd, threshold, fastidious, fastidious_threshold, submethod_type, percentage) 
                              %>% summarise(ave_time = mean(time), sd_time = sd(time), ave_memory = mean(memory), sd_memory = sd(memory)))
  
  # adjust palette & threshold labels depending on set of compared methods
  pal <- c(brewer.pal(8, "Set3"), "darkgrey")
  averaged_data$threshold_label <- paste0("Threshold t = ", averaged_data$threshold, " (", formatC(((100 - averaged_data$threshold) / 100), digits = 2, format = "f"), ")")
  if (all(data$cmd %in% method_levels[c(2, 4)])) { # only Swarm (v2) and GeFaST (scoring function)
    
    pal <- pal[c(1, 4)]
    averaged_data$threshold_label <- paste0("Threshold t = ", averaged_data$threshold)
    
  } else if (all(data$cmd %in% method_levels[c(2, 4, 5, 6, 7, 8, 9)])) { # all except Swarm (v1) and GeFaST (edit distance)
    pal <- pal[c(2, 4, 5, 6, 7, 8, 9)]
  }
  
  # plot time measurements
  facets <- ggplot(averaged_data, aes(x = percentage, y = ave_time, group = cmd, colour = cmd, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_colour_manual(values = pal) +
    xlab("Subset size [%]") +
    ylab("Average time [s]") +
    labs(colour = "Method", shape = "Submethod") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_time - sd_time, ymax = ave_time + sd_time), width = 0.5)
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.2, 0.85), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_time", ".", dev), width = 16, height = 9, device = dev)
  
  # plot memory measurements
  facets <- ggplot(averaged_data, aes(x = percentage, y = ave_memory, group = cmd, colour = cmd, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_colour_manual(values = pal) +
    xlab("Subset size [%]") +
    ylab("Average memory [MiB]") +
    labs(colour = "Method", shape = "Submethod") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_memory - sd_memory, ymax = ave_memory + sd_memory), width = 0.5)
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.2, 0.85), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_memory", ".", dev), width = 16, height = 9, device = dev)
  
  # plot time & memory measurements
  facets1 <- ggplot(averaged_data, aes(x = percentage, y = ave_time, group = cmd, colour = cmd, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_colour_manual(values = pal) +
    scale_x_continuous(labels = NULL) +
    xlab("") +
    ylab("Average time [s]") +
    labs(colour = "Method", shape = "Submethod") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets1 <- facets1 + geom_errorbar(aes(ymin = ave_time - sd_time, ymax = ave_time + sd_time), width = 0.5)
  }
  
  if (with_legend) {
    facets1 <- facets1 + theme(legend.position = c(0.2, 0.75), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets1 <- facets1 + theme(legend.position = "none")
  }
  
  facets2 <- ggplot(averaged_data, aes(x = percentage, y = ave_memory, group = cmd, colour = cmd, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_colour_manual(values = pal) +
    xlab("Subset size [%]") +
    ylab("Average memory [MiB]") +
    labs(colour = "Method", shape = "Submethod") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold"),
      legend.position = "none"
    )
  
  if (with_error_bars) {
    facets2 <- facets2 + geom_errorbar(aes(ymin = ave_memory - sd_memory, ymax = ave_memory + sd_memory), width = 0.5)
  }
  
  ggsave(grid.draw(rbind(ggplotGrob(facets1), ggplotGrob(facets2), size = "last")), file = paste0(plot_file, "_combined", ".", dev), width = 16, height = 11, device = dev)
  
}


# get arguments & create visualisation
args = commandArgs(trailingOnly = T)

metrics_file <- args[1]
plot_file <- args[2]

facet_point_plots(metrics_file, paste0(plot_file, "_with_legend"), T, T)
facet_point_plots(metrics_file, paste0(plot_file, "_without_legend"), F, T)
