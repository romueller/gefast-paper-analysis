#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# Extracts the parameters from the executed GeFaST command
# and creates a new, extended (and more readable) data frame
prepare_gefast_data <- function(data) {
  
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- as.numeric(substr(data$cmd, regexpr("-t ", data$cmd) + 3, regexpr("--use-score", data$cmd) - 2))
  fastidious <- grepl("-sf", data$cmd)
  suppressWarnings(fastidious_threshold <- as.numeric(substr(data$cmd, 
                                                             regexpr("--swarm-fastidious-threshold", data$cmd) + 28, 
                                                             nchar(data$cmd))))
  fastidious_type <- rep("deactivated", nrow(data))
  fastidious_type[grepl("_f1", data$cmd)] <- "t + 1"
  fastidious_type[grepl("_2f", data$cmd)] <- "2 * t"
  mode <- grepl("--use-score", data$cmd)
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, fastidious_type, mode, percentage, repetition, stringsAsFactors = F)
  
}

# Extracts the parameters from the executed Swarm command
# and creates a new, extended (and more readable) data frame
prepare_swarm_data <- function(data) {
  
  cmd <- as.character(substr(data$cmd, regexpr("/", data$cmd) + 1, regexpr(" ", data$cmd) - 1))
  threshold <- as.numeric(substr(data$cmd, regexpr("-d ", data$cmd) + 3, regexpr("-a ", data$cmd) - 2))
  fastidious <- grepl("-f ", data$cmd)
  fastidious_threshold <- rep(NA, nrow(data))
  fastidious_threshold[fastidious] <- 2
  fastidious_type <- rep("deactivated", nrow(data))
  fastidious_type[fastidious] <- "2 * t"
  mode <- rep(T, nrow(data))
  subsample <- substr(data$cmd, regexpr("sub_", data$cmd) + 4, regexpr(".fasta", data$cmd) - 1)
  percentage <- as.numeric(substr(subsample, 0, regexpr("_", subsample) - 1))
  repetition <- as.numeric(substr(subsample, regexpr("_", subsample) + 1, nchar(subsample)))
  
  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, fastidious_type, mode, percentage, repetition, stringsAsFactors = F)
  
}

# Loads and prepares data using above functions.
# Creates visualisations of the development of runtime and 
# memory consumption with increasing data set size.
facet_point_plots <- function(log_file, plot_file, with_legend = T, with_error_bars = F, dev = "pdf") {

  # read & prepare data
  data <- read.csv(log_file, header = T, sep = ",", colClasses = c("numeric", "numeric", "character"))
  gefast_data <- prepare_gefast_data(data[grepl("GeFaST", data$cmd),])
  swarm_data <- prepare_swarm_data(data[!grepl("GeFaST", data$cmd),])
  data <- rbind(gefast_data, swarm_data)
  
  data$cmd[grepl("Swarm", data$cmd)] <- "Swarm (v2)"
  data$cmd <- factor(data$cmd, levels = c("Swarm (v2)", "GeFaST"))
  
  data$fastidious_type <- factor(data$fastidious_type, levels = c("deactivated", "t + 1", "2 * t"))
  
  # average over repetitions
  averaged_data <- data.frame(group_by(data, cmd, threshold, fastidious, fastidious_threshold, fastidious_type, mode, percentage) 
                              %>% summarise(ave_time = mean(time), sd_time = sd(time), ave_memory = mean(memory), sd_memory = sd(memory)))
  
  averaged_data$threshold_label <- paste0("Threshold t = ", averaged_data$threshold)
  
  # plot time measurements
  facets <- ggplot(averaged_data, aes(x = percentage, y = ave_time, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    xlab("Subset size [%]") +
    ylab("Average time [s]") +
    labs(colour = "Method", shape = "Fastidious") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_time - sd_time, ymax = ave_time + sd_time), width = 0.5)
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.1, 0.85))
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_time", ".", dev), width = 16, height = 9, device = dev)
  
  # plot memory measurements
  facets <- ggplot(averaged_data, aes(x = percentage, y = ave_memory, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    xlab("Subset size [%]") +
    ylab("Average memory [MiB]") +
    labs(colour = "Method", shape = "Fastidious") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_memory - sd_memory, ymax = ave_memory + sd_memory), width = 0.5)
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.1, 0.85))
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_memory", ".", dev), width = 16, height = 9, device = dev)
  
  # plot time & memory measurements
  facets1 <- ggplot(averaged_data, aes(x = percentage, y = ave_time, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    scale_x_continuous(labels = NULL) +
    xlab("") +
    ylab("Average time [s]") +
    labs(colour = "Method", shape = "Fastidious") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_error_bars) {
    facets1 <- facets1 + geom_errorbar(aes(ymin = ave_time - sd_time, ymax = ave_time + sd_time), width = 0.5)
  }
  
  if (with_legend) {
    facets1 <- facets1 + theme(legend.position = c(0.1, 0.75))
  } else {
    facets1 <- facets1 + theme(legend.position = "none")
  }
  
  facets2 <- ggplot(averaged_data, aes(x = percentage, y = ave_memory, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    xlab("Subset size [%]") +
    ylab("Average memory [MiB]") +
    labs(colour = "Method", shape = "Fastidious") +
    facet_grid(. ~ threshold_label) + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold"),
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
