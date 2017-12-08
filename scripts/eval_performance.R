#!/usr/bin/env Rscript

library(ggplot2)
library(grid)
library(gridExtra)

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

  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, fastidious_type, stringsAsFactors = F)
  
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

  time <- data$time
  memory <- data$memory / 1024
  
  data.frame(cmd, time, memory, threshold, fastidious, fastidious_threshold, fastidious_type, stringsAsFactors = F)
  
}

# Loads and prepares data using above functions.
# Creates visualisations of the development of runtime and 
# memory consumption with increasing threshold.
facet_point_plots <- function(metrics_file, plot_file, thresholds, with_legend = T, dev = "pdf") {

  data <- read.csv(metrics_file, header = T, sep = ",", colClasses = c("numeric", "numeric", "character"))
  gefast_data <- prepare_gefast_data(data[grepl("GeFaST", data$cmd),])
  swarm_data <- prepare_swarm_data(data[!grepl("GeFaST", data$cmd),])
  data <- rbind(gefast_data, swarm_data)

  data$cmd[grepl("Swarm", data$cmd)] <- "Swarm (v2)"
  data$cmd <- factor(data$cmd, levels = c("Swarm (v2)", "GeFaST"))
  
  data$fastidious_type <- factor(data$fastidious_type, levels = c("deactivated", "t + 1", "2 * t"))

  
  # plot time measurements
  facet <- ggplot(data, aes(x = threshold, y = time, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) +  
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL) +
    xlab("Threshold t") +
    ylab("Time [s]") +
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_legend) {
    facet <- facet + theme(legend.position = c(0.08, 0.85), legend.direction = "vertical", legend.box = "vertical")
  } else {
    facet <- facet + theme(legend.position = "none")
  }
  
  ggsave(facet, file = paste0(plot_file, "_time", ".", dev), width = 16, height = 9, device = dev)
  
  # plot memory measurements
  facet <- ggplot(data, aes(x = threshold, y = memory, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) +  
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL) +
    xlab("Threshold t") +
    ylab("Memory [MiB]") +
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_legend) {
    facet <- facet + theme(legend.position = c(0.08, 0.85), legend.direction = "vertical", legend.box = "vertical")
  } else {
    facet <- facet + theme(legend.position = "none")
  }
  
  ggsave(facet, file = paste0(plot_file, "_memory", ".", dev), width = 16, height = 9, device = dev)
  
  # plot time & memory measurements
  facet1 <- ggplot(data, aes(x = threshold, y = time, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = NULL) +
    xlab("") +
    ylab("Time [s]") +
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_legend) {
    facet1 <- facet1 + theme(legend.position = c(0.08, 0.75), legend.direction = "vertical", legend.box = "vertical")
  } else {
    facet1 <- facet1 + theme(legend.position = "none")
  }
  
  facet2 <- ggplot(data, aes(x = threshold, y = memory, group = cmd, colour = cmd, shape = fastidious_type)) + 
    geom_point(cex = 2.5, stroke = 1) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL) +
    xlab("Threshold t") +
    ylab("Memory [MiB]") +
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(grid.draw(rbind(ggplotGrob(facet1), ggplotGrob(facet2), size = "last")), file = paste0(plot_file, "_combined", ".", dev), width = 16, height = 11, device = dev)
  
}


# get arguments & create visualisation
args = commandArgs(trailingOnly = T)

metrics_file <- args[1]
plot_file <- args[2]
min_t <- as.numeric(args[3])
max_t <- as.numeric(args[4])

facet_point_plots(metrics_file, paste0(plot_file, "_with_legend"), min_t:max_t, T)
facet_point_plots(metrics_file, paste0(plot_file, "_without_legend"), min_t:max_t, F)
