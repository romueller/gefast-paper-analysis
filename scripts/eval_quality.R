#!/usr/bin/env Rscript

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# Creates a faceted plot showing recall, precision and adjusted Rand index 
# of Swarm (v1.2.3 & v2) and GeFaST (edit distance & scoring function) 
# for the specified thresholds.
facet_plot_jitter <- function(metrics_file, plot_file, thresholds, with_legend = T, 
                              col_names = c("method", "threshold", "recall", "precision", 
                                            "nmi", "randindex", "adjrandindex"), 
                              dev = "pdf") {

  data <- read.table(metrics_file, header = T, dec = ".", sep = ",")
  colnames(data) <- col_names
  data <- melt(data[, c(1, 2, 3, 4, 7)], id = c("method", "threshold"))
  levels(data$variable) <- c("Recall", "Precision", "Adjusted Rand index")

  fastidious <- rep("deactivated", length(data$method))
  fastidious[grepl("f1", data$method)] <- "t + 1" 
  fastidious[grepl("2f", data$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  data$method <- gsub("-2f", "", gsub("-f1", "", data$method))
  data$method <- gsub("swarm-v1", method_levels[1], data$method)
  data$method <- gsub("swarm-v2", method_levels[2], data$method)
  data$method <- gsub("gefast-e", method_levels[3], data$method)
  data$method <- gsub("gefast-s", method_levels[4], data$method)
  data$method <- factor(data$method, levels = method_levels)

  data <- cbind(data, fastidious)
  remove(fastidious, method_levels)
  
  facets <- ggplot(data, aes(x = threshold, y = value)) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL) + 
    geom_jitter(aes(colour = method, shape = fastidious), cex = 2.5, stroke = 1, width = 0.25, height = 0) +
    xlab("Threshold t") +
    ylab("Metric values") + 
    facet_grid(. ~ variable) + 
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.5, 0.25), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets <- facets + theme(legend.position = "none")
  }

  
  ggsave(facets, file = paste0(plot_file, ".", dev), width = 16, height = 5, device = dev)
  
}

# Creates a faceted plot showing recall, precision and adjusted Rand index 
# of Swarm (v1.2.3 & v2) and GeFaST (edit distance & scoring function) 
# for the specified thresholds on two data sets.
facet_plot_jitter_both <- function(metrics_file1, metrics_file2, plot_file, thresholds, with_legend = T, 
                              col_names = c("method", "threshold", "recall", "precision", 
                                            "nmi", "randindex", "adjrandindex"), 
                              dev = "pdf") {
  
  # data set 1
  data <- read.table(metrics_file1, header = T, dec = ".", sep = ",")
  colnames(data) <- col_names
  data <- melt(data[, c(1, 2, 3, 4, 7)], id = c("method", "threshold"))
  levels(data$variable) <- c("Recall", "Precision", "Adjusted Rand index")
  
  fastidious <- rep("deactivated", length(data$method))
  fastidious[grepl("f1", data$method)] <- "t + 1" 
  fastidious[grepl("2f", data$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  data$method <- gsub("-2f", "", gsub("-f1", "", data$method))
  data$method <- gsub("swarm-v1", method_levels[1], data$method)
  data$method <- gsub("swarm-v2", method_levels[2], data$method)
  data$method <- gsub("gefast-e", method_levels[3], data$method)
  data$method <- gsub("gefast-s", method_levels[4], data$method)
  data$method <- factor(data$method, levels = method_levels)
  
  data <- cbind(data, fastidious)
  remove(fastidious, method_levels)
  
  facets1 <- ggplot(data, aes(x = threshold, y = value)) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = NULL) + 
    geom_jitter(aes(colour = method, shape = fastidious), cex = 2.5, stroke = 1, width = 0.25, height = 0) +
    xlab("") +
    ylab("Metric values") + 
    facet_grid(. ~ variable) + 
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold")
    )
  
  if (with_legend) {
    facets1 <- facets1 + theme(legend.position = c(0.5, 0.25), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets1 <- facets1 + theme(legend.position = "none")
  }
  
  
  # data set 2
  data <- read.table(metrics_file2, header = T, dec = ".", sep = ",")
  colnames(data) <- col_names
  data <- melt(data[, c(1, 2, 3, 4, 7)], id = c("method", "threshold"))
  levels(data$variable) <- c("Recall", "Precision", "Adjusted Rand index")
  
  fastidious <- rep("deactivated", length(data$method))
  fastidious[grepl("f1", data$method)] <- "t + 1" 
  fastidious[grepl("2f", data$method)] <- "2 * t"
  fastidious <- factor(fastidious, levels = c("deactivated", "t + 1", "2 * t"))
  
  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)")
  data$method <- gsub("-2f", "", gsub("-f1", "", data$method))
  data$method <- gsub("swarm-v1", method_levels[1], data$method)
  data$method <- gsub("swarm-v2", method_levels[2], data$method)
  data$method <- gsub("gefast-e", method_levels[3], data$method)
  data$method <- gsub("gefast-s", method_levels[4], data$method)
  data$method <- factor(data$method, levels = method_levels)
  
  data <- cbind(data, fastidious)
  remove(fastidious, method_levels)
  
  facets2 <- ggplot(data, aes(x = threshold, y = value)) + 
    scale_shape_manual(values = c(1, 0, 5)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL) + 
    geom_jitter(aes(colour = method, shape = fastidious), cex = 2.5, stroke = 1, width = 0.25, height = 0) +
    xlab("Threshold t") +
    ylab("Metric values") + 
    facet_grid(. ~ variable) + 
    labs(colour = "Method", shape = "Fastidious") + 
    theme(
      axis.title=element_text(size=20, face = "bold"),
      axis.text=element_text(size=15),
      legend.title=element_text(size=13, face = "bold"), 
      legend.text=element_text(size=13),
      strip.text.x=element_text(size=13, face = "bold"),
      legend.position = "none"
    )

  ggsave(grid.arrange(facets1, facets2, ncol = 1), file = paste0(plot_file, ".", dev), width = 16, height = 10, device = dev)
  
}


# get arguments & create visualisation
args = commandArgs(trailingOnly = T)

metrics_file <- args[1]
plot_file <- args[2]
min_t <- as.numeric(args[3])
max_t <- as.numeric(args[4])

facet_plot_jitter(metrics_file, paste0(plot_file, "_with_legend"), min_t:max_t, T)
facet_plot_jitter(metrics_file, paste0(plot_file, "_without_legend"), min_t:max_t, F)
