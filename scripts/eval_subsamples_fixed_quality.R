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

  submethod_type <- rep("-nf", nrow(data))
  submethod_type[grepl("-[es]-f1", data$method)] <- "-f1"
  submethod_type[grepl("-[es]-2f", data$method)] <- "-2f"
  
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
  submethod_type <- rep("-nf", nrow(data))

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
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[fastidious] <- "-2f"

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
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[grepl("-fast-length", data$method)] <- "-fast-length"
  submethod_type[grepl("-fast-abund", data$method)] <- "-fast-abund"
  submethod_type[grepl("-small-length", data$method)] <- "-small-length"
  submethod_type[grepl("-small-abund", data$method)] <- "-small-abund"

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
  submethod_type <- rep("-nf", nrow(data))
  submethod_type[grepl("-fast-length", data$method)] <- "-fast-length"
  submethod_type[grepl("-size-abund", data$method)] <- "-size-abund"
  submethod_type[grepl("-small-length", data$method)] <- "-small-length"
  submethod_type[grepl("-small-abund", data$method)] <- "-small-abund"

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
  submethod_type <- rep("-nf", nrow(data))
  
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
  submethod_type <- rep("-nf", nrow(data))
  
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
  submethod_type <- rep("-nf", nrow(data))
  
  data.frame(method, threshold, fastidious, fastidious_threshold, submethod_type, 
             data$rep, data$recall, data$precision, data$adjrandindex, stringsAsFactors = F)
  
}

# Loads and prepares data using above functions.
# Creates visualisations of clustering quality for different thresholds.
quality_facet_point_plots <- function(metrics_file, plot_file, with_legend = T, with_error_bars = F, 
                                      col_names = c("method", "threshold", "rep", "recall", "precision", 
                                                    "nmi", "randindex", "adjrandindex"),
                                      dev = "pdf") {

  # read & prepare data
  data <- read.csv(metrics_file, header = T, sep = ",", colClasses = c("character", rep("numeric", length(col_names) - 1)))
  
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

  method_levels = c("Swarm (v1.2.3)", "Swarm (v2)", "GeFaST (edit distance)", "GeFaST (scoring function)", 
                    "USEARCH", "VSEARCH", "CD-HIT", "DNACLUST", "Sumaclust")
  data$method[grepl("swarm-v1", data$method)] <- method_levels[1]
  data$method[grepl("swarm-v2", data$method)] <- method_levels[2]
  data$method[grepl("gefast-e", data$method)] <- method_levels[3]
  data$method[grepl("gefast-s", data$method)] <- method_levels[4]
  data$method[grepl("usearch", data$method)] <- method_levels[5]
  data$method[grepl("vsearch", data$method)] <- method_levels[6]
  data$method[grepl("cd-hit", data$method)] <- method_levels[7]
  data$method[grepl("dnaclust", data$method)] <- method_levels[8]
  data$method[grepl("sumaclust", data$method)] <- method_levels[9]
  data$method <- factor(data$method, levels = method_levels)

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
  averaged_data <- data.frame(group_by(data, method, threshold, fastidious, fastidious_threshold, submethod_type) 
                              %>% summarise(ave_recall = mean(recall), sd_recall = sd(recall), 
                                            ave_precision = mean(precision), sd_precision = sd(precision),
                                            ave_adjrandindex = mean(adjrandindex), sd_adjrandindex = sd(adjrandindex)))
  
  # adjust palette & threshold labels depending on set of compared methods
  pal <- c(brewer.pal(8, "Set3"), "darkgrey")
  averaged_data$threshold_label <- paste0("Threshold t = ", averaged_data$threshold, " (", formatC(((100 - averaged_data$threshold) / 100), digits = 2, format = "f"), ")")
  if (all(data$method %in% method_levels[c(2, 4)])) { # only Swarm (v2) and GeFaST (scoring function)
    
    pal <- pal[c(1, 4)]
    averaged_data$threshold_label <- paste0("Threshold t = ", averaged_data$threshold)
    
  } else if (all(data$method %in% method_levels[4:9])) { # all except Swarm (v1 & v2) and GeFaST (edit distance)
    pal <- pal[4:9]
  } else if (all(data$method %in% method_levels[c(2, 4, 5, 6, 7, 8, 9)])) { # all except Swarm (v1) and GeFaST (edit distance)
    pal <- pal[c(2, 4, 5, 6, 7, 8, 9)]
  } 

  thresholds <- unique(averaged_data$threshold)
  scale <- function(x) {sprintf("%.2f", x)}
  
  # plot recall measurements
  facets <- ggplot(averaged_data, aes(x = threshold, y = ave_recall, group = interaction(method, submethod_type), colour = method, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = paste0(thresholds, "\n", formatC(((100 - thresholds) / 100), digits = 2, format = "f"))) + 
    scale_y_continuous(labels = scale) +
    scale_colour_manual(values = pal) +
    xlab("Threshold t") +
    ylab("Recall") +
    labs(colour = "Method", shape = "Submethod") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_recall - sd_recall, ymax = ave_recall + sd_recall), width = 0.5, position = position_dodge(width = 0.6))
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.5, 0.5), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_recall", ".", dev), width = 16, height = 9, device = dev)
  
  # plot precision measurements
  facets <- ggplot(averaged_data, aes(x = threshold, y = ave_precision, group = interaction(method, submethod_type), colour = method, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = paste0(thresholds, "\n", formatC(((100 - thresholds) / 100), digits = 2, format = "f"))) + 
    scale_y_continuous(labels = scale) +
    scale_colour_manual(values = pal) +
    xlab("Threshold t") +
    ylab("Precision") +
    labs(colour = "Method", shape = "Submethod") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_precision - sd_precision, ymax = ave_precision + sd_precision), width = 0.5, position = position_dodge(width = 0.6))
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.5, 0.5), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_precision", ".", dev), width = 16, height = 9, device = dev)
  
  # plot adjusted Rand index measurements
  facets <- ggplot(averaged_data, aes(x = threshold, y = ave_adjrandindex, group = interaction(method, submethod_type), colour = method, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = paste0(thresholds, "\n", formatC(((100 - thresholds) / 100), digits = 2, format = "f"))) + 
    scale_y_continuous(labels = scale) +
    scale_colour_manual(values = pal) +
    xlab("Threshold t") +
    ylab("Adjusted Rand index") +
    labs(colour = "Method", shape = "Submethod") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets <- facets + geom_errorbar(aes(ymin = ave_adjrandindex - sd_adjrandindex, ymax = ave_adjrandindex + sd_adjrandindex), width = 0.5, position = position_dodge(width = 0.6))
  }
  
  if (with_legend) {
    facets <- facets + theme(legend.position = c(0.5, 0.3), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets <- facets + theme(legend.position = "none")
  }
  
  ggsave(facets, file = paste0(plot_file, "_adjrandindex", ".", dev), width = 16, height = 9, device = dev)
  
  # plot all metrics
  facets1 <- ggplot(averaged_data, aes(x = threshold, y = ave_recall, group = interaction(method, submethod_type), colour = method, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = NULL) +
    scale_y_continuous(labels = scale) +
    xlab("") +
    ylab("Recall") +
    labs(colour = "Method", shape = "Submethod") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold"),
      legend.position = "none"
    )
  
  if (with_error_bars) {
    facets1 <- facets1 + geom_errorbar(aes(ymin = ave_recall - sd_recall, ymax = ave_recall + sd_recall), width = 0.5, position = position_dodge(width = 0.6))
  }
  
  facets2 <- ggplot(averaged_data, aes(x = threshold, y = ave_precision, group = interaction(method, submethod_type), colour = method, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = NULL) +
    scale_y_continuous(labels = scale) +
    xlab("") +
    ylab("Precision") +
    labs(colour = "Method", shape = "Submethod") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold")
    )
  
  if (with_error_bars) {
    facets2 <- facets2 + geom_errorbar(aes(ymin = ave_precision - sd_precision, ymax = ave_precision + sd_precision), width = 0.5, position = position_dodge(width = 0.6))
  }
  
  if (with_legend) {
    facets2 <- facets2 + theme(legend.position = c(0.25, 0.3), legend.direction = "vertical", legend.box = "horizontal")
  } else {
    facets2 <- facets2 + theme(legend.position = "none")
  }
  
  facets3 <- ggplot(averaged_data, aes(x = threshold, y = ave_adjrandindex, group = interaction(method, submethod_type), colour = method, shape = submethod_type)) + 
    geom_point(cex = 2.5, stroke = 1, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(1, 0, 5, 4, 7, 3, 12, 10)) + 
    scale_x_continuous(breaks = thresholds, minor_breaks = NULL, labels = paste0(thresholds, "\n", formatC(((100 - thresholds) / 100), digits = 2, format = "f"))) + 
    scale_y_continuous(labels = scale) +
    scale_colour_manual(values = pal) +
    xlab("Threshold t") +
    ylab("Adjusted Rand index") +
    labs(colour = "Method", shape = "Submethod") +
    theme(
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 13, face = "bold"), 
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 13, face = "bold"),
      legend.position = "none"
    )
  
  if (with_error_bars) {
    facets3 <- facets3 + geom_errorbar(aes(ymin = ave_adjrandindex - sd_adjrandindex, ymax = ave_adjrandindex + sd_adjrandindex), width = 0.5, position = position_dodge(width = 0.6))
  }
  
  ggsave(grid.draw(rbind(ggplotGrob(facets1), ggplotGrob(facets2), ggplotGrob(facets3), size = "last")), file = paste0(plot_file, "_combined", ".", dev), width = 16, height = 18, device = dev)
  
}


# get arguments & create visualisation
args = commandArgs(trailingOnly = T)

metrics_file <- args[1]
plot_file <- args[2]

if (length(args) > 2) {

  cols = strsplit(args[3], ",")[[1]]
  quality_facet_point_plots(metrics_file, paste0(plot_file, "_with_legend"), T, T, cols)
  quality_facet_point_plots(metrics_file, paste0(plot_file, "_without_legend"), F, T, cols)

} else {

  quality_facet_point_plots(metrics_file, paste0(plot_file, "_with_legend"), T, T)
  quality_facet_point_plots(metrics_file, paste0(plot_file, "_without_legend"), F, T)

}
