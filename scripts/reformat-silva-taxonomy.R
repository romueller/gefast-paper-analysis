#!/usr/bin/env Rscript



# remove species information
removeSpecies <- function(vec) {
  posSpec <- regexpr(";D_6__", vec[2])
  if (posSpec == -1) {
    posSpec <- regexpr(";Ambiguous_taxa$", vec[2])
  }
  vec[2] <- substr(vec[2], 1, posSpec - 1)
  
  paste0(vec, collapse = "\t")
  
}



# get arguments & create visualisation
args = commandArgs(trailingOnly = T)

in_file <- args[1]
out_file <- args[2]

lines <- scan(file = in_file, what = "", sep = "\n", quiet = TRUE)
lines <- sapply(strsplit(lines, split="\t"), removeSpecies)
lines <- lines[!grepl("(Ambiguous_taxa|D_5__uncultured)", lines)]

write(lines, out_file)
