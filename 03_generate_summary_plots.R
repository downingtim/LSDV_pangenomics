#!/usr/bin/env Rscript
#
# PURPOSE: Reads various statistics files (SNPS, GAMSTATS, BAMSTATS, TSTV)
#          and generates numerous summary boxplots comparing performance
#          across different variant callers and mapping strategies.
#          This script consolidates the logic from the original 03, 04, and 08.
#
# USAGE:
# Rscript 03_generate_summary_plots.R --input_dir <path/to/metrics> \
#                                     --output_dir <path/to/plots>

# --- Load Libraries ---
library(optparse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# --- Command Line Argument Parsing ---
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=".", help="Directory containing the metrics files (e.g., SNPS.ILLUMINA.txt).", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="summary_plots", help="Directory to save output PDF plots.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dir.create(opt$output_dir, showWarnings = FALSE)

# --- Helper Functions ---

# Function to rename sources for consistent plotting
rename_sources <- function(dataset) {
  source_map <- c(
    'LSDVVG' = 'VG-MAP_1', 'LSDVVG_3' = 'VG-MAP_3', 'LSDVVG_6' = 'VG-MAP_6', 'LSDVVG_ALL' = 'VG-MAP_All',
    'LSDVG' = 'Giraffe_1', 'LSDVG_3' = 'Giraffe_3', 'LSDVG_6' = 'Giraffe_6', 'LSDVG_ALL' = 'Giraffe_All',
    'LSDV1' = 'Minimap2_1', 'LSDV1_3' = 'Minimap2_3', 'LSDV1_6' = 'Minimap2_6', 'LSDV1_ALL' = 'Minimap2_All',
    'M' = 'Minimap2_1', '1' = 'VG-MAP_1', '3' = 'VG-MAP_3', '6' = 'VG-MAP_6', 'ALL' = 'VG-MAP_All',
    '1_GBWT' = 'Giraffe_1', '3_GBWT' = 'Giraffe_3', '6_GBWT' = 'Giraffe_6', 'ALL_GBWT' = 'Giraffe_All'
  )
  dataset$Source <- recode(dataset$Source, !!!source_map)
  return(dataset)
}

# General plotting function
generate_plots <- function(dataset, base_name, y_label, y_limit, log_scale = FALSE) {
  
  # Clean and prepare dataset
  dataset <- rename_sources(dataset)
  if (log_scale) {
    dataset <- dataset %>% filter(Rate > 0)
    dataset$PlotValue <- log10(dataset$Rate)
  } else {
    dataset$PlotValue <- dataset$Rate
  }
  
  # Plot 1: Boxplot by Caller
  p_caller <- ggplot(dataset, aes(x=Caller, y=PlotValue, color=Caller)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    geom_jitter(size=1, alpha=0.4) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=10, angle=90)) +
    labs(y=y_label) +
    ylim(y_limit)
  ggsave(file.path(opt$output_dir, paste0(base_name, "_by_caller.pdf")), p_caller, width=5, height=5)

  # Plot 2: Boxplot by Source
  p_source <- ggplot(dataset, aes(x=Source, y=PlotValue, color=Source)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA) +
    geom_jitter(size=1, alpha=0.4) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=10, angle=90)) +
    labs(y=y_label) +
    ylim(y_limit)
  ggsave(file.path(opt$output_dir, paste0(base_name, "_by_source.pdf")), p_source, width=8, height=5)
  
  # Plot 3: Faceted by Source
  p_facet_source <- p_caller + facet_wrap(~Source, ncol = 4)
  ggsave(file.path(opt$output_dir, paste0(base_name, "_caller_faceted_by_source.pdf")), p_facet_source, width=10, height=8)

  # Plot 4: Faceted by Caller
  p_facet_caller <- p_source + facet_wrap(~Caller, ncol = 3)
  ggsave(file.path(opt$output_dir, paste0(base_name, "_source_faceted_by_caller.pdf")), p_facet_caller, width=12, height=5)
}

# --- Main Logic ---

# 1. Process SNP data
print("Processing SNP data...")
try({
  snps_illumina <- read.csv(file.path(opt$input_dir, "SNPS.ILLUMINA.txt"), sep="\t", header=FALSE, col.names=c("Sample", "Source", "Caller", "Rate"))
  snps_ont <- read.csv(file.path(opt$input_dir, "SNPS.ONT.txt"), sep="\t", header=FALSE, col.names=c("Sample", "Caller", "Source", "Rate"))
  generate_plots(snps_illumina, "snps_illumina", "Log10(Number of SNPs)", c(1, 6), log_scale=TRUE)
  generate_plots(snps_ont, "snps_ont", "Log10(Number of SNPs)", c(1, 6), log_scale=TRUE)
}, silent=TRUE)

# 2. Process Ti/Tv data
print("Processing Ti/Tv data...")
try({
  tstv_files <- list.files(opt$input_dir, pattern="^TSTV.*\\.txt$", full.names=TRUE)
  tstv_data <- do.call(rbind, lapply(tstv_files, function(f) {
      df <- read.csv(f, sep="\t", header=FALSE, col.names=c("Sample", "Caller", "Source", "Rate"))
      df$Library <- gsub(".*\\.(.*?)\\.txt", "\\1", basename(f))
      return(df)
  }))
  generate_plots(tstv_data, "tstv_all_libs", "Ti/Tv Ratio", c(0, 8))
}, silent=TRUE)


# 3. Process GAM/BAM stats
print("Processing GAM and BAM stats...")
for (stat_type in c("GAMSTATS", "BAMSTATS")) {
  try({
    files <- list.files(opt$input_dir, pattern=paste0("^", stat_type, ".*\\.txt$"), full.names=TRUE)
    data <- do.call(rbind, lapply(files, function(f) {
      df <- read.csv(f, sep="\t", header=FALSE, col.names=c("Sample", "Source", "Rate"))
      df$Library <- gsub(".*\\.(.*?)\\.txt", "\\1", basename(f))
      return(df)
    }))
    data <- data %>% filter(Rate != "NA") %>% mutate(Rate = as.numeric(Rate))
    
    # Add a dummy "Caller" column for the function to work
    data$Caller <- "Aligner"

    y_label <- ifelse(stat_type == "GAMSTATS", "% Reads Mapped (GAM)", "% Reads Mapped (BAM)")
    base_name <- tolower(stat_type)
    
    generate_plots(data, base_name, y_label, c(0, 101))
  }, silent=TRUE)
}

print(paste("Script finished. All plots saved in:", opt$output_dir))
