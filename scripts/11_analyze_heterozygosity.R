#!/usr/bin/env Rscript
#
# PURPOSE: Analyses and visualises allele frequencies from heterozygosity data.
#          It processes .vcf.het files from multiple experiments and generates
#          summary plots. This script consolidates the logic from the original
#          11_heterzygosity.R and 18_heterozygosity1.R files.
#
# USAGE:
# Rscript 11_analyze_heterozygosity.R --metadata metadata2.csv \
#                                     --input_dir ./ \
#                                     --output_dir ./het_plots

# --- Load Libraries ---
library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Command Line Argument Parsing ---
option_list <- list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="Path to the metadata CSV/TSV file.", metavar="character"),
  make_option(c("-i", "--input_dir"), type="character", default=".", help="Parent directory containing SCREENED_MERGED_* folders.", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="het_plots", help="Directory to save output PDF plots.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Metadata file must be supplied.", call.=FALSE)
}

# Create output directory if it doesn't exist
dir.create(opt$output_dir, showWarnings = FALSE)

# --- Helper Function ---
# Function to load heterozygosity data from a given folder
load_het_data <- function(folder_path, metadata_df) {
  all_sample_data <- list()
  
  for (i in 1:nrow(metadata_df)) {
    sample_name <- metadata_df$Sample[i]
    library_type <- metadata_df$Library_type[i]
    
    # Define paths for different callers
    file_bcf <- file.path(folder_path, paste0("BCF_", sample_name, ".vcf.het"))
    file_vg <- file.path(folder_path, paste0("VG_", sample_name, ".vcf.het"))
    file_fb <- file.path(folder_path, paste0("FB_", sample_name, ".vcf.het"))
    
    # Read data safely
    data_bcf <- if (file.exists(file_bcf)) try(read.table(file_bcf, header = FALSE), silent=TRUE) else NULL
    data_vg <- if (file.exists(file_vg)) try(read.table(file_vg, header = FALSE), silent=TRUE) else NULL
    data_fb <- if (file.exists(file_fb)) try(read.table(file_fb, header = FALSE), silent=TRUE) else NULL
    
    # Combine data from callers
    data <- bind_rows(
      if (!is.null(data_bcf) && !inherits(data_bcf, "try-error")) data_bcf,
      if (!is.null(data_vg) && !inherits(data_vg, "try-error")) data_vg,
      if (!is.null(data_fb) && !inherits(data_fb, "try-error")) data_fb
    )
    
    if (nrow(data) > 0) {
      data$Sample <- sample_name
      data$Library_type <- library_type
      all_sample_data[[i]] <- data
    }
  }
  
  combined_data <- bind_rows(all_sample_data)
  if (nrow(combined_data) > 0) {
    colnames(combined_data) <- c("Caller", "Ref", "Position", "Value", "Depth", "Sample", "Library_type")
  }
  return(combined_data)
}

# --- Main Analysis ---

# 1. Read metadata
metadata <- read.csv(opt$metadata, sep="\t", header=FALSE, col.names=c("Sample", "Library_type"))

# 2. Find all SCREENED_MERGED folders
exp_folders <- list.dirs(path = opt$input_dir, full.names = TRUE, recursive = FALSE)
screened_folders <- exp_folders[grep("SCREENED_MERGED_", basename(exp_folders))]

print(paste("Found", length(screened_folders), "SCREENED_MERGED folders to process."))

# 3. Loop through each folder, load data, and generate plots
for (folder in screened_folders) {
  folder_name <- basename(folder)
  print(paste("Processing folder:", folder_name))
  
  # Load data for the current folder
  het_data <- load_het_data(folder, metadata)
  
  if (nrow(het_data) == 0) {
    print(paste("No data found in", folder_name, ". Skipping."))
    next
  }
  
  # Filter data for plotting (as in original script)
  unwanted_prefixes <- c("paired", "449", "SRR233", "SRR255")
  het_data_filtered <- het_data %>%
    filter(Library_type %in% c("WGS", "AMPLICON", "Metagenomic")) %>%
    filter(!grepl(paste(unwanted_prefixes, collapse="|"), Sample)) %>%
    filter(Value > 0.05 & Value < 0.95) %>%
    distinct()
    
  if (nrow(het_data_filtered) == 0) {
    print(paste("No data remaining after filtering in", folder_name, ". Skipping."))
    next
  }

  # Plot 1: Allele Frequency Histogram
  plot_hist_path <- file.path(opt$output_dir, paste0("freq_histogram_", gsub("SCREENED_MERGED_", "", folder_name), ".pdf"))
  
  p_hist <- ggplot(het_data_filtered, aes(x = Value)) +
    geom_histogram(data = ~subset(., Caller=="FB"), alpha = 0.4, binwidth = 0.02, position = "identity", fill="blue") +
    geom_histogram(data = ~subset(., Caller=="BCF"), alpha = 0.4, binwidth = 0.02, position = "identity", fill="red") +
    geom_histogram(data = ~subset(., Caller=="VG"), alpha = 0.4, binwidth = 0.02, position = "identity", fill="green") +
    geom_density() +
    facet_wrap(~ Sample, scales = "free") +
    theme_bw() +
    labs(x = "Allele Frequency", y = "Count") +
    xlim(0, 1)

  ggsave(plot_hist_path, p_hist, width = 18, height = 10)
  
  # Plot 2: Allele Frequency vs. Position
  plot_scatter_path <- file.path(opt$output_dir, paste0("freq_vs_position_", gsub("SCREENED_MERGED_", "", folder_name), ".pdf"))

  p_scatter <- ggplot(het_data_filtered, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.3, size = 0.3) +
    facet_wrap(~Sample, ncol=4) +
    theme_minimal() +
    labs(x = "Genomic Position", y = "Allele Frequency", color = "Caller") +
    xlim(0, 150000) +
    ylim(0, 1)

  ggsave(plot_scatter_path, p_scatter, width = 20, height = 15)
}

print(paste("Script finished. All plots saved in:", opt$output_dir))
