#!/usr/bin/env Rscript
#
# PURPOSE: Reads raw coverage data from multiple mapping experiments, merges it
#          with sample metadata, and generates quality control plots.
#
# USAGE:
# Rscript 01_qc.R --samples_list acc_list.sim.no_ont \
#                 --metadata metadata2.csv \
#                 --coverage_dir ./ \
#                 --out_csv coverage_summary.csv \
#                 --out_prefix qc_plots

# Load required libraries
library(optparse)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# --- Command Line Argument Parsing ---
option_list <- list(
  make_option(c("-s", "--samples_list"), type="character", default=NULL, help="Path to a file with one sample name per line.", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="Path to the metadata CSV/TSV file.", metavar="character"),
  make_option(c("-c", "--coverage_dir"), type="character", default=".", help="Parent directory containing experiment subfolders (e.g., LSDV1, LSDVG).", metavar="character"),
  make_option(c("-o", "--out_csv"), type="character", default="coverage_summary.csv", help="Output path for the summary CSV file.", metavar="character"),
  make_option(c("-p", "--out_prefix"), type="character", default="qc_plots", help="Prefix for the output PDF plot files.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$samples_list) || is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Samples list and metadata file must be supplied.", call.=FALSE)
}

# --- Main Analysis ---

# 1. Read input files
samples <- readLines(opt$samples_list)
metadata <- read.csv(opt$metadata, sep="\t", header = FALSE, col.names = c("Sample", "Library_type"))

# Dynamically find experiment subfolders within the coverage directory
subfolders <- list.dirs(path = opt$coverage_dir, full.names = FALSE, recursive = FALSE)
# Filter for expected folder names if necessary, for now we use all directories found
print(paste("Found experiment folders:", paste(subfolders, collapse=", ")))


# 2. Loop through folders and samples to read coverage data
coverage_data <- list()
for (folder in subfolders) {
  for (sample in samples) {
    file_path <- file.path(opt$coverage_dir, folder, "COVERAGE", paste0(sample, ".coverage.txt"))
    if (file.exists(file_path)) {
      cov_data <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE) %>%
        select(rname = `#rname`, meandepth, meanbaseq, meanmapq) %>%
        rename(depth = meandepth, BQ = meanbaseq, MQ = meanmapq) %>%
        mutate(
          rname = ifelse(rname == "KX894507", "KX894508", rname),
          Sample = sample,
          Folder = folder
        ) %>%
        group_by(rname, Sample, Folder) %>%
        summarise(across(c(depth, BQ, MQ), sum, na.rm = TRUE), .groups = 'drop')
      
      coverage_data[[paste0(folder, "_", sample)]] <- cov_data
    }
  }
}

# 3. Combine and process data
comparison_data <- bind_rows(coverage_data) %>%
  merge(metadata, by = "Sample", all.x = TRUE)

# Create a lookup for folder names to make renaming cleaner
folder_map <- c(
  "LSDVG_3" = "Giraffe 3", "LSDVVG_3" = "VG-MAP 3",
  "LSDVG_6" = "Giraffe 6", "LSDVVG_6" = "VG-MAP 6",
  "LSDVG_ALL" = "Giraffe all", "LSDVVG_ALL" = "VG-MAP all",
  "LSDV1" = "Minimap2", "LSDVG" = "Giraffe 1", "LSDVVG" = "VG-MAP 1"
)
comparison_data$Folder <- folder_map[comparison_data$Folder]

comparison_data$depth <- as.numeric(comparison_data$depth)
comparison_data$BQ <- as.numeric(comparison_data$BQ)
comparison_data$MQ <- as.numeric(comparison_data$MQ)

# Write summary CSV
write.csv(comparison_data, opt$out_csv, row.names = FALSE)
print(paste("Summary data written to:", opt$out_csv))

# 4. Generate Plots
# Subset data for plotting to match original script's logic
comparison_data1 <- subset(comparison_data, Folder %in% c("Minimap2", "Giraffe 1", "VG-MAP 1", "Giraffe 3", "VG-MAP 3", "Giraffe 6", "VG-MAP 6"))

# Depth vs. BQ
pdf(paste0(opt$out_prefix, "_Depth_BQ.pdf"), width=10, height=4)
p_bq <- ggplot(comparison_data1, aes(x = log10(depth + 1.1), y = BQ, colour = Library_type)) +
  geom_point(alpha = 0.5, size = 1.1) +
  facet_wrap(~Folder, nrow = 1) +
  ylim(24, 101) +
  labs(x = "Log10-scaled read depth", y = "Base Quality (BQ)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8), strip.text = element_text(size = 8)) +
  scale_fill_manual(values = brewer.pal(7, "Set2")) +
  geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "dashed", alpha = 0.5) +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = 'left', label.y.npc = 'top')
print(p_bq)
dev.off()

# Depth vs. MQ
pdf(paste0(opt$out_prefix, "_Depth_MQ.pdf"), width=10, height=4)
p_mq <- ggplot(comparison_data1, aes(x = log10(depth + 1.1), y = MQ, colour = Library_type)) +
  geom_point(alpha = 0.5, size = 1.1) +
  facet_wrap(~Folder, nrow = 1) +
  labs(x = "Log10-scaled read depth", y = "Mapping Quality (MQ)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8), strip.text = element_text(size = 8)) +
  scale_fill_manual(values = brewer.pal(7, "Set2")) +
  geom_smooth(method = "lm", se = TRUE, color = "grey", linetype = "dashed", alpha = 0.5) +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = 'left', label.y.npc = 'top')
print(p_mq)
dev.off()

print("Script finished successfully. Plots are saved with the prefix:")
print(opt$out_prefix)
