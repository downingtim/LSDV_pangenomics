#!/usr/bin/env Rscript

# =============================================================================
# GENETIC VARIANT DENSITY ANALYSIS AND VISUALIZATION
# =============================================================================
#
# DESCRIPTION:
# This script analyzes genetic variants from a VCF file and creates visualizations
# showing the distribution of mutations across a genome. It generates a density plot
# of mutations per genomic window alongside a visual representation of coding sequences
# (CDS) from a GenBank file.
#
# INPUT DATA:
# 1. VCF file: "vcf/gfavariants.vcf" - Contains genetic variants (SNPs, indels, etc.)
# 2. GenBank file: "KX894508.gb" - Contains genome annotation including CDS features
#
# OUTPUT DATA:
# 1. "bin_counts.csv" - CSV file with mutation counts per genomic window
# 2. "PVG_paper_Figure2.pdf" - Combined plot with mutation density and CDS annotations
#
# PARAMETERS:
# - genome_length: 151,000 bp (total genome length)
# - window_size: 400 bp (size of bins for density analysis)
#
# =============================================================================

library(VariantAnnotation)  # For reading and processing VCF files
library(ggplot2)           # For creating plots
library(dplyr)             # For data manipulation
library(scales)            # For formatting plot scales (comma formatting)
library(genbankr)          # For reading GenBank files
library(patchwork)         # For combining multiple plots

# =============================================================================
# READ AND PROCESS VCF FILE
# =============================================================================

# Read VCF file containing genetic variants
vcf <- readVcf("vcf/gfavariants.vcf")

# Filter for SNVs (Single Nucleotide Variants) only
vcf_snv <- vcf[isSNV(vcf)]

# Override filter to include all mutations (not just SNVs)
# Comment: This line overwrites the SNV filtering - decide if you want only SNVs or all variants
vcf_snv <- vcf # all mutations

# Extract genomic positions of variants
positions <- start(rowRanges(vcf_snv))

# =============================================================================
# BINNING ANALYSIS
# =============================================================================

# Define genome parameters
genome_length <- 151000  # Total genome length in base pairs
window_size <- 400       # Size of each bin in base pairs

# Create bins for the entire genome
# cut() divides positions into bins of specified size
bins <- cut(positions, breaks = seq(0, genome_length, by = window_size), right = FALSE)

# Count number of variants in each bin
bin_counts <- as.data.frame(table(bins))

# Parse bin information to get start positions
# Remove brackets and extract start coordinate from bin names
bin_counts$start <- as.numeric(gsub("\\[|\\)|\\]", "", sapply(strsplit(as.character(bin_counts$bins), ","), "[[", 1)))

# Calculate end position and midpoint for each bin
bin_counts$end <- bin_counts$start + window_size
bin_counts$midpoint <- bin_counts$start + window_size / 2

# Rename count column for clarity
colnames(bin_counts)[2] <- "count"

# Calculate summary statistics
median_val <- median(bin_counts$count)  # Median number of mutations per bin
top5_val <- quantile(bin_counts$count, 0.95)  # 95th percentile (top 5%)

# Export bin counts to CSV for further analysis
write.csv(bin_counts, "bin_counts.csv")

# =============================================================================
# CREATE MUTATION DENSITY PLOT WITH HIGHLIGHTED REGIONS
# =============================================================================

# Define regions of interest to highlight
regions <- data.frame(
  region = c("Region1", "Region2", "Region3"),
  start = c(5000, 7950, 136000),
  end = c(6700, 8360, 141000),
  color = c("#FFA500", "#FFA500", "#FFA500"),  # Orange, Sky Blue, Pale Green
  label_x = c(5500, 8230, 135500),  # Center positions for labels
  label_y = c(max(bin_counts$count) * 0.87,  # Position labels near top
              max(bin_counts$count) * 0.96,
              max(bin_counts$count) * 0.95) )

# Create main plot showing mutation density across genome
tree_plot <- ggplot(bin_counts, aes(x = midpoint, y = count)) +  
	  geom_rect(data = regions, 
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = color),
          alpha = 0.3, inherit.aes = FALSE) + scale_fill_identity() +
  geom_col(fill = "steelblue") +
  # Add horizontal reference lines for median and 95th percentile
  geom_hline(yintercept = median_val, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = top5_val, linetype = "dashed", color = "black", alpha = 0.5) +
  # Add region labels
  geom_text(data = regions,
            aes(x = label_x, y = label_y, label = region, color=color),
            color = "darkorange", fontface = "bold", size = 5,
            inherit.aes = FALSE) +
scale_x_continuous(name = "Genome Position (Kb)", 
                  breaks = seq(0, genome_length, by = 5000),
                  minor_breaks = seq(0, genome_length, by = 2500),
                  labels = seq(0, genome_length, by = 5000) / 1000) +
		  scale_color_identity() + 
  scale_y_continuous(name = "Mutations per Kb") +  
  # Apply minimal theme
  theme_minimal() +
  # Customize grid lines
  theme(panel.grid.minor = element_line(size = 0.3),
        panel.grid.major = element_line(size = 0.6),
        panel.grid.major.x = element_line(color = "gray70"),
        panel.grid.minor.x = element_line(color = "gray85"))

# Note: Individual density plot creation removed - only combined plot will be generated

# =============================================================================
# READ AND PROCESS GENBANK ANNOTATION
# =============================================================================

# Read GenBank file containing genome annotations
gb <- readGenBank("KX894508.gb")

# Extract CDS (Coding Sequence) features from GenBank file
cds_data <- gb@cds %>%  
  as.data.frame() %>%
  mutate(start = start,         # Start position of CDS
         end = end,             # End position of CDS
         strand = as.character(strand),  # Strand information (+ or -)
         # Set y-coordinate based on strand for visual separation
         y = ifelse(strand == "+", 1, -1))

# =============================================================================
# CREATE CDS ANNOTATION PLOT WITH HIGHLIGHTED GENES
# =============================================================================

# Define genes of interest to highlight (corresponding to the regions above)
highlighted_genes <- c("LD008", "LD009", "LD011", "LD012",
                   "LD144", "LD145", "LD146", "LD147")

# Assuming gene names are stored in a column (you may need to adjust this
# based on your actual GenBank file structure)
# Add a column to identify highlighted genes
cds_data <- cds_data %>%
  mutate(highlight = ifelse(gene %in% highlighted_genes, "orange", "grey50"))

# Create plot showing CDS features
cds_plot <- ggplot(cds_data, aes(xmin = start, xmax = end, ymin = 0, ymax = y)) +
  # Add orange transparent boxes for regions of interest (same as upper plot)
  geom_rect(data = regions, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "orange", alpha = 0.3, inherit.aes = FALSE) +
  # Draw rectangles for each CDS with conditional coloring
  geom_rect(aes(fill = highlight), color = "black", linewidth = 0.2) +  
  # Add manual color scale for highlighting specific genes
  scale_fill_identity() +
  # Add horizontal line at y=0 to separate strands
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  # Apply minimal theme with smaller base font size
  theme_minimal(base_size = 11) +
  # Set axis labels
  labs(x = "Genome Position (bp)", y = "CDS Strand") +
  # Match x-axis formatting with main plot
scale_x_continuous(name = "Genome Position (Kb)", 
                  breaks = seq(0, genome_length, by = 5000),
                  minor_breaks = seq(0, genome_length, by = 2500),
                  labels = seq(0, genome_length, by = 5000) / 1000) +
  # Set y-axis with custom labels for strands
  scale_y_continuous(breaks = c(-1, 1), labels = c("-", "+")) +
  # Customize theme to hide y-axis elements
  theme(axis.text.y = element_blank(),      # Hide y-axis text
        axis.ticks.y = element_blank(),     # Hide y-axis ticks
        axis.title.y = element_blank(),     # Hide y-axis title
        # Match grid formatting with main plot
        panel.grid.minor = element_line(size = 0.3),
        panel.grid.major = element_line(size = 0.6),
        panel.grid.major.x = element_line(color = "gray70"),
        panel.grid.minor.x = element_line(color = "gray85"))

# =============================================================================
# COMBINE PLOTS AND SAVE FINAL OUTPUT
# =============================================================================

# Combine mutation density plot with CDS annotation plot
# Heights argument controls relative size (4:1 ratio)
final_plot <- tree_plot / cds_plot + plot_layout(heights = c(4, 1))

# Save combined plot as PDF
ggsave("PVG_paper_Figure2.pdf", final_plot, width = 16, height = 5)

# Outputs
# 1. bin_counts.csv - Contains the raw data for mutation counts per bin
# 2. PVG_paper_Figure2.pdf - Combined mutation density and CDS annotation plot
#    with highlighted regions and genes