library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(gridExtra)
library(grid)
library(cowplot)
library(vcfR)  # Added for reading VCF files

# Create plot theme for consistency
my_theme <- theme_minimal() +
  theme(    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")  )

file_paths <- c(  "Unique.1.minimap2.txt",
  "Unique.13.minimap2.txt",  "Unique.minimap2.1.txt",  "Unique.minimap2.13.txt" )

vcf_folder <- "FINAL_SNPS_2"
vcf_files <- list.files(vcf_folder, pattern = "\\.giraffe13\\.vcf\\.gz$", full.names = T)

# Function to read VCF file and extract positions
read_vcf_file <- function(file_path) {
  # Extract sample name from filename
  sample_name <- str_extract(basename(file_path), "^[^.]+")
  
  # Read VCF file
  vcf <- read.vcfR(file_path, verbose = FALSE)
  
  # Extract chromosome, position, and convert to data frame
  vcf_data <- data.frame( Chrom = getCHROM(vcf),    Pos = getPOS(vcf),    Sample = sample_name )
  return(vcf_data) }

# Read all VCF files and combine them
vcf_all_data <- map_df(vcf_files, read_vcf_file)
vcf_filtered <- vcf_all_data %>%
  filter(Chrom %in% c("OQ511520", "KX894508", "KX764645"))
vcf_counts <- vcf_filtered %>%  group_by(Chrom, Pos) %>%
  summarise(Count = n(), .groups = "drop")

# Create plot A - a faceted plot for the three chromosomes
chroms <- c("OQ511520", "KX894508", "KX764645")
my_chrom_colors <- c("purple", "darkgrey", "darkred")
names(my_chrom_colors) <- chroms

p_chrom <- ggplot(vcf_counts, aes(x = Pos, y = Count, color = Chrom)) +
  geom_point(alpha = 0.5, size = 1.3) +
  facet_wrap(~ Chrom, ncol = 1, scales = "free_x") +
  scale_color_manual(values = my_chrom_colors) +
  labs(x = "Position", y = "Frequency")+  my_theme +
  theme(    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),    legend.position = "none"  ) +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)/10000) * 10000, ceiling(max(x)/10000) * 10000, by = 10000),
    minor_breaks = function(x) seq(floor(min(x)/10000) * 10000 + 5000, ceiling(max(x)/10000) * 10000, by = 10000) ) +
    annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -1, vjust = 2, size = 6, fontface = "bold") 

regions <- list(  region1 = c(5200, 6500),  region2 = c(8055, 8231),
  region3 = c(136200, 140400) )
region_widths <- c(4, 2.5, 5)  # Proportional widths

# Define custom names for the file sources
file_names <- c(
  "Unique.1.minimap2.txt" = "Giraffe_1 not in Minimap2",
  "Unique.13.minimap2.txt" = "Giraffe_1&3 not in Minimap2",
  "Unique.minimap2.1.txt" = "Minimap2 not in Giraffe_1",
  "Unique.minimap2.13.txt" = "Minimap2 not in Giraffe_1&3" )

# Function to read and process each file
read_snp_file <- function(file_path) {
  data <- read_tsv(file_path,
                  col_names = c("Chrom", "Pos", "Count", "Allele1", "Allele2"),
                  col_types = cols(
                    Chrom = col_character(),
                    Pos = col_integer(),
                    Count = col_integer(),
                    Allele1 = col_character(),
                    Allele2 = col_character()                  ))
  # Add source file information
  data$source_file <- basename(file_path)
  return(data) }

# Read all files and combine them
all_data <- map_df(file_paths, read_snp_file)
all_data <- all_data %>%  mutate(source_name = file_names[source_file])

# Print summary of data
if (nrow(all_data) > 0) {
  cat("Data summary:\n")
  print(summary(all_data))
  cat("\nNumber of records per file:\n")
  print(table(all_data$source_file))
} else {
  cat("No data was loaded. Please check file paths.\n")
  quit() }

# Set a distinct color palette
my_colors <- c("cyan", "navy", "red", "orange")
names(my_colors) <- unique(all_data$source_name)

# Main overlaid plot with transparent points and smaller point size
p_main <- ggplot(all_data, aes(x = Pos, y = Count, color = source_name)) +
  geom_jitter(alpha = 0.5, size = 2, width=1) + ylim(0,31) +
  scale_color_manual(values = my_colors, name = "Data Source") +
  labs(x = "Position", y = "Frequency") +  my_theme +
scale_x_continuous(
    breaks = function(x) seq(floor(min(x)/10000) * 10000, ceiling(max(x)/10000) * 10000, by = 10000),
    minor_breaks = function(x) seq(floor(min(x)/10000) * 10000 + 5000, ceiling(max(x)/10000) * 10000, by = 10000)  )

# Create individual region plots
region_plots <- list()

for (i in seq_along(regions)) {
  region_name <- names(regions)[i]
  region_start <- regions[[i]][1]
  region_end <- regions[[i]][2]

  region_data <- all_data %>%  filter(Pos >= region_start & Pos <= region_end)


  if (nrow(region_data) > 0) {
    p_region <- ggplot(region_data, aes(x = Pos, y = Count, color = source_name)) +
      geom_jitter(alpha = 0.5, size = 2, width=1) +
      scale_color_manual(values = my_colors) + ylim(0,15)+
      labs(        x = "Position", y = "Count",
        title = paste(region_start, "-", region_end)      ) +
      my_theme +      theme(legend.position = "none",
            plot.title = element_text(size = 10, hjust = 0.5)) +
                  scale_x_continuous(limits = c(region_start, region_end))

    region_plots[[i]] <- p_region

    region_summary <- region_data %>%      group_by(source_name) %>%
      summarize(        count_n = n(),
        count_min = min(Count),        count_q1 = quantile(Count, 0.25),
        count_median = median(Count),        count_mean = mean(Count),
        count_q3 = quantile(Count, 0.75),        count_max = max(Count),
        count_sum = sum(Count)      )

    cat(paste0("\nNumerical summary for region ", region_start, "-", region_end, ":\n"))
    print(region_summary)

    # Save region summary
    summary_filename <- paste0("summary_region_", region_start, "-", region_end, ".csv")
    write_csv(region_summary, summary_filename)
  } else {
    cat(paste0("\nNo data found in region ", region_start, "-", region_end, "\n"))
    # Create empty plot with message if no data
    p_region <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = paste("No data in region", region_start, "-", region_end)) +
      theme_void() +
      theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

    region_plots[[i]] <- p_region
  }
}

# Create the complete layout with all plots (A, B, C, D, E)
# A is the chromosome panel plot
# B is the main overlay plot
# C, D, E are the region plots

for (i in seq_along(regions)) {
  region_name <- names(regions)[i]
  region_start <- regions[[i]][1]
  region_end <- regions[[i]][2]
  label <- LETTERS[i+2]  # C, D, E
  
  region_data <- all_data %>% filter(Pos >= region_start & Pos <= region_end)
  
  if (nrow(region_data) > 0) {
     
    p_region <- ggplot(region_data, aes(x = Pos, y = Count, color = source_name)) +
      geom_jitter(alpha = 0.5, size = 2, width=1) +
      scale_color_manual(values = my_colors) +  ylim(0,15) +
      labs(x = "Position", y = "Count",
           title = paste(region_start, "-", region_end)) +
      my_theme +      theme(legend.position = "none",
            plot.title = element_text(size = 10, hjust = 0.5)) +
      scale_x_continuous(limits = c(region_start, region_end)) +
      # Add label directly
      annotate("text", x = -Inf, y = Inf, label = label, 
               hjust = -1, vjust = 2, size = 6, fontface = "bold")
    
    region_plots[[i]] <- p_region
  } else {
    p_region <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label =
        paste("No data in region", region_start, "-", region_end)) +
      theme_void() +
      theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
      annotate("text", x = 0, y = 1, label = label, 
               hjust = -1, vjust = 2, size = 6, fontface = "bold")
    region_plots[[i]] <- p_region }
}      

# Extract the legend from the main plot
legend <- get_legend(p_main)

p_main_no_legend <- p_main + 
  theme(legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -1, vjust = 2, size = 6, fontface = "bold")
           
# Arrange the region plots in the bottom row with specified relative widths
bottom_row <- plot_grid(
  plotlist = region_plots,  nrow = 1,  rel_widths = region_widths )
             
# Arrange with chromosome plot (A) on top, then main plot (B), then region plots (C,D,E)
composite <- plot_grid(  p_chrom,  p_main_no_legend,   bottom_row,  ncol = 1,
  rel_heights = c(4, 2, 1)  ) # Adjust heights as needed

# Add the legend to the right
composite_with_legend <- plot_grid(  composite,  legend, rel_widths = c(8, 2)  )
ggsave("snp_variation_composite.pdf", composite_with_legend, width = 12, height = 10)
ggsave("chromosome_distribution.pdf", p_chrom, width = 10, height = 8)

# Get top positions with highest counts
top_positions <- all_data %>%  group_by(Pos) %>%
  summarize(    total_count = sum(Count),    avg_count = mean(Count),
    n_files = n_distinct(source_file),
    files_present = paste(sort(unique(source_name)), collapse = ", ")  ) %>%
  arrange(desc(total_count)) %>%  head(10)

cat("\nTop 10 positions with highest total counts:\n")
print(top_positions)
write_csv(all_data, "combined_snp_data.csv")
write_csv(vcf_counts, "vcf_position_counts.csv")