library(ggplot2)
library(dplyr)
library(ggrepel)
# install.packages("vcfR")
# install.packages("R.utils") 
library(vcfR) 
library(R.utils)

# Read the sample list from 'metadata2.csv' with tab-separated format
samples <- read.csv("metadata2.csv", sep = "\t", header = F, stringsAsFactors =F)
colnames(samples) <- c("Sample", "Library_type")
samples <- subset(samples, Library_type != "ONT")

# Define folders to process in order
folders <- c("SCREENED_MERGED_M/", "SCREENED_MERGED_1_GBWT/", "SCREENED_MERGED_3_GBWT/",
             "SCREENED_MERGED_6_GBWT/", "SCREENED_MERGED_ALL_GBWT/", "SCREENED_MERGED_1/",
             "SCREENED_MERGED_3/", "SCREENED_MERGED_6/", "SCREENED_MERGED_ALL/")

# Initialize an empty dataframe to store the results
results <- data.frame(Sample = character(),
                      Library_type = character(),
                      Folder = character(),
                      Total_Changes = integer(),
                      A_to_G = integer(),
                      A_to_C = integer(),
                      A_to_T = integer(),
                      G_to_C = integer(),
                      G_to_T = integer(),
                      T_to_C = integer(),
                      Ti = integer(), Tv = integer(),
                      Ti_Tv_Ratio = numeric(),
                      stringsAsFactors = FALSE)

# Define the function to parse mutation changes
parse_mutations <- function(file_path) {
  if (!file.exists(file_path) || file.info(file_path)$size == 0) {
    warning(paste("File", file_path, "not found or empty. Skipping."))
    return(NULL) }
  mutations <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(mutations) == 0) {
    warning(paste("File", file_path, "has no data. Skipping."))
    return(NULL) }
  
  colnames(mutations) <- c("SeqID", "Position", "Ref", "Alt", "Quality")
  
  A_to_G <- nrow(mutations %>% filter(Ref == "A" & Alt == "G"))
  A_to_C <- nrow(mutations %>% filter(Ref == "A" & Alt == "C"))
  A_to_T <- nrow(mutations %>% filter(Ref == "A" & Alt == "T"))
  G_to_C <- nrow(mutations %>% filter(Ref == "G" & Alt == "C"))
  G_to_T <- nrow(mutations %>% filter(Ref == "G" & Alt == "T"))
  T_to_C <- nrow(mutations %>% filter(Ref == "T" & Alt == "C"))
  
  total_changes <- A_to_G + A_to_C + A_to_T + G_to_C + G_to_T + T_to_C
  Ti <- A_to_G + T_to_C
  Tv <- A_to_C + A_to_T + G_to_C + G_to_T
  Ti_Tv_Ratio <- ifelse(Ti > 0 & Tv > 0, Ti / Tv, NA)
  
  return(list(Total_Changes = total_changes, A_to_G = A_to_G, A_to_C = A_to_C, A_to_T = A_to_T,
              G_to_C = G_to_C, G_to_T = G_to_T, T_to_C = T_to_C, Ti = Ti, Tv = Tv, Ti_Tv_Ratio = Ti_Tv_Ratio))
}

# Loop through each folder
for (folder in folders) {
  # Loop through each sample in the list and process its corresponding file
  for (i in 1:nrow(samples)) {
    sample_name <- samples$Sample[i]
    library_type <- samples$Library_type[i] 
    if (startsWith(sample_name, "paired_dat")) {  next  } 
    file_path <- file.path(folder, paste0(sample_name, ".sites")) 
    mutation_data <- parse_mutations(file_path) 
    if (is.null(mutation_data)) {next } 
    results <- rbind(results, data.frame(Sample = sample_name,
                                         Library_type = library_type,
                                         Folder = folder,
                                         Total_Changes = mutation_data$Total_Changes,
                                         A_to_G = mutation_data$A_to_G,
                                         A_to_C = mutation_data$A_to_C,
                                         A_to_T = mutation_data$A_to_T,
                                         G_to_C = mutation_data$G_to_C,
                                         G_to_T = mutation_data$G_to_T,
                                         T_to_C = mutation_data$T_to_C,
                                         Ti = mutation_data$Ti, Tv = mutation_data$Tv,
                                         Ti_Tv_Ratio = mutation_data$Ti_Tv_Ratio,
                                         stringsAsFactors = FALSE))
  }
}

# Create a summary dataframe that sums values across Library_type and Folder
summary_results <- results %>%
  group_by(Folder, Library_type) %>%
  summarise(Total_Changes = sum(Total_Changes),
            A_to_G = sum(A_to_G),
            A_to_C = sum(A_to_C),
            A_to_T = sum(A_to_T),
            G_to_C = sum(G_to_C),
            G_to_T = sum(G_to_T),
            T_to_C = sum(T_to_C),
            Ti = sum(Ti),
            Tv = sum(Tv)) %>%
  mutate(Ti_Tv_Ratio = ifelse(Ti > 0 & Tv > 0, Ti / Tv, NA))

results$Folder <- gsub("SCREENED_MERGED_", "", gsub("/", "", results$Folder))
summary_results$Folder <- gsub("SCREENED_MERGED_", "",
                               gsub("/", "", summary_results$Folder))

results$APOBEC = results$A_to_G/results$Total_Changes
# Print the results and summary dataframes
print(results)
print(summary_results)

write.csv(results, "ratios.results.csv")

results2 <-  subset(results,  Folder=="M" | Folder=="3_GBWT" |Folder=="1_GBWT" | Folder=="1_GBWT")
             
# Replace "M" with "Minimap2" in the Folder column
results2$Folder <- ifelse(results2$Folder == "M", "Minimap2", results2$Folder)

# Filter out rows where Ti_Tv_Ratio is NA and select non-NA data
filtered_results <- results2 %>% filter(!is.na(Ti_Tv_Ratio)) 
data1 = filtered_results %>% filter(Ti_Tv_Ratio > 6)

pdf("titv.plot.pdf", width=5, height=6)
ggplot(filtered_results, aes(x = Ti_Tv_Ratio, fill = Library_type)) +
  geom_histogram(binwidth = 0.4, position = "dodge") +
  labs(  x = "Ti/Tv Ratio", y = "Count") +
  theme_minimal() + scale_fill_brewer(palette = "Set2") +
  facet_wrap(~ Folder, ncol = 1) + geom_text_repel(data = data1, 
                  aes(x = Ti_Tv_Ratio, y = 0, label = Sample), 
                  size = 3, 
                  color = "red", 
                  box.padding = 1.0,  # Increase space around the label
                  point.padding = 1.0, # Increase distance between label and bar
                  max.overlaps = Inf,  # Allow all labels to be displayed
                  nudge_y = 5,         # Adjust vertical position of text
                  force = 10,            direction = "y")     # Ensure  
dev.off()

samples_all_greater_than_5 <- filtered_results %>%
  group_by(Sample) %>%
  summarize(all_greater_than_5 = all(Ti_Tv_Ratio > 5)) %>%
  filter(all_greater_than_5) %>%
  pull(Sample)  # Extract only the sample names
  
################## repeat on gzipped files for BCF and FB and VG

# Function to parse VCF files and compute Ti/Tv ratios
parse_vcf <- function(file_path, sample_name,library_type, category,folder1) {
  # Load the VCF file using vcfR
  vcf_data <- read.vcfR(file_path, verbose = F)
  
  # Extract the ALT and REF alleles from the VCF file
  ref <- vcf_data@fix[, "REF"]
  alt <- vcf_data@fix[, "ALT"]
  
  # Count transitions (Ti) and transversions (Tv)
  A_to_G <- sum(ref == "A" & alt == "G")
  G_to_A <- sum(ref == "G" & alt == "A")
  T_to_C <- sum(ref == "T" & alt == "C")
  C_to_T <- sum(ref == "C" & alt == "T")
  
  # Transversions (Tv)
  A_to_C <- sum(ref == "A" & alt == "C")
  A_to_T <- sum(ref == "A" & alt == "T")
  G_to_C <- sum(ref == "G" & alt == "C")
  G_to_T <- sum(ref == "G" & alt == "T")
  T_to_A <- sum(ref == "T" & alt == "A")
  T_to_G <- sum(ref == "T" & alt == "G")
  C_to_A <- sum(ref == "C" & alt == "A")
  C_to_G <- sum(ref == "C" & alt == "G")
  
  # Compute Ti and Tv
  Ti <- A_to_G + G_to_A + T_to_C + C_to_T
  Tv <- A_to_C + A_to_T + G_to_C + G_to_T + T_to_A + T_to_G + C_to_A + C_to_G
  
  # Calculate Ti/Tv ratio
  Ti_Tv_Ratio <- ifelse(Ti > 0 & Tv > 0, Ti / Tv, NA) 
  return(data.frame(Sample = sample_name, Library_type = library_type,
                    Category = category, Ti = Ti, Tv = Tv,
                    Ti_Tv_Ratio = Ti_Tv_Ratio, Folder=folder1))
}


# Read metadata file (adjust for your specific file structure)
samples <- read.csv("metadata2.csv", sep="\t", header=F, stringsAsFactors=F)
colnames(samples) <- c("Sample", "Library_type") 
samples <- samples %>% filter(!grepl("^paired_dat", Sample)) 

# List of folders to process
folders <- c("SCREENED_MERGED_M/", 
             "SCREENED_MERGED_1_GBWT/", 
             "SCREENED_MERGED_3_GBWT/", 
             "SCREENED_MERGED_6_GBWT/", 
             "SCREENED_MERGED_ALL_GBWT/", 
             "SCREENED_MERGED_1/", 
             "SCREENED_MERGED_3/", 
             "SCREENED_MERGED_6/", 
             "SCREENED_MERGED_ALL/")

results <- data.frame()

# Loop through each folder and process VCF files
for (folder in folders) {
  for (i in 1:nrow(samples)) {
    sample_name <- samples$Sample[i]
    library_type <- samples$Library_type[i] 
    bcf_file <- file.path(folder, paste0("BCF_", sample_name, ".vcf.gz"))
    if (file.exists(bcf_file)) {
      result_bcf <- parse_vcf(bcf_file, sample_name, library_type, "BCF", folder)
      results <- rbind(results,  result_bcf)  } 
    fb_file <- file.path(folder, paste0("FB_", sample_name, ".vcf.gz"))
    if (file.exists(fb_file)) {
      result_fb <- parse_vcf(fb_file, sample_name, library_type, "FB", folder)
      results <- rbind(results, result_fb) }
    vg_file <- file.path(folder, paste0("VG_", sample_name, ".vcf.gz"))
    if (file.exists(vg_file)) {
      result_vg <- parse_vcf(vg_file, sample_name, library_type, "VG", folder)
      results <- rbind(results, result_vg) }
  }
}

results$Folder <- gsub("SCREENED_MERGED_", "", gsub("/", "", results$Folder)) 
summary_df <- results %>%
  group_by(Folder, Library_type, Category) %>%
  summarize(Total_Ti = sum(Ti, na.rm = T), Total_Tv = sum(Tv, na.rm = T),
  Ti_Tv_Ratio=ifelse(sum(Tv,na.rm=T)>0,sum(Ti,na.rm=T)/sum(Tv,na.rm=T),NA)) 
print(summary_df)

# List of invalid SNP folders to process

# Function to parse VCF files and compute Ti/Tv ratios
parse_vcf2 <- function(file_path,sample_name,library_type, category,folder1) {
  # Load the VCF file using vcfR
  vcf_data <- read.csv(file_path, sep="\t", header=F, stringsAsFactors=F)
  
  # Extract the ALT and REF alleles from the VCF file
  ref <- vcf_data[,4]
  alt <- vcf_data[,5]
  
  # Count transitions (Ti) and transversions (Tv)
  A_to_G <- sum(ref == "A" & alt == "G")
  G_to_A <- sum(ref == "G" & alt == "A")
  T_to_C <- sum(ref == "T" & alt == "C")
  C_to_T <- sum(ref == "C" & alt == "T")
  
  # Transversions (Tv)
  A_to_C <- sum(ref == "A" & alt == "C")
  A_to_T <- sum(ref == "A" & alt == "T")
  G_to_C <- sum(ref == "G" & alt == "C")
  G_to_T <- sum(ref == "G" & alt == "T")
  T_to_A <- sum(ref == "T" & alt == "A")
  T_to_G <- sum(ref == "T" & alt == "G")
  C_to_A <- sum(ref == "C" & alt == "A")
  C_to_G <- sum(ref == "C" & alt == "G")
  
  # Compute Ti and Tv
  Ti <- A_to_G + G_to_A + T_to_C + C_to_T
  Tv <- A_to_C + A_to_T + G_to_C + G_to_T + T_to_A + T_to_G + C_to_A + C_to_G
  
  # Calculate Ti/Tv ratio
  Ti_Tv_Ratio <- ifelse(Ti > 0 & Tv > 0, Ti / Tv, NA) 
  return(data.frame(Sample = sample_name, Library_type = library_type,
                    Category = category, Ti = Ti, Tv = Tv,
                    Ti_Tv_Ratio = Ti_Tv_Ratio, Folder=folder1))
}

folders <- c("UNSCREENED_MERGED_M/", 
             "UNSCREENED_MERGED_1_GBWT/", 
             "UNSCREENED_MERGED_3_GBWT/", 
             "UNSCREENED_MERGED_6_GBWT/", 
             "UNSCREENED_MERGED_ALL_GBWT/", 
             "UNSCREENED_MERGED_1/", 
             "UNSCREENED_MERGED_3/", 
             "UNSCREENED_MERGED_6/", 
             "UNSCREENED_MERGED_ALL/") 
        
results_invalid <- data.frame()

print("invalid")
# Loop through each folder and process VCF files
for (folder in folders) {
    print(folder)
  for (i in 1:nrow(samples)) {
    sample_name <- samples$Sample[i]
    library_type <- samples$Library_type[i] 
    bcf_file <- file.path(folder, paste0( sample_name, ".BCF.sites"))
    if (file.exists(bcf_file)) {
      result_bcf <- parse_vcf2(bcf_file,sample_name,library_type,"BCF",folder)
      results_invalid <- rbind(results_invalid,  result_bcf)  } 
    fb_file <- file.path(folder, paste0(sample_name, ".FB.sites"))
    if (file.exists(fb_file)) {
      result_fb <- parse_vcf2(fb_file,sample_name, library_type,"FB",folder)
      results_invalid <- rbind(results_invalid, result_fb) }
    }
} 

results_invalid$Folder <- gsub("SCREENED_MERGED_", "",
   gsub("/", "", results_invalid$Folder)) 
summary_df_invalid <- results_invalid %>%
  group_by(Folder, Library_type, Category) %>%
  summarize(Total_Ti = sum(Ti, na.rm = T), Total_Tv = sum(Tv, na.rm = T),
  Ti_Tv_Ratio=ifelse(sum(Tv,na.rm=T)>0,sum(Ti,na.rm=T)/sum(Tv,na.rm=T),NA)) 
print(summary_df_invalid)
write.csv(summary_df_invalid, "summary_df_invalid.csv")
write.csv(summary_df, "summary_df_valid.csv")