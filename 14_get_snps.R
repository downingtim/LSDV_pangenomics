library(ggrepel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse) 

# Function to process each text file
process_file <- function(file_path) {
  lines <- readLines(file_path)
  
  # Remove lines matching the specified patterns
  lines <- lines[!grepl("Note: -w option not given, printing list of sites..", lines)] 
  lines <- lines[!grepl("uses Picard LiftoverVcf potential SNPs in", lines)]
  lines <- lines[!grepl("get valid SNP sites for all genome reference coordinates", lines)]
  
  # Extract the file stem
  file_stem <- sub("^snps\\.lsdv", "", sub("\\.txt$", "", basename(file_path)))
  
  # Initialize an empty data frame to store the results
  result <- data.frame(Sample = character(), Method = character(), SNPs = integer(), Indels = integer(),
  Invalid_SNPs = integer(), Invalid_Indels = integer(), File = character(), stringsAsFactors =F)
  
  # Variables to store the current sample and method data
  current_sample <- NULL
  
  # Iterate through the lines to extract data
  for (i in 1:length(lines)) {
   if ((grepl("^[^:,\\-]+$", lines[i])) &&(!( grepl("VCF_FILE", lines[i])))){
      current_sample <- trimws(lines[i])
    } else if (grepl("SNPs", lines[i]) && !(grepl("VG -", lines[i])) ) {  # Lines SNP data
      method_data <- strsplit(lines[i], ",")[[1]]
      method <- strsplit(method_data[1], "-")[[1]][1]
      
      snps <- as.integer(strsplit(method_data[1], "\\s+")[[1]][4])
      indels <- as.integer(strsplit(method_data[2], "\\s+")[[1]][3])
      invalid_snps <- as.integer(strsplit(method_data[3], "\\s+")[[1]][4])
      invalid_indels <- as.integer(strsplit(method_data[4], "\\s+")[[1]][4])
      
      result <- rbind(result,
          data.frame(Sample = current_sample, Method = method, 
               SNPs = snps, Indels = indels, 
               Invalid_SNPs = invalid_snps, Invalid_Indels = invalid_indels, 
               File = file_stem))   }  }
  
  # Add NA for VG data if the file is snps.lsdv_minimap2.txt
  # because VG has no SNPs there
  if ((file_stem == "_minimap2")|(file_stem == "_minimap3")|
      (file_stem == "_minimap6")|(file_stem == "_minimapAll")) {
    unique_samples <- unique(result$Sample)
  }
  return(result)
}

# File paths
file_paths <- list.files(pattern = "^snps\\.lsdv.*\\.txt$", full.names = T)

# Process all files and combine the results
all_data <- do.call(rbind, lapply(file_paths, process_file))

# Read the metadata file (tab-delimited, no header)
metadata <- read.csv("metadata2.csv", sep = "\t", header =F,
                        col.names = c("sample_name", "group_id"))

# Combine the data with the metadata
final_data <- merge(all_data, metadata, by.x = "Sample", by.y = "sample_name")

# Save the final long table to a CSV file
write.csv(final_data, "final_data.csv", row.names = F)

#################new part to find valid SNPs that did not transfer
 
count_entries_per_sample <- function(directory) {
  # Initialize an empty data frame to store the counts
  count_data <- data.frame(Sample = character(),  
  UNSCREENED_MERGED_M = integer(),   UNSCREENED_MERGED_3 = integer(),
                   UNSCREENED_MERGED_6 = integer(),  UNSCREENED_MERGED_ALL = integer(),
                   UNSCREENED_MERGED_3_GBWT = integer(),
                   UNSCREENED_MERGED_6_GBWT = integer(), UNSCREENED_MERGED_ALL_GBWT = integer(),
                           stringsAsFactors = FALSE)
  
  subfolders <- c(    "UNSCREENED_MERGED_3",
                  "UNSCREENED_MERGED_6", "UNSCREENED_MERGED_ALL",
                    "UNSCREENED_MERGED_3_GBWT",
                  "UNSCREENED_MERGED_6_GBWT", "UNSCREENED_MERGED_ALL_GBWT")
  
  for (subfolder in subfolders) {
    file_list <- list.files(path = subfolder, pattern = "\\.sites$", full.names = TRUE)
        
      file_list <- file_list[!grepl("\\.(FB|BCF|any|VG|bf)", file_list)]
      file_list <- file_list[grepl("\\.(reject)", file_list)]  
    
    for (file_path in file_list) {
      # Read the lines of the file
      lines <- readLines(file_path)
      sample_name <- sub(".*/([^/]+)\\.sites$", "\\1", file_path) 
      entry_count <- length(lines)
      
      if (sample_name %in% count_data$Sample) {
        count_data[count_data$Sample == sample_name, subfolder] <- entry_count
      } else {
        new_row <- data.frame(Sample = sample_name,   UNSCREENED_MERGED_3 = 0, 
                              UNSCREENED_MERGED_6 = 0,    UNSCREENED_MERGED_ALL = 0,
                                 UNSCREENED_MERGED_3_GBWT = 0,
                              UNSCREENED_MERGED_6_GBWT = 0,    UNSCREENED_MERGED_ALL_GBWT = 0,
                              stringsAsFactors = FALSE)
        new_row[[subfolder]] <- entry_count
        count_data <- rbind(count_data, new_row)       }    }
  }
  return(count_data)
}

# Missed mutations
entry_count_data <- count_entries_per_sample(".")  # contains Mongolia
entry_count_data$Sample  <- gsub(".reject", "",  entry_count_data$Sample)   # false mutations
colnames(entry_count_data) <- c("Sample",   "Missed_VG-MAP_3", "Missed_VG-MAP_6",
                            "Missed_VG-MAP_ALL",   "Missed_Giraffe_3",
                            "Missed_Giraffe_6", "Missed_Giraffe_ALL")                    
                            
############# VALID mutations for any site first

# For valid detected SNPs
count_entries_per_sample2 <- function(directory) {
  # Initialize an empty data frame to store the counts
  count_data <- data.frame(Sample = character(), 
                            SCREENED_MERGED_1 = integer(),   SCREENED_MERGED_3 = integer(),
                           SCREENED_MERGED_6 = integer(),  SCREENED_MERGED_ALL = integer(),
                           SCREENED_MERGED_1_GBWT = integer(),   SCREENED_MERGED_3_GBWT = integer(),
                           SCREENED_MERGED_6_GBWT = integer(),  SCREENED_MERGED_ALL_GBWT = integer(),
                           stringsAsFactors = FALSE)
  
  subfolders <- c("SCREENED_MERGED_1", "SCREENED_MERGED_3",
                  "SCREENED_MERGED_6", "SCREENED_MERGED_ALL",
                  "SCREENED_MERGED_1_GBWT", "SCREENED_MERGED_3_GBWT",
                  "SCREENED_MERGED_6_GBWT", "SCREENED_MERGED_ALL_GBWT")
  
  for (subfolder in subfolders) {
    file_list <- list.files(path = subfolder, pattern = "\\.sites$", full.names = TRUE)
    file_list <- file_list[!grepl("\\.(FB|BCF|any|VG|bf)", file_list)]  
    
    for (file_path in file_list) {
      # Read the lines of the file
      lines <- readLines(file_path)
      sample_name <- sub(".*/([^/]+)\\.sites$", "\\1", file_path)
      
      # Count the number of entries in the file (each line is an entry)
      entry_count <- length(lines)
      
      # Check if the sample is already in the data frame
      if (sample_name %in% count_data$Sample) {
        count_data[count_data$Sample == sample_name, subfolder] <- entry_count
      } else {
        # Add a new row for the sample
        new_row <- data.frame(Sample = sample_name, 
                                 SCREENED_MERGED_1 = 0,   SCREENED_MERGED_3 = 0,
                              SCREENED_MERGED_6 = 0,   SCREENED_MERGED_ALL = 0,
                              SCREENED_MERGED_1_GBWT = 0,  SCREENED_MERGED_3_GBWT = 0,
                              SCREENED_MERGED_6_GBWT = 0,  SCREENED_MERGED_ALL_GBWT = 0,
                              stringsAsFactors = FALSE)
        new_row[[subfolder]] <- entry_count
        count_data <- rbind(count_data, new_row)
      }
    }
  }
  return(count_data)
}

# Valid detected SNPs
entry_count_data2 <- count_entries_per_sample2(".")
colnames(entry_count_data2) <- c("Sample", "Valid_VG-MAP_1", "Valid_VG-MAP_3", "Valid_VG-MAP_6",
                            "Valid_VG-MAP_ALL", "Valid_Giraffe_1", "Valid_Giraffe_3",
                            "Valid_Giraffe_6", "Valid_Giraffe_ALL")

# # # # # # # # # # # # # # # # # # #  For BF only
                            
# Merge the two data frames by the 'Sample' column
merged_data <- merge(entry_count_data2, entry_count_data, by = "Sample", all = T) # valid + false mutations 

# Read the metadata file (tab-delimited, no header)
metadata <- read.csv("metadata2.csv", sep = "\t", header = FALSE, col.names = c("sample_name", "group_id"))

# Combine the data with the metadata
merged_data2 <- merge(merged_data, metadata, by.x = "Sample", by.y = "sample_name") 

# ##################### Get lifted SNP numbers  #####################

# Function to count the number of entries per sample per subfolder
count_entries_per_sample_lifted <- function(directory) {
  # Initialize an empty data frame to store the counts
  count_data <- data.frame(Sample = character(), 
                           LIFTOVER_3 = integer(),
                           LIFTOVER_6 = integer(),
                           LIFTOVER_ALL = integer(),
                           LIFTOVER_3_GBWT = integer(),
                           LIFTOVER_6_GBWT = integer(),
                           LIFTOVER_ALL_GBWT = integer(),
                           stringsAsFactors = FALSE)
  
  # List of subfolders
  subfolders <- c( "LIFTOVER_3", "LIFTOVER_6", "LIFTOVER_ALL", "LIFTOVER_3_GBWT",
                  "LIFTOVER_6_GBWT", "LIFTOVER_ALL_GBWT")
  
  for (subfolder in subfolders) {
    file_list <- list.files(path = subfolder, pattern = "\\.sites$", full.names = TRUE)
    
    # Exclude files containing .FB, .BCF, .any, or .VG
    file_list <- file_list[!grepl("\\.(FB|BCF|any|VG|bf)", file_list)]
    
    for (file_path in file_list) {
      # Read the lines of the file
      lines <- readLines(file_path)
      
      # Extract the sample name from the file path
      sample_name <- sub(".*/([^/]+)\\.sites$", "\\1", file_path)
      
      # Count the number of entries in the file (each line is an entry)
      entry_count <- length(lines)
      
      # Check if the sample is already in the data frame
      if (sample_name %in% count_data$Sample) {
        count_data[count_data$Sample == sample_name, subfolder] <- entry_count
      } else {
        # Add a new row for the sample
        new_row <- data.frame(Sample = sample_name, 
                              LIFTOVER_3 = 0,
                              LIFTOVER_6 = 0,
                              LIFTOVER_ALL = 0,
                              LIFTOVER_3_GBWT = 0,
                              LIFTOVER_6_GBWT = 0,
                              LIFTOVER_ALL_GBWT = 0,
                              stringsAsFactors = FALSE)
        new_row[[subfolder]] <- entry_count
        count_data <- rbind(count_data, new_row)   }
    }
  }
  return(count_data)
}

entry_count_data_lifted <- count_entries_per_sample_lifted(".")
# Remove Mongolia duplicates and A4 and VDM and B1
entry_count_data_lifted <- entry_count_data_lifted[-c(26:27),]
entry_count_data_lifted <- entry_count_data_lifted[-c(7:8,133),] 

# Merge the two data frames by the 'Sample' column
merged_data33 <- merge(merged_data2, entry_count_data_lifted, by = "Sample", all = T)

remove1 <- function(col1){ col2  <- gsub("LIFTOVER", "Lifted",  col1)
   return(col2)  }
colnames(merged_data33) <- remove1(colnames(merged_data33) ) 
 
############### BCF only ###############

count_entries_BFonly_function <- function(directory) {
  # Initialize an empty data frame to store the counts
  count_data <- data.frame(Sample = character(), 
                       SCREENED_MERGED_M = integer(),  SCREENED_MERGED_M3 = integer(),
                    SCREENED_MERGED_M6 = integer(),  SCREENED_MERGED_MALL = integer(),
                    SCREENED_MERGED_1 = integer(),   SCREENED_MERGED_3 = integer(),
                    SCREENED_MERGED_6 = integer(),  SCREENED_MERGED_ALL = integer(),
                   SCREENED_MERGED_1_GBWT =integer(), SCREENED_MERGED_3_GBWT=integer(),
               SCREENED_MERGED_6_GBWT =integer(), SCREENED_MERGED_ALL_GBWT=integer(),
                           stringsAsFactors = FALSE)
  
  subfolders <- c("SCREENED_MERGED_M", "SCREENED_MERGED_M3",
                     "SCREENED_MERGED_M6", "SCREENED_MERGED_MALL",
                  "SCREENED_MERGED_1", "SCREENED_MERGED_3",
                  "SCREENED_MERGED_6", "SCREENED_MERGED_ALL",
                  "SCREENED_MERGED_1_GBWT", "SCREENED_MERGED_3_GBWT",
                  "SCREENED_MERGED_6_GBWT", "SCREENED_MERGED_ALL_GBWT")
  
  for (subfolder in subfolders) {
    file_list <- list.files(path = subfolder, pattern = "\\.vcf.gz$", full.names = T)
    file_list <- file_list[grepl("BCF", file_list)] 
    
    for (file_path in file_list) { 
      lines <- readLines(gzfile(file_path), warn = F)  # gz instead
      sample_name <- sub(".*/([^/]+)\\.vcf.gz$", "\\1", file_path)
      sample_name <- sub("^BCF_", "", sample_name)
      entry_count <- length(lines)
      
      if (sample_name %in% count_data$Sample) {
        count_data[count_data$Sample == sample_name, subfolder] <- entry_count
      } else {
        new_row <- data.frame(Sample = sample_name, 
                              SCREENED_MERGED_M = 0, SCREENED_MERGED_M3 = 0,
                              SCREENED_MERGED_M6 = 0, SCREENED_MERGED_MALL = 0,
                              SCREENED_MERGED_1 = 0,   SCREENED_MERGED_3 = 0,
                              SCREENED_MERGED_6 = 0,   SCREENED_MERGED_ALL = 0,
                              SCREENED_MERGED_1_GBWT = 0,  SCREENED_MERGED_3_GBWT = 0,
                              SCREENED_MERGED_6_GBWT = 0,  SCREENED_MERGED_ALL_GBWT = 0,
                              stringsAsFactors = F)
        new_row[[subfolder]] <- entry_count
        count_data <- rbind(count_data, new_row)
      }
    }
  }
  return(count_data)
}

# VCFtools mutations only
entry_count_BFonly <- count_entries_BFonly_function(".") 
colnames(entry_count_BFonly)<-c("Sample", "BCF_M_1", "BCF_M_3", "BCF_M_6", "BCF_M_ALL",
                            "BCF_VG-MAP_1",   "BCF_VG-MAP_3", "BCF_VG-MAP_6",
                            "BCF_VG-MAP_ALL", "BCF_Giraffe_1", "BCF_Giraffe_3",
                            "BCF_Giraffe_6",  "BCF_Giraffe_ALL")  
                            
# Merge the two data frames by the 'Sample' column
merged_data5 <- merge(entry_count_BFonly, merged_data33, by = "Sample", all = T)  

############### FB only ###############

count_entries_FBonly_function <- function(directory) {
  # Initialize an empty data frame to store the counts
  count_data <- data.frame(Sample = character(), 
                           SCREENED_MERGED_M = integer(),  SCREENED_MERGED_M3 = integer(),
                           SCREENED_MERGED_M6 = integer(),  SCREENED_MERGED_MALL = integer(),         
                           SCREENED_MERGED_1 = integer(),   SCREENED_MERGED_3 = integer(),
                           SCREENED_MERGED_6 = integer(),  SCREENED_MERGED_ALL = integer(),
                           SCREENED_MERGED_1_GBWT =integer(), SCREENED_MERGED_3_GBWT=integer(),
                           SCREENED_MERGED_6_GBWT =integer(), SCREENED_MERGED_ALL_GBWT=integer(),
                           stringsAsFactors = FALSE)
  
  subfolders <- c("SCREENED_MERGED_M", "SCREENED_MERGED_M3",
                   "SCREENED_MERGED_M6", "SCREENED_MERGED_MALL",
                  "SCREENED_MERGED_1", "SCREENED_MERGED_3",
                  "SCREENED_MERGED_6", "SCREENED_MERGED_ALL",
                  "SCREENED_MERGED_1_GBWT", "SCREENED_MERGED_3_GBWT",
                  "SCREENED_MERGED_6_GBWT", "SCREENED_MERGED_ALL_GBWT")
  
  for (subfolder in subfolders) {
    file_list <- list.files(path = subfolder, pattern = "\\.vcf.gz$", full.names = T)
    file_list <- file_list[grepl("FB", file_list)] 
    
    for (file_path in file_list) { 
      lines <- readLines(gzfile(file_path), warn = F)  # gz instead
      sample_name <- sub(".*/([^/]+)\\.vcf.gz$", "\\1", file_path)
      sample_name <- sub("^FB_", "", sample_name)
      entry_count <- length(lines)
      
      if (sample_name %in% count_data$Sample) {
        count_data[count_data$Sample == sample_name, subfolder] <- entry_count
      } else {
        new_row <- data.frame(Sample = sample_name, 
                              SCREENED_MERGED_M = 0, SCREENED_MERGED_M3 = 0,
                              SCREENED_MERGED_M6 = 0, SCREENED_MERGED_MALL = 0,
                              SCREENED_MERGED_1 = 0,   SCREENED_MERGED_3 = 0,
                              SCREENED_MERGED_6 = 0,   SCREENED_MERGED_ALL = 0,
                              SCREENED_MERGED_1_GBWT = 0,  SCREENED_MERGED_3_GBWT = 0,
                              SCREENED_MERGED_6_GBWT = 0,  SCREENED_MERGED_ALL_GBWT = 0,
                              stringsAsFactors = F)
        new_row[[subfolder]] <- entry_count
        count_data <- rbind(count_data, new_row)
      }
    }
  }
  return(count_data)
}

# VCFtools mutations only
entry_count_FBonly <- count_entries_FBonly_function(".")
entry_count_FBonly$Sample  <- gsub(".vcf.gz", "",  entry_count_FBonly$Sample)  
entry_count_FBonly$Sample  <- gsub("FB_", "",  entry_count_FBonly$Sample) 
colnames(entry_count_FBonly) <- c("Sample", "FB_M_1", "FB_M_3",  "FB_M_6", "FB_M_ALL",
                            "FB_VG-MAP_1",   "FB_VG-MAP_3", "FB_VG-MAP_6",
                            "FB_VG-MAP_ALL", "FB_Giraffe_1", "FB_Giraffe_3",
                            "FB_Giraffe_6",  "FB_Giraffe_ALL")  
                            
# Merge the two data frames by the 'Sample' column
merged_data6 <- merge(merged_data5, entry_count_FBonly, by = "Sample", all = T)  
 merged_data6 <- merged_data6[, c(1, 28, 2:27, 29:46)]

############### VG only ###############

#count_entries_VGonly_function <- function(directory) {
  # Initialize an empty data frame to store the counts
 # count_data <- data.frame(Sample = character(), 
 #                            SCREENED_MERGED_1 = integer(),   SCREENED_MERGED_3 = integer(),
 #                          SCREENED_MERGED_6 = integer(),  SCREENED_MERGED_ALL = integer(),
 #                          SCREENED_MERGED_1_GBWT =integer(), SCREENED_MERGED_3_GBWT=integer(),
 #                          SCREENED_MERGED_6_GBWT =integer(), SCREENED_MERGED_ALL_GBWT=integer(),
#                           stringsAsFactors = FALSE)
  
 # subfolders <- c( "SCREENED_MERGED_1", "SCREENED_MERGED_3",
 #                 "SCREENED_MERGED_6", "SCREENED_MERGED_ALL",
 #                 "SCREENED_MERGED_1_GBWT", "SCREENED_MERGED_3_GBWT",
 #                 "SCREENED_MERGED_6_GBWT", "SCREENED_MERGED_ALL_GBWT")
  
 # for (subfolder in subfolders) {
 #   file_list <- list.files(path = subfolder, pattern = "\\.vcf.gz$", full.names = T)
 #   file_list <- file_list[grepl("VG", file_list)] 
    
 #   for (file_path in file_list) { 
 #     lines <- readLines(gzfile(file_path), warn = F)  # gz instead
 #     sample_name <- sub(".*/([^/]+)\\.vcf.gz$", "\\1", file_path)
 #     sample_name <- sub("^VG_", "", sample_name)
 #     entry_count <- length(lines)
 #     
 #     if (sample_name %in% count_data$Sample) {
 #       count_data[count_data$Sample == sample_name, subfolder] <- entry_count
 #     } else {
 #       new_row <- data.frame(Sample = sample_name, 
 #                                SCREENED_MERGED_1 = 0,   SCREENED_MERGED_3 = 0,
 #                             SCREENED_MERGED_6 = 0,   SCREENED_MERGED_ALL = 0,
 ##                             SCREENED_MERGED_1_GBWT = 0,  SCREENED_MERGED_3_GBWT = 0,
  #                            SCREENED_MERGED_6_GBWT = 0,  SCREENED_MERGED_ALL_GBWT = 0,
  #                            stringsAsFactors = F)
   #     new_row[[subfolder]] <- entry_count
   #     count_data <- rbind(count_data, new_row)   }
   # }
 # }
 # return(count_data) }

# VCFtools mutations only
#entry_count_VGonly <- count_entries_VGonly_function(".")
#entry_count_VGonly$Sample  <- gsub(".vcf.gz", "",  entry_count_VGonly$Sample)  
#entry_count_VGonly$Sample  <- gsub("VG_", "",  entry_count_VGonly$Sample) 
#colnames(entry_count_VGonly) <- c("Sample", 
#                            "VG_VG-MAP_1",   "VG_VG-MAP_3", "VG_VG-MAP_6",
#                            "VG_VG-MAP_ALL", "VG_Giraffe_1", "VG_Giraffe_3",
#                            "VG_Giraffe_6",  "VG_Giraffe_ALL")  
                            
# Merge the two data frames by the 'Sample' column
# merged_data7 <- merge(merged_data6, entry_count_VGonly, by = "Sample", all = T)  
# 
# rownames(merged_data7) <- merged_data7$Sample
# merged_data7 <- merged_data7[,-1]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
merged_data7 <- merged_data6
merged_data7 <- merged_data7 %>% mutate(across(-c(  group_id), as.numeric)) 
merged_data7$group_id <- ifelse(merged_data7$group_id == "AMPLICON_SE", "Amplicon", merged_data7$group_id)
merged_data7$group_id <- ifelse(merged_data7$group_id == "AMPLICON_PE", "Amplicon", merged_data7$group_id)

# Save the entry count data to a CSV file
write.csv(merged_data7, "All_SNPs_with_lifted.csv", row.names = T)
write.csv(merged_data7, "All_SNPs_with_lifted.summary.csv", row.names = T)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
#  Now plot the data BCF_VG-MAP_ALL

colnames(merged_data7) <- sub("_M_1", "_Minimap2_1", colnames(merged_data7))
colnames(merged_data7) <- sub("_M_3", "_Minimap2_3", colnames(merged_data7))
colnames(merged_data7) <- sub("_M_6", "_Minimap2_6", colnames(merged_data7))
colnames(merged_data7) <- sub("_M_ALL", "_Minimap2_All", colnames(merged_data7))
colnames(merged_data7) <- sub("_ALL", "-All", colnames(merged_data7))
colnames(merged_data7) <- sub("_All", "-All", colnames(merged_data7))
colnames(merged_data7) <- sub("_6", "-6", colnames(merged_data7))
colnames(merged_data7) <- sub("_3", "-3", colnames(merged_data7))
colnames(merged_data7) <- sub("_1", "-1", colnames(merged_data7))
colnames(merged_data7) <- sub("_BF", "-BF", colnames(merged_data7))

long_data <- merged_data7 %>%
   select(group_id, starts_with("BCF"), starts_with("FB")) %>%
  pivot_longer(cols = -group_id,  names_to = c("Type", "Mapper"),
               names_sep = "_",     values_to = "Value")
long_data2 <- subset(long_data, group_id != "NA")   

long_data2 <- subset(long_data2, group_id != "SIMULATED" & Mapper != "Minimap2-3" & Mapper != "Minimap2-6" & Mapper != "Minimap2-All" & group_id != "Metagenomic") 
pdf("Figure7_Mutation.plot3.pdf", height=6, width=7)    
ggplot(long_data2, aes(Mapper, log10(Value), fill=group_id, colour=group_id)) +
  geom_jitter(alpha=0.2, width=0.4, size=0.3) +geom_boxplot(alpha=0.2)+
  facet_grid( Type ~ group_id, margins=F, switch = "x") +
  theme_minimal()   +  labs(  x = "", y = "Log10-scaled number of mutations") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text.x = element_text(margin = margin(0, 0, 0, 0)),
     strip.background = element_rect(fill="lightgrey", color="grey", linewidth=0),
      strip.placement = "outside")  + ylim(0,4) +
      scale_fill_manual(values=c("#999999", "red", "#E69F00", "#56B4E9", "black"))
dev.off()

long_data3 <-  long_data2  
pdf("Figure7_Mutation.plot4.pdf", height=5, width=6)    
ggplot(long_data3, aes(Mapper, log10(Value), fill=group_id, colour=group_id)) +
  geom_jitter(alpha=0.2, width=0.4, size=0.3) +geom_boxplot(alpha=0.2)+
  facet_grid( Type ~ group_id, margins=F, switch = "x") +
  theme_minimal()   +  labs(  x = "", y = "Log10-scaled number of mutations") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text.x = element_text(margin = margin(0, 0, 0, 0)),
     strip.background = element_rect(fill="lightgrey", color="grey", linewidth=0),
      strip.placement = "outside")  +  ylim(0,4) +
      scale_fill_manual(values=c("#999999", "red", "#E69F00", "#56B4E9", "black"))
dev.off()

# Add rownames as a column named "sample_id" 
merged_data7 <- merged_data7[-127,]
merged_data7 <- merged_data7 %>% rownames_to_column(var = "sample_id")

# Identify columns to pivot, excluding "group_id" and "sample_id"
columns_to_pivot <- merged_data7 %>% select(-group_id, -sample_id) %>% names()

# Convert the identified columns to numeric
merged_data7[columns_to_pivot] <- lapply(merged_data7[columns_to_pivot], as.numeric)

# Identify columns with an underscore in the name
columns_with_separator <- columns_to_pivot[grepl("_", columns_to_pivot)]

# Select only the columns "Sample", "sample_id", and those with an underscore
merged_data7 <- merged_data7 %>% select(group_id, sample_id, all_of(columns_with_separator))

long_data <- merged_data7[,c(1:14,32:43)]  %>%
  pivot_longer(
    cols = -c(group_id, sample_id),  # Exclude 'group_id' and 'sample_id' from pivoting
    names_to = c("Type", "Mapper"),  # New column names for the split components
    names_sep = "_",                 # Separator for splitting the column names
    values_to = "Value"       )       # Name of the new column containing the values
str(long_data)

long_data <- subset(long_data, group_id!="NA")
long_data <- subset(long_data, group_id != "SIMULATED" & group_id != "ONT" & group_id != "Metagenomic" & Mapper != "Minimap2-3" & Mapper != "Minimap2-6" & Mapper != "Minimap2-All" ) 
pdf("Figure_S_mutation_detection_per_sample_arrayed.pdf", height=21, width=13)  
print( ggplot(long_data, aes(x=Mapper, y=log10(Value), fill=group_id, colour=group_id)) +
  geom_jitter(alpha = 0.2, width = 0.4, size = 0.3) +
  geom_point( ) +
  facet_wrap(~Type+group_id+sample_id,ncol=21,scales="free_x", strip.position="bottom") +
  theme_minimal() +
  labs(x = " ", y = "Log10-scaled number of mutations") +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text.x = element_text(margin = margin(0, 0, 0, 0)),
    strip.background = element_rect(fill = "lightgrey", color = "grey", linewidth = 0),
    strip.placement = "outside"  ) )  
dev.off()

pdf("Figure_S_mutation_detection_per_sample_wide.pdf", height=5, width=100)  
print( ggplot(long_data, aes(x=Mapper, y=log10(Value), fill=group_id, colour=group_id)) +
  geom_jitter(alpha = 0.2, width = 0.4, size = 0.3) +
  geom_point( ) +
  facet_grid(Type ~group_id+sample_id, margins=F, switch = "x")  +
  theme_minimal() +
  labs(x = " ", y = "Log10-scaled number of mutations") +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text.x = element_text(margin = margin(0, 0, 0, 0)),
    strip.background = element_rect(fill = "lightgrey", color = "grey", linewidth = 0),
    strip.placement = "outside"  ) )
dev.off()
