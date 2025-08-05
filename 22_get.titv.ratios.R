# Load necessary library
library(dplyr)

# Read the sample list from 'acc_list.sim.no_ont'
# Assuming no header and first column is sample name, second is library type
samples <- read.table("metadata2.csv", header = FALSE, stringsAsFactors = F)
colnames(samples) <- c("Sample", "Library_type")

# Initialize an empty dataframe to store the results
results <- data.frame(Sample = character(),
                      Total_Changes = integer(),
                      A_to_G = integer(),
                      A_to_C = integer(),
                      A_to_T = integer(),
                      G_to_C = integer(),
                      G_to_T = integer(),
                      T_to_C = integer(),
                      Ti = integer(),
                      Tv = integer(),
                      Ti_Tv_Ratio = numeric(),
                      stringsAsFactors = FALSE)

# Define the function to parse mutation changes
parse_mutations <- function(file_path) {
  # Read the mutation file (assuming no header)
  mutations <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(mutations) <- c("SeqID", "Position", "Ref", "Alt", "Quality")
  
  # Calculate mutation types
  A_to_G <- nrow(mutations %>% filter(Ref == "A" & Alt == "G"))
  A_to_C <- nrow(mutations %>% filter(Ref == "A" & Alt == "C"))
  A_to_T <- nrow(mutations %>% filter(Ref == "A" & Alt == "T"))
  G_to_C <- nrow(mutations %>% filter(Ref == "G" & Alt == "C"))
  G_to_T <- nrow(mutations %>% filter(Ref == "G" & Alt == "T"))
  T_to_C <- nrow(mutations %>% filter(Ref == "T" & Alt == "C"))
  
  # Total number of changes
  total_changes <- A_to_G + A_to_C + A_to_T + G_to_C + G_to_T + T_to_C
  
  # Calculate Ti (Transitions) and Tv (Transversions)
  Ti <- A_to_G + T_to_C
  Tv <- A_to_C + A_to_T + G_to_C + G_to_T
  
  # Calculate Ti/Tv ratio
  Ti_Tv_Ratio <- ifelse(Ti > 0 & Tv > 0, Ti / Tv, NA)
  
  return(list(Total_Changes = total_changes,
              A_to_G = A_to_G, A_to_C = A_to_C, A_to_T = A_to_T,
              G_to_C = G_to_C, G_to_T = G_to_T, T_to_C = T_to_C,
              Ti = Ti, Tv = Tv, Ti_Tv_Ratio = Ti_Tv_Ratio))
}

# Loop through each sample in the list and process its corresponding file
for (i in 1:nrow(samples)) {
  sample_name <- samples$Sample[i]
  # Construct the path to the mutation file (adjust the path as needed)
  file_path <- file.path("SCREENED_MERGED_3_GBWT", paste0(sample_name, ".sites"))
  
  # Parse the mutations for the current sample
  mutation_data <- parse_mutations(file_path)
  
  # Add the data to the results dataframe
  results <- rbind(results, data.frame(Sample = sample_name,
                                       Total_Changes = mutation_data$Total_Changes,
                                       A_to_G = mutation_data$A_to_G,
                                       A_to_C = mutation_data$A_to_C,
                                       A_to_T = mutation_data$A_to_T,
                                       G_to_C = mutation_data$G_to_C,
                                       G_to_T = mutation_data$G_to_T,
                                       T_to_C = mutation_data$T_to_C,
                                       Ti = mutation_data$Ti,
                                       Tv = mutation_data$Tv,
                                       Ti_Tv_Ratio = mutation_data$Ti_Tv_Ratio,
                                       stringsAsFactors = FALSE))
}

# Display the results dataframe
print(results)

# Optionally, save the results to a CSV file
write.csv(results, "mutation_analysis_results.csv", row.names = FALSE)
