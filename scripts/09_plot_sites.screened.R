library(dplyr)
library(ggplot2)
library(ggforce)

# Read metadata file
metadata <- read.csv("metadata2.csv", sep = "\t", header = F, col.names = c("sample_name", "group_id"))

# Select samples belonging to the "AMPLICON_SE" group
lista = c("AMPLICON", "WGS",  "Metagenomic",  "SIMULATED", "all")

empty1 <- data.frame( ) # create empty dataframe
#selected_group <- metadata # set defaults 
bin_rate = 1 # for visualising SNP density
#reject_sites_dir ="UNSCREENED_MERGED_3_GBWT"
 
for (filtername in lista){ 

if(filtername == "AMPLICON"){
   selected_group <- metadata %>% filter(group_id == filtername | group_id == "AMPLICON")
   } else if (filtername == "all") {
   selected_group <- metadata 
   } else { selected_group <- metadata %>% filter(group_id == filtername ) }
   
# Function to read reject sites file
read_reject_sites <- function(file_path) {
  data <- read.table(file_path, header=F,
  col.names = c("genome", "site", "col3", "col4", "col5"))
  data <- data %>% select(genome, site)  # Only keep relevant columns
  return(data) }

# Directory containing the "SCREENED_MERGED_3" files
list1 = c("SCREENED_MERGED_3_GBWT", "SCREENED_MERGED_3",
          "SCREENED_MERGED_6_GBWT", "SCREENED_MERGED_6",
          "SCREENED_MERGED_ALL_GBWT", "SCREENED_MERGED_ALL")
for (reject_sites_dir in list1){
print(reject_sites_dir)
print(filtername)

if(reject_sites_dir =="SCREENED_MERGED_3_GBWT"){ h1=6 }
if(reject_sites_dir =="SCREENED_MERGED_3"){ h1=6 }
if(reject_sites_dir =="SCREENED_MERGED_6_GBWT"){ h1=12 }
if(reject_sites_dir =="SCREENED_MERGED_6"){  h1=12 }
if(reject_sites_dir =="SCREENED_MERGED_ALL_GBWT"){  h1=40 }
if(reject_sites_dir =="SCREENED_MERGED_ALL"){  h1=40 }

# Read all reject sites files for the selected samples
reject_sites_list <- lapply(selected_group$sample_name, function(sample) {
  file_path <- file.path(reject_sites_dir, paste0(sample, ".any.sites"))
  if (file.exists(file_path)) {
    return(read_reject_sites(file_path))
  } else { return(data.frame(genome = character(), site = integer()))  }
})

# Combine all data into one data frame
reject_sites_list <- Filter(function(df) nrow(df) > 0 && ncol(df) > 0, reject_sites_list)
reject_sites_combined <- bind_rows(reject_sites_list)

# Count the frequency of each site per genome
site_counts <- reject_sites_combined %>%
  group_by(genome, site) %>% summarise(frequency = n(), .groups = 'drop')

# Calculate the frequency as a percentage
total_samples <- nrow(selected_group)
site_counts <- site_counts %>%
  mutate(frequency_percentage = (frequency / total_samples) * 100)

# Ensure sites fit within the 1-150000 range and bin the data
site_counts <- site_counts %>%
  filter(site <= 149000) %>% filter(site > 2000) %>%
  mutate(bin = cut(site, breaks = seq(0, 150000, by = bin_rate),
    include.lowest =T, labels =F))

# Summarize the frequency percentages within each bin
binned_site_counts <- site_counts %>%
  group_by(genome, bin) %>%
  summarise(frequency_percentage = sum(frequency_percentage),
     .groups = 'drop') %>%
  mutate(start = (bin - 1) *bin_rate + 1, end = bin *bin_rate) %>%
  mutate(site_range = paste0(start, "-", end))

# Add a column to represent the number of bins per genome
binned_site_counts <- binned_site_counts %>%
  group_by(genome) %>% mutate(n_bins = n()) %>% ungroup()

rd=gsub("SCREENED_MERGED_", "", reject_sites_dir)
empty1 <- rbind(empty1, c(rd, filtername, dim(reject_sites_combined),
             length(unique(binned_site_counts$genome)) ) )

print(rd)
print(filtername)
print(sort(table(reject_sites_combined$genome)))

# Create bar plot
pdf(paste0("site_frequencies_",rd,"_",filtername,".pdf"), height = h1, width =11)
print( ggplot(binned_site_counts, aes(x=start, y=frequency_percentage, fill=genome)) +
  geom_bar(stat = "identity", width =500*bin_rate, alpha = 0.2) +
  labs(x = "Genome position", y = "Mutation frequency") +
  scale_x_continuous(breaks = seq(0, 150000, by = 5000)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        panel.grid.major.x = element_line(color="grey80", linewidth= 0.5, linetype = "dashed"),
        legend.position = "bottom") + facet_wrap(~genome, ncol=1, strip.position = "bottom") )
dev.off()}
}
colnames(empty1) = c("Source", "Library_type", "New_mutations", "Other_genomes")
write.csv(empty1, "Rates_new_mutations.csv", row.names=F)
