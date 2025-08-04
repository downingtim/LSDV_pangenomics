# Required libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# 1. Read the list of samples
samples <- readLines("acc_list.sim.no_ont") 
subfolders <- c("LSDV1", "LSDVG", "LSDVVG",   "LSDVG_3", "LSDVVG_3",
               "LSDVG_6", "LSDVVG_6",   "LSDVG_ALL", "LSDVVG_ALL") 
coverage_data <- list() 
for (folder in subfolders) {
  for (sample in samples) { 
    file_path <- file.path(folder, "COVERAGE", paste0(sample, ".coverage.txt")) 
    if (file.exists(file_path)) { 
      cov_data <- read_tsv(file_path, col_names = T, show_col_types = F)
      cov_data <- cov_data %>%
        select(rname = `#rname`, meandepth, meanbaseq, meanmapq) %>%
        rename(depth = meandepth, BQ = meanbaseq, MQ = meanmapq)
      # Sum data for rname 'KX894507' and 'KX894508' and label them as 'KX894508'
      cov_data <- cov_data %>%
        mutate(rname = ifelse(rname == "KX894507", "KX894508", rname)) %>% 
        group_by(rname) %>%
        summarize(across(c(depth, BQ, MQ), sum, na.rm = T)) 
     cov_data <- cov_data %>%
        mutate(Sample = sample, Folder = folder)
      coverage_data[[paste0(folder, "_", sample)]] <- cov_data }  
  }
}
   
comparison_data <- bind_rows(coverage_data)
metadata <- read.csv("metadata2.csv", sep="\t", header = F) 
colnames(metadata) <- c("Sample", "Library_type") 
comparison_data <- merge(comparison_data, metadata, by = "Sample", all.x = T)
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVG_3", "Giraffe 3", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVVG_3", "VG-MAP 3", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVG_6", "Giraffe 6", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVVG_6", "VG-MAP 6", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVG_ALL", "Giraffe all", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVVG_ALL", "VG-MAP all", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDV1", "Minimap2", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVG", "Giraffe 1", x) }))
comparison_data <- data.frame(lapply(comparison_data, function(x) { gsub("LSDVVG", "VG-MAP 1", x) })) 
comparison_data$depth <- as.numeric( comparison_data$depth)
comparison_data$BQ <- as.numeric( comparison_data$BQ)
comparison_data$MQ <- as.numeric( comparison_data$MQ)
write.csv(comparison_data, "coverage1.csv")
color_palette <- brewer.pal(3, "Set2")


comparison_data1<-subset(comparison_data, Folder=="Minimap2" | Folder=="Giraffe 1" | Folder=="VG-MAP 1" | Folder=="Giraffe 3"   | Folder=="VG-MAP 3"| Folder=="Giraffe 6"   | Folder=="VG-MAP 6"  |
  Folder=="Giraffe all" |Folder=="VG-MAP all"  )
  
summed_depths <- comparison_data1 %>%
  group_by(Folder, Sample ) %>%         # Group by Folder and Sample
  summarise(total_depth = sum(depth, na.rm = T))
print(summed_depths %>% group_by(Folder) %>% 
   summarise(median_depth = median(total_depth) ) )

comparison_data1 %>% filter (rname=="KX894508") %>% group_by(Folder, Sample) %>% 
     group_by(Folder) %>% summarise( BQ = median(BQ))   
comparison_data1 %>% filter (rname=="KX894508") %>% group_by(Folder, Sample) %>% 
     group_by(Folder) %>% summarise( MQ = median(MQ)) 
     
#  ignore all samples, not useful
     
comparison_data1<-subset(comparison_data, Folder=="Minimap2" | Folder=="Giraffe 1" | Folder=="VG-MAP 1" | Folder=="Giraffe 3"   | Folder=="VG-MAP 3"| Folder=="Giraffe 6"   | Folder=="VG-MAP 6"   )

summed_depths2 <- comparison_data1  %>% 
  group_by(Folder, Sample, Library_type )  %>% 
  summarise_at(.vars = vars(depth,BQ,MQ), .funs = c(median="median"))

pdf("Depth_BQ_1.pdf", width=10, height=4)
ggplot(comparison_data1, aes(x = log10(depth + 1.1), y = BQ, colour = Library_type)) +
  geom_point(alpha = 0.5, size = 1.1) +    
  facet_wrap(~Folder, nrow = 1) +   ylim(24, 101) + 
  labs(x = "Log10-scaled read depth", y = "BQ") +   theme_minimal() + 
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 1),
        strip.text = element_text(size = 8)) + 
  scale_fill_manual(values = brewer.pal(7, "Set2")) + 
  geom_smooth(method = "lm", se =T, color = "grey", linetype = "dashed", alpha = 0.5) + 
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = 'left', label.y.npc = 'top')
dev.off()
  
pdf("Depth_MQ_1.pdf", width=10, height=4)
ggplot(comparison_data1, aes(x = log10(depth + 1.1), y =MQ, colour = Library_type)) +
  geom_point(alpha = 0.5, size = 1.1) +    
  facet_wrap(~Folder, nrow = 1) +    
  labs(x = "Log10-scaled read depth", y = "MQ") +   theme_minimal() + 
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 1),
        strip.text = element_text(size = 8)) + 
  scale_fill_manual(values = brewer.pal(7, "Set2")) + 
  geom_smooth(method = "lm", se =T, color = "grey", linetype = "dashed", alpha = 0.5) + 
  stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = 'left', label.y.npc = 'top')
dev.off()


  # Plot for depth
  pdf("Depth_1.pdf", width=30, height=13)
   ggplot(comparison_data1, aes(x = Folder, y = depth, fill = rname)) +
  geom_bar(stat = "identity", position = "stack") +    
  facet_wrap(~ Library_type + Sample, ncol = 12, scales = "free_y") +   
  labs(x = "Folder", y = "Mean Read Depth") +   theme_minimal() + 
  theme( axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
    strip.text = element_text(size = 8) ) + 
  scale_fill_manual(values = brewer.pal(6, "Set2")) 
  dev.off()

  # Plot for BQ
  pdf( "BQ_1.pdf", width=23, height=33)
   ggplot(comparison_data1, aes(x = Folder, y = BQ, fill = rname)) +
  geom_bar(stat = "identity") +    
  facet_wrap(~ Library_type + Sample +rname, ncol = 24, scales = "free_y") +   
  labs(x = "Folder", y = "Mean BQ") +   theme_minimal() + 
  theme( axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
    strip.text = element_text(size = 5)  ) + 
  scale_fill_manual(values = brewer.pal(6, "Set2")) 
  dev.off()

  # Plot for MQ
  pdf( "MQ_1.pdf", width=23, height=33)
   ggplot(comparison_data1, aes(x = Folder, y = MQ, fill = rname)) +
  geom_bar(stat = "identity" ) +    
  facet_wrap(~ Library_type + Sample +rname, ncol = 24, scales = "free_y") +   
  labs(x = "Folder", y = "Mean MQ") +   theme_minimal() + 
  theme( axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
    strip.text = element_text(size = 5)  ) + 
  scale_fill_manual(values = brewer.pal(6, "Set2")) 
  dev.off()
 



