library(ggplot2)
library(ggExtra)
library(dplyr)
library(tidyr)
 
metadata <- read.csv("metadata2.csv", sep="\t") # Read metadata
colnames(metadata)=c("Sample", "Library_type")
all_data <- list()

files1 <- c("1", "3", "6", "ALL", "1_GBWT", "3_GBWT", "6_GBWT", "ALL_GBWT", "M") 
for (f in files1){ 
folder <- paste0("SCREENED_MERGED_",f)
print(folder)
all_data <- c()
all_data_het <- c() # het SNPs only
# Loop through each sample in the metadata

for (i in 1:nrow(metadata)) {
  sample_name <- metadata$Sample[i]
  library_type <- metadata$Library[i]
  file_bcf <- paste0(folder,"/BCF_", sample_name, ".vcf.het")
  file_fb <- paste0(folder,"/FB_", sample_name, ".vcf.het")
  file_vg <- paste0(folder,"/VG_", sample_name, ".vcf.het")
  data_bcf = c()
  data_fb = c()
    data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header =F))}
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header =F))  } 
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header =F))  } 
    data <- rbind(data_bcf, data_fb,     data_vg ) 

 if(length(data$V1)!=0){
   data$Sample= rep(sample_name, nrow(data))
   data$Library_type= rep(library_type, nrow(data)) # print(str(data))
   all_data[[i]] <- data # all SNPs 
   all_data_het[[i]] <- subset(data, V4<0.95 & V4 > 0.05) # het SNPs
   } 
} 

combined_data <- bind_rows(all_data)
# dim(combined_data)
colnames(combined_data) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type")
plot1 <- paste0("combined_plots_449_",f,".pdf")
combined_data2 <- combined_data[grep("^449", combined_data$Sample), ]
combined_data3 <- subset(combined_data2, Value>0.05 & Value<0.95)

pdf(plot1, width =3, height = 12)
print( ggplot(combined_data3, aes(x = Position, y = Value, color =Caller)) +
  geom_point(alpha = 0.5, size=1) +
  labs(x="Genomic position", y = "Read-depth allele frequency", color = "Type"  ) +
  xlim(0,150000) + ylim(0,1) + theme_minimal() + facet_wrap(~Sample,ncol=1) )
dev.off() }

 
folder1 <-  "SCREENED_MERGED_1_GBWT"
all_data <- c() 
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$Sample[i]
  library_type <- metadata$Library[i]
  file_bcf <- paste0(folder1,"/BCF_", sample_name, ".vcf.het")
  file_fb <- paste0(folder1,"/FB_", sample_name, ".vcf.het")
  file_vg <- paste0(folder1,"/VG_", sample_name, ".vcf.het")
   data_bcf = c()
  data_fb = c()
  data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header =F))}
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header =F))  } 
  if (file.exists(file_vg )){ try(data_vg <- read.table(file_vg, header =F))  } 
    data <- rbind(data_bcf, data_fb, data_vg)  
  if(length(data$V1)!=0){
   data$Sample= rep(sample_name, nrow(data))
   data$Library_type= rep(library_type, nrow(data)) # print(str(data))
   all_data[[i]] <- data } # all SNPs    
}

folder2 <-  "SCREENED_MERGED_3_GBWT"
all_data2 <- c() 
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$Sample[i]
  library_type <- metadata$Library[i] 
  file_bcf <- paste0(folder2,"/BCF_", sample_name, ".vcf.het")
  file_fb <- paste0(folder2,"/FB_", sample_name, ".vcf.het")
  file_vg <- paste0(folder2,"/VG_", sample_name, ".vcf.het")
   data_bcf = c()
  data_fb = c()
  data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header =F))}
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header =F))  } 
  if (file.exists(file_vg )){ try(data_vg <- read.table(file_vg, header =F))  } 
    data2 <- rbind(data_bcf, data_fb, data_vg)  
  if(length(data2$V1)!=0){
   data2$Sample= rep(sample_name, nrow(data2))
   data2$Library_type= rep(library_type, nrow(data2)) # print(str(data))
   all_data2[[i]] <- data2 } # all SNPs      
}

folder3 <-  "SCREENED_MERGED_6_GBWT"
all_data3 <- c() 
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$Sample[i]
  library_type <- metadata$Library[i] 
  file_bcf <- paste0(folder3,"/BCF_", sample_name, ".vcf.het")
  file_fb <- paste0(folder3,"/FB_", sample_name, ".vcf.het")
  file_vg <- paste0(folder3,"/VG_", sample_name, ".vcf.het")
   data_bcf = c()
  data_fb = c()
  data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header =F))}
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header =F))  } 
  if (file.exists(file_vg )){ try(data_vg <- read.table(file_vg, header =F))  } 
    data2 <- rbind(data_bcf, data_fb, data_vg)  
  if(length(data2$V1)!=0){
   data2$Sample= rep(sample_name, nrow(data2))
   data2$Library_type= rep(library_type, nrow(data2)) # print(str(data))
   all_data3[[i]] <- data2 } # all SNPs      
}

folder4 <-  "SCREENED_MERGED_M"
all_data4 <- c() 
for (i in 1:nrow(metadata)) {
  sample_name <- metadata$Sample[i]
  library_type <- metadata$Library[i] 
  file_bcf <- paste0(folder4,"/BCF_", sample_name, ".vcf.het")
  file_fb <- paste0(folder4,"/FB_", sample_name, ".vcf.het")
  file_vg <- paste0(folder4,"/VG_", sample_name, ".vcf.het")
   data_bcf = c()
  data_fb = c()
  data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header =F))}
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header =F))  } 
  if (file.exists(file_vg )){ try(data_vg <- read.table(file_vg, header =F))  } 
    data2 <- rbind(data_bcf, data_fb, data_vg)  
  if(length(data2$V1)!=0){
   data2$Sample= rep(sample_name, nrow(data2))
   data2$Library_type= rep(library_type, nrow(data2)) # print(str(data))
   all_data4[[i]] <- data2 } # all SNPs      
}

lim1 = 0.05
lim2 = 0.95

combined_data <- bind_rows(all_data, all_data2)  # 1 and 3 GWT 
colnames(combined_data) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type") 
combined_data_6 <- bind_rows(all_data, all_data3) # 1 and 6 GBWT
colnames(combined_data_6) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type")                 
combined_data_m <- bind_rows(all_data4) # Minimap2
colnames(combined_data_m) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type")               
combined_data_all_pvg <- bind_rows(all_data, all_data2)  #  1+ all gbwt
colnames(combined_data_all_pvg) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type") 
dim (unique( subset( combined_data[order(combined_data$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 
dim( unique( subset( combined_data_6[order(combined_data_6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) # 
dim( unique( subset( combined_data_m[order(combined_data_m$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 
dim( unique( subset( combined_data_all_pvg[order(combined_data_all_pvg$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 

combined_data2 <- unique(combined_data[grep("^449", combined_data$Sample), ])
combined_data3 <- subset(combined_data2, Value>lim1 & Value<lim2)

combined_data2_6 <- unique(combined_data_6[grep("^449", combined_data_6$Sample), ])
combined_data3_6 <- subset(combined_data2_6, Value>lim1 & Value<lim2)

combined_data2_m <- combined_data_m[grep("^449", combined_data_m$Sample), ]
combined_data3_m <- subset(combined_data2_m, Value>lim1 & Value<lim2)
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "449_05k_S1", "WGS" ) # dummy dat afor plotting
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "449_07k_S2", "WGS" ) # dummy dat afor plotting
combined_data3_m$Position <- as.numeric(combined_data3_m$Position)
combined_data3_m$Value <- as.numeric(combined_data3_m$Value)
combined_data3_m$Depth <- as.numeric(combined_data3_m$Depth)

combined_data3$Dataset <- "GBWT_1_3"
combined_data3_6$Dataset <- "GBWT_1_6"
combined_data3_m$Dataset <- "Minimap2"
combined_all <- bind_rows(combined_data3, combined_data3_6, combined_data3_m)

# Plot with multi-panel layout
pdf("combined_plots_449_1_3_6_M_GBWT.pdf", width =7, height =7)   
print(ggplot(combined_all, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.5, size = 1) +
    labs(x = "Genomic position", y="Read-depth allele frequency", color="Type") +
    xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
    facet_wrap(~Sample + Dataset, ncol = 3)  )
dev.off()

# count SNPs 449
a3 <- unique(combined_data[grep("^449", combined_data$Sample), ]) 
a6 <- unique(combined_data_6[grep("^449", combined_data_6$Sample), ]) 
am <- unique(combined_data_m[grep("^449", combined_data_m$Sample), ])  
dim (unique( subset( a3[order(a3$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 
dim( unique( subset( a6[order(a6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) # 
dim( unique( subset( am[order(am$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 

# paired data

combined_data <- bind_rows(all_data, all_data2)  # 1 and 3 GWT 
colnames(combined_data) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type") 
combined_data_6 <- bind_rows(all_data, all_data3) # 1 and 6 GBWT
colnames(combined_data_6) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type")                 
combined_data_m <- bind_rows(all_data4) # Minimap2
colnames(combined_data_m) <-c("Caller","Ref","Position","Value","Depth",
                           "Sample", "Library_type") 
                           
combined_data2 <- unique(combined_data[grep("^paired", combined_data$Sample), ])
combined_data3 <- subset(combined_data2, Value>lim1 & Value<lim2) 

combined_data2_6 <- unique(combined_data_6[grep("^paired", combined_data_6$Sample), ])
combined_data3_6 <- subset(combined_data2_6, Value>lim1 & Value<lim2)

combined_data2_m <- unique(combined_data_m[grep("^paired", combined_data_m$Sample), ])
combined_data3_m <- subset(combined_data2_m, Value>lim1 & Value<lim2)
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "paired_dat_e5", "WGS" ) # dummy dat afor plotting
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "paired_dat_e8", "WGS" ) # dummy dat afor plotting
combined_data3_m$Position <- as.numeric(combined_data3_m$Position)
combined_data3_m$Value <- as.numeric(combined_data3_m$Value)
combined_data3_m$Depth <- as.numeric(combined_data3_m$Depth) 

combined_data3$Dataset <- "GBWT_1_3"
combined_data3_6$Dataset <- "GBWT_1_6"
combined_data3_m$Dataset <- "Minimap2"
combined_all_p <- bind_rows(combined_data3, combined_data3_6, combined_data3_m)

pdf("combined_plots_paired_1_3_6_M_GBWT.pdf", width =10, height =11)   
print(ggplot(combined_all_p, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.5, size = 0.5) +
    labs(x = "Genomic position", y="Read-depth allele frequency", color="Type") +
    xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
    facet_wrap(~Sample + Dataset, ncol = 3)  )
dev.off()

# count SNPs, include simulated  
a3 <- unique(combined_data[grep("^paired", combined_data$Sample), ]) 
a6 <- unique(combined_data_6[grep("^paired", combined_data_6$Sample), ]) 
am <- unique(combined_data_m[grep("^paired", combined_data_m$Sample), ])  
dim (unique( subset( a3[order(a3$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 
dim( unique( subset( a6[order(a6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) # 
dim( unique( subset( am[order(am$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) # 

#####################  examples ####################

# SRR26747370
#  SRR26747386

plot_sample_data <- function(sample_name, lim1 = 0.05) { 
  combined_data2 <- unique(combined_data[grep(sample_name, combined_data$Sample), ])
  if (nrow(combined_data2) == 0) {
    stop(paste("No data available in combined_data for sample:", sample_name))}
  combined_data3 <- subset(combined_data2, Value > lim1 & Value < 1)
  combined_data2_6 <- unique(combined_data_6[grep(sample_name, combined_data_6$Sample), ])
  if (nrow(combined_data2_6) == 0) {
    warning(paste("No data available in combined_data_6 for sample:", sample_name))  }
  combined_data3_6 <- subset(combined_data2_6, Value > lim1 & Value < 1)
  combined_data2_m <- unique(combined_data_m[grep(sample_name, combined_data_m$Sample), ])
  if (nrow(combined_data2_m) == 0) {
    warning(paste("No data available in combined_data_m for sample:", sample_name)) }
  combined_data3_m <- subset(combined_data2_m, Value > lim1 & Value < 1) 
  if (nrow(combined_data3) > 0) {
    combined_data3 <- combined_data3 %>% mutate(Dataset = "GBWT_1_3") } 
  if (nrow(combined_data3_6) > 0) {
    combined_data3_6 <- combined_data3_6 %>% mutate(Dataset = "GBWT_1_6") } 
  if (nrow(combined_data3_m) > 0) {
    combined_data3_m <- combined_data3_m %>% mutate(Dataset = "Minimap2") } 
  combined_all_eg <- bind_rows(combined_data3, combined_data3_6, combined_data3_m) 
  if (nrow(combined_all_eg) == 0) {
    stop(paste("No valid data after filtering for sample:", sample_name))  } 
  pdf(paste0(sample_name, ".pdf"), width = 9, height =5)
  print(ggplot(combined_all_eg, aes(x = Position, y = Value, color = Caller)) +
          geom_point(alpha = 0.3, size = 0.1) +
          labs(x = "Genomic position", y = "Read-depth allele frequency", color = "Type") +
          xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
          facet_wrap(~Sample + Dataset, ncol = 1))
  dev.off()
}

# Example usage of the function
plot_sample_data("SRR26747370")
plot_sample_data("SRR12588176")
plot_sample_data("LSDV_NEETHLING_S2")

##################### mixed data #####################

combined_data2 <- unique(combined_data[grep("SRR19090746|SRR19090747|SRR19090748", combined_data$Sample), ])
combined_data3 <- subset(combined_data2, Value>lim1 & Value<lim2) 

combined_data2_6 <- unique(combined_data_6[grep("SRR19090746|SRR19090747|SRR19090748", combined_data_6$Sample), ])
combined_data3_6 <- subset(combined_data2_6, Value>lim1 & Value<lim2)

combined_data2_m <- unique(combined_data_m[grep("SRR19090746|SRR19090747|SRR19090748", combined_data_m$Sample), ])
combined_data3_m <- subset(combined_data2_m, Value>lim1 & Value<lim2)

combined_data3$Dataset <- "GBWT_1_3"
combined_data3_6$Dataset <- "GBWT_1_6"
combined_data3_m$Dataset <- "Minimap2"
combined_all_m <- bind_rows(combined_data3, combined_data3_6, combined_data3_m)

pdf("combined_plots_mix_1_3_6_M_GBWT.pdf", width =18, height =6)   
print(ggplot(combined_all_m, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.4, size = 0.1) +
    labs(x = "Genomic position", y="Read-depth allele frequency", color="Type") +
    xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
    facet_wrap(~Sample + Dataset, ncol = 3)  )
dev.off()

a3 <- unique(combined_data[grep("^SRR19090746|SRR19090747|SRR19090748", combined_data$Sample), ]) 
a6 <- unique(combined_data_6[grep("^SRR19090746|SRR19090747|SRR19090748", combined_data_6$Sample), ]) 
am <- unique(combined_data_m[grep("^SRR19090746|SRR19090747|SRR19090748", combined_data_m$Sample), ])  
dim (unique( subset( a3[order(a3$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  2639
dim( unique( subset( a6[order(a6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  3777
dim( unique( subset( am[order(am$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  825
a_all <- unique(combined_data_all_pvg[grep("^SRR19090746|SRR19090747|SRR19090748", combined_data_all_pvg$Sample), ])  
dim( unique( subset( a_all[order(a_all$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  825

# SRR19090746 only
a3 <- unique(combined_data[grep("^SRR19090746", combined_data$Sample), ]) 
a6 <- unique(combined_data_6[grep("^SRR19090746", combined_data_6$Sample), ]) 
am <- unique(combined_data_m[grep("^SRR19090746", combined_data_m$Sample), ])  
dim (unique( subset( a3[order(a3$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  2271
dim( unique( subset( a6[order(a6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  3363
dim( unique( subset( am[order(am$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  448
a_all <- unique(combined_data_all_pvg[grep("^SRR19090746", combined_data_all_pvg$Sample), ])  
dim( unique( subset( a_all[order(a_all$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  825

# SRR19090747 only
a3 <- unique(combined_data[grep("^SRR19090747", combined_data$Sample), ]) 
a6 <- unique(combined_data_6[grep("^SRR19090747", combined_data_6$Sample), ]) 
am <- unique(combined_data_m[grep("^SRR19090747", combined_data_m$Sample), ])  
dim (unique( subset( a3[order(a3$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  660
dim( unique( subset( a6[order(a6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  799
dim( unique( subset( am[order(am$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  620
a_all <- unique(combined_data_all_pvg[grep("^SRR19090747", combined_data_all_pvg$Sample), ])  
dim( unique( subset( a_all[order(a_all$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  825

# SRR19090748 only
a3 <- unique(combined_data[grep("^SRR19090748", combined_data$Sample), ]) 
a6 <- unique(combined_data_6[grep("^SRR19090748", combined_data_6$Sample), ]) 
am <- unique(combined_data_m[grep("^SRR19090748", combined_data_m$Sample), ])  
dim (unique( subset( a3[order(a3$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  721
dim( unique( subset( a6[order(a6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  1065
dim( unique( subset( am[order(am$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  354
a_all <- unique(combined_data_all_pvg[grep("^SRR19090748", combined_data_all_pvg$Sample), ])  
dim( unique( subset( a_all[order(a_all$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  825

# amplicon data

combined_data2 <- unique(subset(combined_data, Library_type=="AMPLICON"))
combined_data3 <- subset(combined_data2, Value>lim1 & Value<lim2) 

combined_data2_6 <- unique(subset(combined_data_6, Library_type=="AMPLICON"))
combined_data3_6 <- subset(combined_data2_6, Value>lim1 & Value<lim2)

combined_data2_m <- unique(subset(combined_data_m, Library_type=="AMPLICON"))
combined_data3_m <- subset(combined_data2_m, Value>lim1 & Value<lim2)

combined_data3$Dataset <- "GBWT_1_3"
combined_data3_6$Dataset <- "GBWT_1_6"
combined_data3_m$Dataset <- "Minimap2"
combined_all_m <- bind_rows(combined_data3, combined_data3_6, combined_data3_m)

pdf("combined_plots_amplicon_1_3_6_M_GBWT.pdf", width =11, height =33)   
print(ggplot(combined_all_m, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.4, size = 1) +
    labs(x = "Genomic position", y="Read-depth allele frequency", color="Type") +
    xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
    facet_wrap(~Sample + Dataset, ncol = 3)  )
dev.off()
  
dim (unique( subset( combined_data2[order(combined_data2$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  2245
dim( unique( subset( combined_data2_6[order(combined_data2_6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  2249
dim( unique( subset( combined_data2_m[order(combined_data2_m$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  1834


# Metagenomic data

combined_data2 <- unique(subset(combined_data, Library_type=="Metagenomic"))
combined_data3 <- subset(combined_data2, Value>lim1 & Value<lim2) 

combined_data2_6 <- unique(subset(combined_data_6, Library_type=="Metagenomic"))
combined_data3_6 <- subset(combined_data2_6, Value>lim1 & Value<lim2)

combined_data2_m <- unique(subset(combined_data_m, Library_type=="Metagenomic"))
combined_data3_m <- subset(combined_data2_m, Value>lim1 & Value<lim2)

combined_data3$Dataset <- "GBWT_1_3"
combined_data3_6$Dataset <- "GBWT_1_6"
combined_data3_m$Dataset <- "Minimap2"
combined_all_m <- bind_rows(combined_data3, combined_data3_6, combined_data3_m)

pdf("combined_plots_Metagenomic_1_3_6_M_GBWT.pdf", width =8, height =5)   
print(ggplot(combined_all_m, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.8, size = 1) +
    labs(x = "Genomic position", y="Read-depth allele frequency", color="Type") +
    xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
    facet_wrap(~Sample + Dataset, ncol = 3)  )
dev.off()

dim (unique( subset( combined_data2[order(combined_data2$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  57
dim( unique( subset( combined_data2_6[order(combined_data2_6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  62
dim( unique( subset( combined_data2_m[order(combined_data2_m$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  23

# WGS data

combined_data22 <- unique( subset(combined_data, Library_type=="WGS") )
combined_data2 <- unique(combined_data22[!grepl("^paired|^449|SRR19090746|SRR19090747|SRR19090748", combined_data22$Sample), ])
combined_data3 <- subset(combined_data2, Value>lim1 & Value<lim2) 

combined_data22_6 <- unique(subset(combined_data_6, Library_type=="WGS"))
combined_data2_6 <- unique(combined_data22_6[!grepl("^paired|^449|SRR19090746|SRR19090747|SRR19090748", combined_data22_6$Sample), ])
combined_data3_6 <- subset(combined_data2_6, Value>lim1 & Value<lim2)

combined_data22_m <- unique(subset(combined_data_m, Library_type=="WGS"))
combined_data2_m <- unique(combined_data22_m[!grepl("^paired|^449|SRR19090746|SRR19090747|SRR19090748", combined_data22_m$Sample), ])
combined_data3_m <- subset(combined_data2_m, Value>lim1 & Value<lim2)
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "Mongolia2", "WGS" ) # dummy dat afor plotting
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "SRR12021190", "WGS" ) # dummy dat afor plotting
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "SRR18028321", "WGS" ) # dummy dat afor plotting
combined_data3_m[dim(combined_data3_m)[1]+1, ]=c("BCF","KX894508", 1,1.1 ,11,
                       "SRR18612540", "WGS" ) # dummy dat afor plotting 
combined_data3_m$Position <- as.numeric(combined_data3_m$Position)
combined_data3_m$Value <- as.numeric(combined_data3_m$Value)
combined_data3_m$Depth <- as.numeric(combined_data3_m$Depth) 

combined_data3$Dataset <- "GBWT_1_3"
combined_data3_6$Dataset <- "GBWT_1_6"
combined_data3_m$Dataset <- "Minimap2"
combined_all_m <- bind_rows(combined_data3, combined_data3_6, combined_data3_m)

pdf("combined_plots_WGS_1_3_6_M_GBWT.pdf", width =8, height =44)   
print(ggplot(combined_all_m, aes(x = Position, y = Value, color = Caller)) +
    geom_point(alpha = 0.5, size = 1) +
    labs(x = "Genomic position", y="Read-depth allele frequency", color="Type") +
    xlim(0, 150000) + ylim(0, 1) + theme_minimal() +
    facet_wrap(~Sample + Dataset, ncol = 3)  )
dev.off()

dim (unique( subset( combined_data2[order(combined_data2$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  2279
dim( unique( subset( combined_data2_6[order(combined_data2_6$Position, decreasing =F), ], select=c(Ref,Position) ) ) ) #  2090
dim( unique( subset( combined_data2_m[order(combined_data2_m$Position, decreasing =F), ], select=c(Ref,Position) ) )  ) #  1917