library(ggplot2)
library(ggExtra)
library(dplyr)
library(tidyr)
 
# Read metadata
metadata <- read.csv("metadata2.csv", sep="\t")
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
  file_vg <- paste0(folder,"/VG_", sample_name, ".vcf.het")
  file_fb <- paste0(folder,"/FB_", sample_name, ".vcf.het")
  data_bcf = c()
  data_fb = c()
   data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header = FALSE))}
   if (file.exists(file_vg )){ try(data_vg <- read.table(file_vg, header = FALSE))  }
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header = FALSE))  } 
   data <- rbind(data_bcf, data_vg, data_fb)  
  if(length(data$V1)!=0){
   data$Sample= rep(sample_name, nrow(data))
   data$Library_type= rep(library_type, nrow(data))
  # print(str(data))
   all_data[[i]] <- data # all SNPs 
   all_data_het[[i]] <- subset(data, V4<1.95 & V4 > 0.05) # all SNPs 
   # all_data_het[[i]] <- subset(data, V4<0.95 & V4 > 0.05) # het SNPs
   } 
}

combined_data_het <- bind_rows(all_data_het)
colnames(combined_data_het) <-c("Caller","Ref","Position","Value","Depth", "Sample", "Library_type") #  data to plot 
combined_data2h <- subset(combined_data_het, Library_type == "WGS" |
               Library_type== "AMPLICON" | Library_type == "Metagenomic")
unwanted_prefixesh <- c("paired", "449", "SRR233", "SRR255")
cc2h <- combined_data2h[!Reduce(`|`, lapply(unwanted_prefixesh, startsWith,
                  x=combined_data2h$Sample)), ]
combined_data3h <- subset(cc2h, Value>0.05 & Value<0.95)
combined_data3h <- distinct(combined_data3h)

plot11 <- paste0("combined_freq_",f,".pdf")
print(plot11)
pdf(plot11, width =18, height = 10)
print( ggplot(combined_data3h, aes(x = Value)) +
  geom_histogram(data=subset(combined_data3h, Caller=="FB"), alpha = 0.4,
           binwidth =0.02, position = "identity", fill="blue") +
  geom_histogram(data=subset(combined_data3h, Caller=="BCF"), alpha = 0.4,
           binwidth =0.02, position = "identity", fill="red") +
    geom_histogram(data=subset(combined_data3h, Caller=="VG"), alpha = 0.4,
           binwidth =0.02, position = "identity", fill="green") +
   geom_density() +  facet_wrap(~ Sample, scales = "free") +  theme_bw() +
  labs( x = "Value",  y = "Frequency") + xlim(0,1) + 
  theme(legend.position = "bottom") )
dev.off()

combined_data <- bind_rows(all_data)
dim(combined_data)
colnames(combined_data) <-c("Caller","Ref","Position","Value","Depth", "Sample", "Library_type")
plot1 <- paste0("combined_plots_",f,".pdf")
combined_data2 <- subset(combined_data, Library_type == "WGS" |
	       Library_type== "AMPLICON" | Library_type == "Metagenomic")
unwanted_prefixes <- c("paired", "449", "SRR233", "SRR255")
cc2 <- combined_data2[!Reduce(`|`, lapply(unwanted_prefixes, startsWith,
 		  x=combined_data2$Sample)), ]
combined_data3 <- subset(cc2, Value>0.05 & Value<0.95)

print(plot1)
pdf(plot1, width =20, height = 15)
print( ggplot(combined_data3, aes(x = Position, y = Value, color =Caller)) +
  geom_point(alpha = 0.3, size=0.3) +
  labs(    x = "Position", y = "Value", color = "Type"  ) +
  xlim(0,150000) + ylim(0,1) + theme_minimal() + facet_wrap(~Sample,ncol=4) +
  geom_text(data = combined_data3 %>% distinct(Sample, Library_type),
            aes(x = 20000, y = 0.5, label = Library_type),
            hjust = 1, vjust = 1, size = 3, inherit.aes = F) )
dev.off() }


#  Now do het SNPs only

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
  file_vg <- paste0(folder,"/VG_", sample_name, ".vcf.het")
  file_fb <- paste0(folder,"/FB_", sample_name, ".vcf.het")
  data_bcf = c()
  data_fb = c()
   data_vg = c()
  if (file.exists(file_bcf)){ try(data_bcf <- read.table(file_bcf, header = FALSE))}
   if (file.exists(file_vg )){ try(data_vg <- read.table(file_vg, header = FALSE))  }
  if (file.exists(file_fb )){ try(data_fb <- read.table(file_fb, header = FALSE))  } 
   data <- rbind(data_bcf, data_vg, data_fb)  
  if(length(data$V1)!=0){
   data$Sample= rep(sample_name, nrow(data))
   data$Library_type= rep(library_type, nrow(data))
  # print(str(data))
   all_data[[i]] <- data # all SNPs 
all_data_het[[i]] <- subset(data, V4<0.95 & V4 > 0.05) # het SNPs
   } 
}

combined_data_het <- bind_rows(all_data_het)
colnames(combined_data_het) <-c("Caller","Ref","Position","Value","Depth", "Sample", "Library_type") #  data to plot 
combined_data2h <- subset(combined_data_het, Library_type == "WGS" |
               Library_type== "AMPLICON" | Library_type == "Metagenomic")
unwanted_prefixesh <- c("paired", "449", "SRR233", "SRR255")
cc2h <- combined_data2h[!Reduce(`|`, lapply(unwanted_prefixesh, startsWith,
                  x=combined_data2h$Sample)), ]
combined_data3h <- subset(cc2h, Value>0.05 & Value<0.95)
combined_data3h <- distinct(combined_data3h)

plot11 <- paste0("combined_freq_",f,".hetonly.pdf")
print(plot11)
pdf(plot11, width = 18, height = 10)
print( ggplot(combined_data3h, aes(x = Value)) +
  geom_histogram(data=subset(combined_data3h, Caller=="FB"), alpha = 0.4,
           binwidth =0.02, position = "identity", fill="blue") +
  geom_histogram(data=subset(combined_data3h, Caller=="BCF"), alpha = 0.4,
           binwidth =0.02, position = "identity", fill="red") +
    geom_histogram(data=subset(combined_data3h, Caller=="VG"), alpha = 0.4,
           binwidth =0.02, position = "identity", fill="green") +
   geom_density() +  facet_wrap(~ Sample, scales = "free") +  theme_bw() +
  labs( x = "Value",  y = "Frequency") + xlim(0,1) + 
  theme(legend.position = "bottom") )
dev.off()

combined_data <- bind_rows(all_data)
dim(combined_data)
colnames(combined_data) <-c("Caller","Ref","Position","Value","Depth", "Sample", "Library_type")
plot1 <- paste0("combined_plots_",f,".hetonly.pdf")
combined_data2 <- subset(combined_data, Library_type == "WGS" |
	       Library_type== "AMPLICON" | Library_type == "Metagenomic")
unwanted_prefixes <- c("paired", "449", "SRR233", "SRR255")
cc2 <- combined_data2[!Reduce(`|`, lapply(unwanted_prefixes, startsWith,
 		  x=combined_data2$Sample)), ]
combined_data3 <- subset(cc2, Value>0.05 & Value<0.95)

print(plot1)
pdf(plot1, width = 20, height = 15)
print( ggplot(combined_data3, aes(x = Position, y = Value, color =Caller)) +
  geom_point(alpha = 0.3, size=0.3) +
  labs(    x = "Position", y = "Value", color = "Type"  ) +
  xlim(0,150000) + ylim(0,1) + theme_minimal() + facet_wrap(~Sample,ncol=4) +
  geom_text(data = combined_data3 %>% distinct(Sample, Library_type),
            aes(x = 20000, y = 0.5, label = Library_type),
            hjust = 1, vjust = 1, size = 3, inherit.aes = F) )
dev.off() }