library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(ggrepel)

################## now for Giraffe 1+3 combined ##################

data <- read.table('rates.new3.txt', header=T, stringsAsFactors=F)
colnames(data)[3]="PVG2"
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = all - M, 
         Intersection = M + PVG2 - all,
         Minimap2 = all - PVG2)  
long_data <- filtered_data %>%
  gather(key="Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type=factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data2 <- long_data
long_data3 <- long_data
long_data3 <- long_data3 %>% filter(Library_type!='SIMULATED')# Illumian real only

write.csv(long_data2, "long_data2.csv",  row.names=F)

long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1", Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1", Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1", Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c() 

for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
    data_subset$Type <- factor(data_subset$Type,
    levels = c("Minimap2", "PVG", "Intersection"))   
  numbers1= c(numbers1, length(data_subset$PVG) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>% group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(all)) %>%
    arrange(desc(Total_Minimap2)) %>% pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  }  

pdf("snp_overlap.giraffe.13.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
g1 <- annotate_figure(p1, bottom  = text_grob("Giraffe 1- & 3-sample PVGs", 
               color = "blue", face = "bold", size = 14))
print(g1)
dev.off()

long_data3$Diff = long_data3$all - long_data3$PVG2
table(subset(long_data3, Diff>0 & Type=="Minimap2")$Library_type)
table(subset(long_data3, Diff==0 & Type=="Minimap2")$Library_type)

pp1 <- long_data3 %>% filter(Type == "PVG") %>%
  ggplot(aes(x =  (Diff), y =  (Value), color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)   +  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1- & 3-sample PVGs") +
  theme_minimal() +  theme(legend.position = "right") + ylim(0,4000) 

pp2 <- long_data3 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)  +   # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size =3,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.3,   # Padding around text boxes
                  point.padding = 0.4,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1- & 3-sample PVGs") + 
    theme(legend.position = "right")

inset <- long_data3 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + ylim(0,1900) + xlim(0,31) +  
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 3,             # Text size
                  max.overlaps = 6,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.6,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels 
  labs(x = " ", y =  "") + theme(legend.position = "none") + 
    theme(legend.position = "none")
  
pdf("relative_rate_snps.giraffe1_3.pdf", width=7.5, height=5.5)
print(pp2 + inset_element(inset,  0.475, 0.3, 1, 1) )
dev.off()

################## now for Giraffe 1+6 combined ##################

data <- read.table('rates.new6.txt', header=T, stringsAsFactors=F)
colnames(data)[3]="PVG2"
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = all - M, 
         Intersection = M + PVG2 - all,
         Minimap2 = all - PVG2) 
long_data <- filtered_data %>%
  gather(key="Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type=factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample)  
long_data2 <- long_data %>% filter(Library_type!='SIMULATED')# Illumian real only

write.csv(long_data2, "long_data2.giraffe6.csv",  row.names=F)

long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1", Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1", Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1", Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
      data_subset$Type <- factor(data_subset$Type, levels = c("Minimap2", "PVG", "Intersection"))  # Specify the desired order here
  numbers1= c(numbers1, length(data_subset$PVG) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>% group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  }  

pdf("snp_overlap.giraffe16.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
g1 <- annotate_figure(p1, bottom  = text_grob("Giraffe 1- & 6-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(g1)
dev.off()

long_data2$Diff = long_data2$all - long_data2$PVG
table(subset(long_data2, Diff>0 & Type=="Minimap2")$Library_type)
table(subset(long_data2, Diff==0 & Type=="Minimap2")$Library_type)  

pp1 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x =  (Diff), y =  (Value), color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)   +  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1- & 6-sample PVGs") +
  theme_minimal() +  theme(legend.position = "right")

pp2 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + xlim(0,62) +ylim(0,1600) +  
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size =3,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.3,   # Padding around text boxes
                  point.padding = 0.4,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = " ", y =  "") + theme(legend.position = "none")

inset <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 3,             # Text size
                  max.overlaps = 6,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.6,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1- & 6-sample PVGs") +
    theme(legend.position = "right")
  
pdf("relative_rate_snps.giraffe1_6.pdf", width=7.5, height=5.5)
print(inset + inset_element(pp2,  0.45, 0.3, 1, 1) )
dev.off()

################## now for VG-MAP 1 ##################
# Sample  Library_type PVG_BCF  PVG_FB  PVG_VG PVG Minimap2   
data <- read.table('rates.new3.vgmap1.txt', header=T, stringsAsFactors=F)
colnames(data)[6]="PVG2"
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG2 - Total_all,
         Minimap2 = Total_all - PVG2) 
long_data <- filtered_data %>%
  gather(key="Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type=factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data2 <- long_data %>% filter(Library_type!='SIMULATED')

long_data <- long_data %>%
  mutate(Library_type = case_when( 
    Sample=="SRR26747370" ~ "AMPLICON_1",
    Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1",
    Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1",
    Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG2) ) 
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>%
    group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.vgmap1.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
v1 <- annotate_figure(p1, bottom = text_grob("VG-MAP 1-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(v1)
dev.off()

long_data2$Diff = long_data2$Total_all - long_data2$PVG2
table(subset(long_data2, Diff>0 & Type=="Minimap2")$Library_type)
table(subset(long_data2, Diff==0 & Type=="Minimap2")$Library_type)

pp1 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x =  (Diff), y =  (Value), color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)   +  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 25,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1-sample PVG VG-MAP") +
  theme_minimal() +  theme(legend.position = "right")

pp2 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + xlim(0,32) +ylim(0,12)+ # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size =3,             # Text size
                  max.overlaps = 25,    # Allow some overlaps
                  box.padding = 0.3,   # Padding around text boxes
                  point.padding = 0.4,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = " ", y =  "") + theme(legend.position = "none")

inset <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 3,             # Text size
                  max.overlaps = 6,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.6,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1-sample PVG VG-MAP") +
    theme(legend.position = "right")
  
pdf("relative_rate_snps.vgmap1.pdf", width=7.5, height=5.5)
print(inset + inset_element(pp2,  0.475, 0.3, 1, 1) )
dev.off()

################## now for VG-MAP 3 ##################
# Sample  Library_type    PVG_BCF  PVG_FB  PVG Minimap2  Total_BCF   Total_FB  Total_all
#                         pb      pf       p     m       pb+m-pbm   pf+m-pfm   p+m-pm
data <- read.table('rates.new3.vgmap3.txt', header=T, stringsAsFactors=F)
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG - Total_all,
         Minimap2 = Total_all - PVG) 
long_data <- filtered_data %>%
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1",
    Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1",
    Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1",
    Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>%
    group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.vgmap3.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
v3 <- annotate_figure(p1, bottom  = text_grob("VG-MAP 3-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(v3)
dev.off()

################## now for VG-MAP 6 ################## 

data <- read.table('rates.new3.vgmap6.txt', header=T, stringsAsFactors=F)
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG - Total_all,
         Minimap2 = Total_all - PVG) 
long_data <- filtered_data %>%
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1",
    Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1",
    Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1",
    Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>% group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.vgmap6.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
v6 <- annotate_figure(p1, bottom  = text_grob("VG-MAP 6-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(v6)
dev.off()

################## now for VG-MAP all ################## 

data <- read.table('rates.new3.vgmapall.txt', header=T, stringsAsFactors=F)
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG - Total_all,
         Minimap2 = Total_all - PVG) 
long_data <- filtered_data %>%
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1",
    Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1",
    Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1",
    Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>% group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.vgmapall.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
vall <- annotate_figure(p1, bottom  = text_grob("VG-MAP all-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(vall)
dev.off()

################## now for Giraffe 1 ##################

data <- read.table('rates.new3.giraffe1.txt', header=T, stringsAsFactors=F)
colnames(data)[6]="PVG2"
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG2 - Total_all,
         Minimap2 = Total_all - PVG2) 
long_data <- filtered_data %>%
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Minimap2","Intersection")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data2 <- long_data %>% filter(Library_type!='SIMULATED')
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1", Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1", Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1", Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c() 

# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG2) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>%
    group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("cyan", "red", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  }  

pdf("snp_overlap.giraffe1.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
g1 <- annotate_figure(p1, bottom  = text_grob("Giraffe 1-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(g1)
dev.off()

long_data2$Diff = long_data2$Total_all - long_data2$PVG2
table(subset(long_data2, Diff>0 & Type=="Minimap2")$Library_type)
table(subset(long_data2, Diff==0 & Type=="Minimap2")$Library_type)

pp1 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x =  (Diff), y =  (Value), color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)   +  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 35,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1-sample PVG") +
  theme_minimal() +  theme(legend.position = "right")

pp2 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) +  xlim(0,32) +ylim(0,12)+  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 35,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x =" ", y = " ") +   theme(legend.position = "none")

inset <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 3,             # Text size
                  max.overlaps = 6,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.6,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 1-sample PVG") +
  theme(legend.position = "none")
  
pdf("relative_rate_snps.giraffe1.pdf", width=5.5, height=5.5)
print(inset + inset_element(pp2,  0.475, 0.3, 1, 1) )
dev.off() 

################## now for Giraffe 3 ##################

data <- read.table('rates.new3.giraffe3.txt', header=T, stringsAsFactors=F)
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
colnames(filtered_data)[6] = "PVG2"
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG2 - Total_all,
         Minimap2 = Total_all - PVG2) 
long_data <- filtered_data %>%
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data2 <- long_data %>% filter(Library_type!='SIMULATED')
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1", Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1", Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1", Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
    
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG2) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>%
    group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.giraffe3.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
g3 <- annotate_figure(p1, bottom  = text_grob("Giraffe 3-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(g3)
dev.off()

long_data2$Diff = long_data2$Total_all - long_data2$PVG2
table(subset(long_data2, Diff>0 & Type=="Minimap2")$Library_type)
table(subset(long_data2, Diff==0 & Type=="Minimap2")$Library_type) 

pp1 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x =  (Diff), y =  (Value), color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)   +  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 3-sample PVG") +
  theme_minimal() +  theme(legend.position = "right")

pp2 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4)  + xlim(0,100) +ylim(0,2000)+  
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "", y = "") +      theme(legend.position = "none")

inset <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) +  # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 3,             # Text size
                  max.overlaps = 6,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.6,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 3-sample PVG") +
   theme(legend.position = "none")
  
pdf("relative_rate_snps.giraffe3.pdf", width=6, height=6.1)
print(inset + inset_element(pp2,  0.475, 0.3, 1, 1) )
dev.off()
  
################## now for Giraffe 6 ##################

data <- read.table('rates.new3.giraffe6.txt', header=T, stringsAsFactors=F)
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT'
colnames(filtered_data)[6] = "PVG2" 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG2 - Total_all,
         Minimap2 = Total_all - PVG2) 
long_data <- filtered_data %>%  
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data2 <- long_data %>% filter(Library_type!='SIMULATED')
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1", Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1", Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1", Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
    
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG2) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>%
    group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.giraffe6.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
g6 <- annotate_figure(p1, bottom  = text_grob("Giraffe 6-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(g6)
dev.off()

long_data2$Diff = long_data2$Total_all - long_data2$PVG2
table(subset(long_data2, Diff>0 & Type=="Minimap2")$Library_type)
table(subset(long_data2, Diff==0 & Type=="Minimap2")$Library_type) 

pp2 <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + xlim(0,140) +ylim(0,410) +    
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 4,             # Text size
                  max.overlaps = 10,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.5,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "", y = "") +
     theme(legend.position = "none")

inset <- long_data2 %>% filter(Type == "PVG") %>%
  ggplot(aes(x = Diff, y =Value, color = Library_type)) +
  geom_point(size = 2, alpha = 0.4) + # Scatter plot
  geom_text_repel(aes(label = Sample),  # Adding labels for samples
                  size = 3,             # Text size
                  max.overlaps = 6,    # Allow some overlaps
                  box.padding = 0.35,   # Padding around text boxes
                  point.padding = 0.6,  # Padding around points
                  segment.size = 0.2) + # Line size connecting points to labels
  labs(x = "SNPs unique to Minimap2", y = "SNPs unique to 6-sample PVG") +
   theme(legend.position = "none")
  
pdf("relative_rate_snps.giraffe6.pdf", width=6, height=6.1)
print(inset + inset_element(pp2,  0.475, 0.3, 1, 1) )
dev.off()

################## now for Giraffe all ##################

data <- read.table('rates.new3.giraffeall.txt', header=T, stringsAsFactors=F)
filtered_data <- data %>% filter(Library_type != 'ONT'  ) # Filter out 'ONT' 
filtered_data <- filtered_data %>%
  mutate(PVG = Total_all - Minimap2, 
         Intersection = Minimap2 + PVG - Total_all,
         Minimap2 = Total_all - PVG) 
long_data <- filtered_data %>%
  gather(key = "Type", value="Value", PVG,Intersection,Minimap2) %>%
  mutate(Type = factor(Type, levels=c("PVG","Intersection","Minimap2")))
long_data$Sample <- gsub("LSDV_", "", long_data$Sample)
long_data$Sample <- gsub("TENAPI_", "", long_data$Sample) 
long_data <- long_data %>%
  mutate(Library_type = case_when(
    Sample=="SRR26747370" ~ "AMPLICON_1", Sample=="SRR12588176" ~ "AMPLICON_1",
    Sample=="SRR15145273" ~ "AMPLICON_1", Sample=="SRR15145275" ~ "AMPLICON_1",
    Sample=="SRR15145277" ~ "AMPLICON_1", Sample=="SRR15145279" ~ "AMPLICON_1",
    Sample=="SRR19090746" ~ "Mix", Sample=="SRR19090747" ~ "Mix",
    Sample=="SRR19090748" ~ "Mix",
    T ~ Library_type ))
plots <- list()
numbers1=c()
# Loop through each unique Library_type to create individual plots
for(lib_type in unique(long_data$Library_type)) { 
  data_subset <- long_data %>% filter(Library_type == lib_type) 
  numbers1= c(numbers1, length(data_subset$PVG) )
  sample_order <- data_subset %>%
    filter(Type == "Minimap2") %>%
    group_by(Sample) %>%
    summarise(Total_Minimap2 = sum(Total_all)) %>%
    arrange(desc(Total_Minimap2)) %>%
    pull(Sample) 
  data_subset$Sample <- factor(data_subset$Sample, levels = sample_order) 
  lib_type1 = gsub("_1", "", lib_type)
  lib_type1 = gsub("_2", "", lib_type1) 
  p <- ggplot(data_subset, aes(x=Sample, y=Value, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("red", "cyan", "black")) +   
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x="", y="Number of SNPs", title=paste("", lib_type1))  
  plots[[lib_type]] <- p  } 

pdf("snp_overlap.giraffeall.pdf", width=21, height=5) 
p1 <- ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend =T,
          widths = (numbers1/3 + 4), heights = rep(5, 6))
gall <- annotate_figure(p1, bottom  = text_grob("Giraffe all-sample PVG", 
               color = "blue", face = "bold", size = 14))
print(gall)
dev.off()

pdf("snp_overlap.pdf", width=21, height=35) 
grid.arrange( g1,g3,g6,gall,v1,v3,v6,vall, ncol=1, nrow=8)
dev.off()

pdf("snp_overlap.giraffe.pdf", width=20, height=8) 
grid.arrange( g1,g3,  ncol=1, nrow=2)
dev.off()