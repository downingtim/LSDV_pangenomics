library(ggplot2)
library(dplyr)
 
# Function to read data
read_data <- function(file_path) {
  read.table(file_path, header = T, sep = "\t", stringsAsFactors=F)}

# Function to create the plot
create_plot <- function(data, x_var, y_var, output_file, x_label, y_label) {
  jj <- ggplot(data = data, aes(x = Minimap2 +1, y = PVG+1,color=Library_type)) +
    geom_point(alpha = 0.5, size = 2) +
    labs(x = x_label, y = y_label, color = "Library type", shape = "") +
    theme_minimal() +
    geom_abline(slope = 1, intercept = 0, alpha = 0.6, linetype = "dashed", color = "grey") +
    geom_abline(slope = 2, intercept = 0, alpha = 0.6, linetype = "dotdash", color = "grey") +
    geom_abline(slope = 3, intercept = 0, alpha = 0.6, linetype = "longdash", color = "grey") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
    scale_color_manual(values = c("AMPLICON" = "red", "WGS" = "blue", "SIMULATED" = "black",
                            "Metagenomic" = "green", "ONT" = "purple")) +
       scale_shape_manual(values = c("AMPLICON"=19, "WGS"=17, "Metagenomic"=15, "ONT"=13, "SIMULATED"=11)) +
    scale_x_log10(limits = c(1, 2500)) + scale_y_log10(limits = c(1, 2500)) #    +
       guides(mcolor = guide_legend(override.aes = list(shape = c(19, 15, 13, 17, 11))),
         shape = guide_legend(override.aes = list(color = c("red", "green", "purple", "red", "black")))  )
 return(jj)
}

# Function to calculate summary statistics
calculate_summary_stats <- function(data) {
  data %>%
    group_by(Library_type) %>%
    summarize(
      PVG_mean = mean(PVG),
      PVG_median = median(PVG),
      PVG_sd = sd(PVG),
      Minimap2_mean = mean(Minimap2),
      Minimap2_median = median(Minimap2),
      Minimap2_sd = sd(Minimap2),
      Total_SNPs_mean = mean(Total),
      Total_SNPs_median = median(Total),
      Total_SNPs_sd = sd(Total) ) }

# Read the data
rates <- read_data("rates.txt")
vgmap6 <- read_data("rates.new3.vgmap6.txt")
vgmap3 <- read_data("rates.new3.vgmap3.txt")
vgmap1 <- read_data("rates.new3.vgmap1.txt")
vgmapAll <- read_data("rates.new3.vgmapall.txt")
giraffe6 <- read_data("rates.new3.giraffe6.txt")
giraffe3 <- read_data("rates.new3.giraffe3.txt")
giraffe1 <- read_data("rates.new3.giraffe1.txt")
giraffeAll <- read_data("rates.new3.giraffeall.txt")

# Generate plots
pdf("finalise_SNPs.pdf", width = 6.5, height = 5)
print(create_plot(rates, "Minimap2", "PVG", "finalise_SNPs.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - PVG"))
  dev.off()
pdf("finalise_vgmap6.pdf", width = 6.5, height = 5)
print(create_plot(vgmap6, "Minimap2", "PVG", "finalise_vgmap6.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_6"))
  dev.off()
pdf("finalise_vgmap3.pdf", width = 6.5, height = 5)
print(create_plot(vgmap3, "Minimap2", "PVG", "finalise_vgmap3.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_3"))
  dev.off()
pdf("finalise_vgmap1.pdf", width = 6.5, height = 5)
print(create_plot(vgmap1, "Minimap2", "PVG", "finalise_vgmap1.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_1"))
  dev.off()
pdf("finalise_vgmapAll.pdf", width = 6.5, height = 5)
  print(create_plot(vgmapAll, "Minimap2", "PVG", "finalise_vgmapAll.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_All"))
  dev.off()
pdf("finalise_giraffe6.pdf", width = 6.5, height = 5)
 print(create_plot(giraffe6, "Minimap2", "PVG", "finalise_giraffe6.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_6"))
  dev.off()
pdf("finalise_giraffe3.pdf", width = 6.5, height = 5)
  print(create_plot(giraffe3, "Minimap2", "PVG", "finalise_giraffe3.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_3"))
  dev.off()
pdf("finalise_giraffe1.pdf", width = 6.5, height = 5)
  print(create_plot(giraffe1, "Minimap2", "PVG", "finalise_giraffe1.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_1"))
  dev.off()
pdf("finalise_giraffeAll.pdf", width = 6.5, height = 5)
  print(create_plot(giraffeAll, "Minimap2", "PVG", "finalise_giraffeAll.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_All"))
  dev.off()
  
# Calculate summary statistics
#summary_stats_rates <- calculate_summary_stats(rates)
summary_stats_vgmap6 <- calculate_summary_stats(vgmap6)
summary_stats_vgmap3 <- calculate_summary_stats(vgmap3)
summary_stats_vgmap1 <- calculate_summary_stats(vgmap1)
summary_stats_vgmapAll <- calculate_summary_stats(vgmapAll)
summary_stats_giraffeAll <- calculate_summary_stats(giraffeAll)
summary_stats_giraffe6 <- calculate_summary_stats(giraffe6)
summary_stats_giraffe3 <- calculate_summary_stats(giraffe3)
summary_stats_giraffe1 <- calculate_summary_stats(giraffe1)

# Print summary statistics
print("All SNPs"); print(summary_stats_rates)
print("All VG-MAP_6 BCF"); print(summary_stats_vgmap6)
print("All VG-MAP_3 BCF"); print(summary_stats_vgmap3)
print("All VG-MAP_All BCF"); print(summary_stats_vgmapAll)
print("All VG-MAP_1 BCF"); print(summary_stats_vgmap1)
print("All Giraffe-1 BCF"); print(summary_stats_giraffe1)
print("All Giraffe-3 BCF"); print(summary_stats_giraffe3)
print("All Giraffe-6 BCF"); print(summary_stats_giraffe6)
print("All Giraffe-All BCF"); print(summary_stats_giraffeAll)

print("Doing FB next")
# FB

rates <- read_data("rates.fb.txt")
vgmap6 <- read_data("rates.fb.vgmap6.txt")
vgmap3 <- read_data("rates.fb.vgmap3.txt")
vgmap1 <- read_data("rates.fb.vgmap1.txt")
vgmapAll <- read_data("rates.fb.vgmapall.txt")
giraffe6 <- read_data("rates.fb.giraffe6.txt")
giraffe3 <- read_data("rates.fb.giraffe3.txt")
giraffe1 <- read_data("rates.fb.giraffe1.txt")
giraffeAll <- read_data("rates.fb.giraffeall.txt")

# Generate plots
pdf("finalise_SNPs.FB.pdf", width = 6.5, height = 5)
print(create_plot(rates, "Minimap2", "PVG", "finalise_SNPs.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - PVG"))
  dev.off()
pdf("finalise_vgmap6.FB.pdf", width = 6.5, height = 5)
print(create_plot(vgmap6, "Minimap2", "PVG", "finalise_vgmap6.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_6"))
  dev.off()
pdf("finalise_vgmap3.FB.pdf", width = 6.5, height = 5)
print(create_plot(vgmap3, "Minimap2", "PVG", "finalise_vgmap3.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_3"))
  dev.off()
pdf("finalise_vgmap1.FB.pdf", width = 6.5, height = 5)
print(create_plot(vgmap1, "Minimap2", "PVG", "finalise_vgmap1.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_1"))
  dev.off()
pdf("finalise_vgmapAll.FB.pdf", width = 6.5, height = 5)
  print(create_plot(vgmapAll, "Minimap2", "PVG", "finalise_vgmapAll.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_All"))
  dev.off()
pdf("finalise_giraffe6.FB.pdf", width = 6.5, height = 5)
 print(create_plot(giraffe6, "Minimap2", "PVG", "finalise_giraffe6.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_6"))
  dev.off()
pdf("finalise_giraffe3.FB.pdf", width = 6.5, height = 5)
  print(create_plot(giraffe3, "Minimap2", "PVG", "finalise_giraffe3.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_3"))
  dev.off()
pdf("finalise_giraffe1.FB.pdf", width = 6.5, height = 5)
  print(create_plot(giraffe1, "Minimap2", "PVG", "finalise_giraffe1.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_1"))
  dev.off()
pdf("finalise_giraffeAll.FB.pdf", width = 6.5, height = 5)
  print(create_plot(giraffeAll, "Minimap2", "PVG", "finalise_giraffeAll.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_All"))
  dev.off()

# Calculate summary statistics
summary_stats_rates <- calculate_summary_stats(rates)
summary_stats_vgmap6 <- calculate_summary_stats(vgmap6)
summary_stats_vgmap3 <- calculate_summary_stats(vgmap3)
summary_stats_vgmap1 <- calculate_summary_stats(vgmap1)
summary_stats_vgmapAll <- calculate_summary_stats(vgmapAll)
summary_stats_giraffeAll <- calculate_summary_stats(giraffeAll)
summary_stats_giraffe6 <- calculate_summary_stats(giraffe6)
summary_stats_giraffe3 <- calculate_summary_stats(giraffe3)
summary_stats_giraffe1 <- calculate_summary_stats(giraffe1)

# Print summary statistics
print("All SNPs FB"); print(summary_stats_rates)
print("All VG-MAP_6 FB"); print(summary_stats_vgmap6)
print("All VG-MAP_3 FB"); print(summary_stats_vgmap3)
print("All VG-MAP_All FB"); print(summary_stats_vgmapAll)
print("All VG-MAP_1 FB"); print(summary_stats_vgmap1)
print("All Giraffe-1 FB"); print(summary_stats_giraffe1)
print("All Giraffe-3 FB"); print(summary_stats_giraffe3)
print("All Giraffe-6 FB"); print(summary_stats_giraffe6)
print("All Giraffe-All FB"); print(summary_stats_giraffeAll)

print("Doing VG next")

# VG

vgmap6 <- read_data("rates.vg.vgmap6.txt")
vgmap3 <- read_data("rates.vg.vgmap3.txt") 
vgmap1 <- read_data("rates.vg.vgmap1.txt") 
vgmapAll <- read_data("rates.vg.vgmapall.txt")
giraffe6 <- read_data("rates.vg.giraffe6.txt")
giraffe3 <- read_data("rates.vg.giraffe3.txt")
giraffe1 <- read_data("rates.vg.giraffe1.txt")
giraffeAll <- read_data("rates.vg.giraffeall.txt")

print("VG VGMAP-6")
print(calculate_summary_stats(vgmap6))
print("VG VGMAP-3")
print(calculate_summary_stats(vgmap3))
print("VG VGMAP-1")
print(calculate_summary_stats(vgmap1))
print("VG VGMAP-All")
print(calculate_summary_stats(vgmapAll))
print("VG Giraffe-6")
print(calculate_summary_stats(giraffe6))
print("VG Giraffe-3")
print(calculate_summary_stats(giraffe3))
print("VG Giraffe-1")
print(calculate_summary_stats(giraffe1))
print("VG Giraffe-All")
print(calculate_summary_stats(giraffeAll))

# Generate plots

pdf("finalise_vgmap6.VG.pdf", width = 6.5, height = 5)
print(create_plot(vgmap6, "Minimap2", "PVG", "finalise_vgmap6.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_6"))
  dev.off()
pdf("finalise_vgmap3.VG.pdf", width = 6.5, height = 5)
print(create_plot(vgmap3, "Minimap2", "PVG", "finalise_vgmap3.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_3"))
  dev.off()
pdf("finalise_vgmap1.VG.pdf", width = 6.5, height = 5)
print(create_plot(vgmap1, "Minimap2", "PVG", "finalise_vgmap1.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_1"))
  dev.off()
pdf("finalise_vgmapAll.VG.pdf", width = 6.5, height = 5)
  print(create_plot(vgmapAll, "Minimap2", "PVG", "finalise_vgmapAll.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - VG-MAP_All"))
  dev.off() 
pdf("finalise_giraffe6.VG.pdf", width = 6.5, height = 5)
 print(create_plot(giraffe6, "Minimap2", "PVG", "finalise_giraffe6.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_6"))
  dev.off()
pdf("finalise_giraffe3.VG.pdf", width = 6.5, height = 5)
  print(create_plot(giraffe3, "Minimap2", "PVG", "finalise_giraffe3.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_3"))
  dev.off()
pdf("finalise_giraffe1.VG.pdf", width = 6.5, height = 5)
  print(create_plot(giraffe1, "Minimap2", "PVG", "finalise_giraffe1.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_1"))
  dev.off()
pdf("finalise_giraffeAll.VG.pdf", width = 6.5, height = 5)
  print(create_plot(giraffeAll, "Minimap2", "PVG", "finalise_giraffeAll.pdf",
  "Number of SNPs - Minimap2", "Number of SNPs - Giraffe_All"))
  dev.off()