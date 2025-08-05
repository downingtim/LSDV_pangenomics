library(ggplot2)
library(patchwork)

coverage_dir <- '/mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVG_3/COVERAGE'
rates <- read.table('rates.new3.txt', header=T, sep='\t')
rates <- rates[, c("Sample", "PVG", "M", "all")]
rates <- rates[!grepl("paired", rates$Sample), ]
rates <- rates[!grepl("SRR1909074", rates$Sample), ]
rates <- rates[!grepl("SRR2334", rates$Sample), ]
rates <- rates[!grepl("SRR2559", rates$Sample), ]

metadata <- read.table('metadata2.csv', sep='\t', header=FALSE,
    col.names = c("Sample", "Library_type"))
rates <- merge(rates, metadata, by = "Sample")
coverage_data <- matrix(NA, nrow = nrow(rates), ncol = 3)
colnames(coverage_data) <- c("KX894508_meandepth", "OQ511520_meandepth", "KX764645_meandepth")

# Iterate over each sample and extract the meandepth values  
for (i in 1:nrow(rates)) {
  sample <- rates$Sample[i]
  file_path <- file.path(coverage_dir, paste0(sample, ".coverage.txt"))
  print(file_path)
  if ( (file.exists(file_path) && file.info(file_path)$size > 0) ) {
    coverage_file <- read.table(file_path, sep = '\t', header = F)
    meandepth_KX894508 <- coverage_file$V7[coverage_file$V1 == "KX894508"]
    meandepth_OQ511520 <- coverage_file$V7[coverage_file$V1 == "OQ511520"]
    meandepth_KX764645 <- coverage_file$V7[coverage_file$V1 == "KX764645"]
    coverage_data[i, ] <- c(meandepth_KX894508, meandepth_OQ511520,
       meandepth_KX764645) } else {
    coverage_data[i, ] <- c(NA, NA, NA)  } }

final_table <- cbind(rates, coverage_data)
meandepth_sums <- rowSums(final_table[, c("KX894508_meandepth",
   "OQ511520_meandepth", "KX764645_meandepth")], na.rm = TRUE)
final_table$PVG_unique = as.numeric(final_table$all-final_table$M)
final_table$M_unique = as.numeric(final_table$all-final_table$PVG)
final_table$KX894508_normalized <- final_table$KX894508_meandepth / meandepth_sums
final_table$OQ511520_normalized <- final_table$OQ511520_meandepth / meandepth_sums
final_table$KX764645_normalized <- final_table$KX764645_meandepth / meandepth_sums

add_correlation <- function(x, y) {
  corr_value <- round(cor(x, y, use = "complete.obs"), 2)
  return(paste("r=", sprintf("%.2f", corr_value), sep="")) }

# First row plots (PVG)
p1P <- ggplot(final_table, aes(x=KX894508_normalized, y=PVG_unique,
    color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=T, color='black') +
  annotate("text", x=0.95, y=max(final_table$PVG_unique), 
     label=add_correlation(final_table$KX894508_normalized, final_table$PVG_unique)) +
  theme_minimal() +  scale_x_continuous(limits = c(0, 1)) + 
  xlab("KX894508 fraction of reads") +
  ylab("Number of SNPs unique to PVG-based approaches detected") +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))
        
p1M <- ggplot(final_table, aes(x=KX894508_normalized, y=M_unique,
    color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=T, color='black') +
  annotate("text", x=0.95, y=max(final_table$M_unique), 
    label=add_correlation(final_table$KX894508_normalized, final_table$M_unique)) +
  theme_minimal() +  scale_x_continuous(limits = c(0, 1)) + 
  xlab("KX894508 fraction of reads") +
  ylab("Number of SNPs unique to Minimap2 detected") +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))
        
p2P <- ggplot(final_table, aes(x=OQ511520_normalized, y=PVG_unique,
    color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=T, color='black') +
  annotate("text", x=0.95, y=max(final_table$PVG_unique), 
     label=add_correlation(final_table$OQ511520_normalized, final_table$PVG_unique)) +
  theme_minimal() +  scale_x_continuous(limits = c(0, 1)) + 
  xlab("OQ511520 fraction of reads") +
  ylab("Number of SNPs unique to PVG-based approaches detected") +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))
        
p2M <- ggplot(final_table, aes(x=OQ511520_normalized, y=M_unique,
    color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=T, color='black') +
  annotate("text", x=0.95, y=max(final_table$M_unique), 
    label=add_correlation(final_table$OQ511520_normalized, final_table$M_unique)) +
  theme_minimal() +  scale_x_continuous(limits = c(0, 1)) + 
  xlab("OQ511520 fraction of reads") +
  ylab("Number of SNPs unique to Minimap2 detected") +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))
        
p1 <- ggplot(final_table, aes(x=KX894508_normalized, y=PVG, color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  annotate("text", x=0.25, y=max(final_table$PVG),
   label=add_correlation(final_table$KX894508_normalized, final_table$PVG)) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_log10() +
  xlab("KX894508 fraction of reads") +
  ylab("Number of PVG-based SNPs detected") +
  theme(legend.position = c(0.05, 0.5),
        legend.justification = c(0, 1))
        
p2 <- ggplot(final_table, aes(x=OQ511520_normalized, y=PVG, color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  annotate("text", x=0.05, y=max(final_table$PVG),
           label=add_correlation(final_table$OQ511520_normalized, final_table$PVG)) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_log10() +
  xlab("OQ511520 fraction of reads") +
  ylab("Number of PVG-based SNPs detected") +
  theme(legend.position = c(0.55, 0.5),
        legend.justification = c(0, 1))

p3 <- ggplot(final_table, aes(x=KX764645_normalized, y=PVG, color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  annotate("text", x=0.05, y=max(final_table$PVG),
           label=add_correlation(final_table$KX764645_normalized, final_table$PVG)) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_log10() +
  xlab("KX764645 fraction of reads") +
  ylab("Number of PVG-based SNPs detected") +
  theme(legend.position = c(0.55, 0.5),
        legend.justification = c(0, 1))

# Second row plots (Minimap2)
p4 <- ggplot(final_table, aes(x=KX894508_normalized, y=M, color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  annotate("text", x=0.65, y=max(final_table$M),
           label=add_correlation(final_table$KX894508_normalized, final_table$M)) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_log10() +
  xlab("KX894508 fraction of reads") +
  ylab("Number of Minimap2-based SNPs detected") +
  theme(legend.position = c(0.45, 0.5),
        legend.justification = c(1, 1))

p5 <- ggplot(final_table, aes(x=OQ511520_normalized, y=M, color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  annotate("text", x=0.15, y=max(final_table$M),
           label=add_correlation(final_table$OQ511520_normalized, final_table$M)) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_log10() +
  xlab("OQ511520 fraction of reads") +
  ylab("Number of Minimap2-based SNPs detected") +
  theme(legend.position = c(0.5, 0.5),
        legend.justification = c(0, 1))

p6 <- ggplot(final_table, aes(x=KX764645_normalized, y=M, color=Library_type)) +
  geom_point(size=4, alpha=0.5) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  annotate("text", x=0.15, y=max(final_table$M),
           label=add_correlation(final_table$KX764645_normalized, final_table$M)) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_log10() +
  xlab("KX764645 fraction of reads") +
  ylab("Number of Minimap2-based SNPs detected") +
  theme(legend.position = c(0.5, 0.5),
        legend.justification = c(0, 1))

# Combine all plots into one multi-panel plot
combined_plot <- (p1 + p2 + p3) / (p4 + p5 + p6) +
  plot_annotation(tag_levels = 'A')

# Print the combined plot
pdf("coverage_v_SNPs.pdf", width=12, height=8)
print(combined_plot)
dev.off()

print ( (p1P + p1M+ p2P + p2M ) + plot_annotation(tag_levels = 'A') )