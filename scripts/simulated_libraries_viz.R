library(vcfR)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggpubr)

read_snps <- function(vcf_file) {
  if (!file.exists(vcf_file)) {
    warning("Missing file: ", vcf_file)
    return(tibble(chr = character(), pos = integer()))  }
  vcf <- tryCatch(  read.vcfR(vcf_file, verbose = FALSE),
    error = function(e) NULL )
  if (is.null(vcf) || nrow(vcf@fix) == 0) {
    return(tibble(chr = character(), pos = integer())) }
  as_tibble(vcf@fix) %>%  transmute(chr = CHROM, pos = as.integer(POS))}

compare_sets <- function(truth, calls, outdir, sample) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  truth_sites <- paste(truth$chr, truth$pos)
  call_sites  <- paste(calls$chr, calls$pos)
  TP <- intersect(truth_sites, call_sites)
  FN <- setdiff(truth_sites, call_sites)
  FP <- setdiff(call_sites, truth_sites)
  write.table(     tibble(site = FN),
    file = file.path(outdir, paste0(sample, "_FN.tsv")),
    quote = FALSE, row.names = FALSE, sep = "\t" )
  write.table(     tibble(site = FP),
    file = file.path(outdir, paste0(sample, "_FP.tsv")),
    quote = FALSE, row.names = FALSE, sep = "\t" )
  precision <- ifelse(length(TP) + length(FP) == 0, NA,
                      length(TP) / (length(TP) + length(FP)))
  recall <- ifelse(length(TP) + length(FN) == 0, NA,
                   length(TP) / (length(TP) + length(FN)))
  F1 <- ifelse(is.na(precision) | is.na(recall) | precision + recall == 0,
               NA, 2 * precision * recall / (precision + recall))
  tibble(Both = length(TP)/length(truth_sites),
            FP= length(FP)/length(truth_sites), Missed=length(FN)/length(truth_sites), 
      precision = precision,  recall = recall, F1 = F1) }

vcf_path <- function(method, sample) {
  file.path( method, "BCF_VCF_FILES", paste0(sample, ".norm.bcf.vcf.gz") )}

refs <- c("OQ511520", "KX764645")
#refs <- c("OQ511520")#, "KX764645")
#refs <- c("KX764645")
#print(refs)
depths <- seq(10, 200, 10)
samples <- expand.grid(ref = refs, depth = depths) %>%
  mutate(sample = paste0("sim_", ref, "_", depth))

setwd("/mnt/lustre/RDS-ephemeral/downing/LSDV/SIM_TEST/")
#print(getwd())
results <- map_dfr(seq_len(nrow(samples)), function(i) {

  ref   <- samples$ref[i]
  depth <- samples$depth[i]
  sample <- samples$sample[i]

  lsdv1   <- read_snps(vcf_path("../LSDV1", sample))
  lsdvg   <- read_snps(vcf_path("../LSDVG", sample))
  lsdvg3  <- read_snps(vcf_path("../LSDVG_3", sample))
  lsdvg6  <- read_snps(vcf_path("../LSDVG_6", sample))

  lsdvg_3_m <- bind_rows(lsdvg, lsdvg3) %>% distinct()
  lsdvg_6_m <- bind_rows(lsdvg, lsdvg6) %>% distinct()

  sample1 <- str_split_1(sample, "_")[2]

  # Reconstruct truth path correctly
  truth_file <- paste0("KX894508_", sample1, ".vcf")
  truth <- read_snps(truth_file)

  if(nrow(truth) == 0){
    return(NULL)
  }

  bind_rows(
    compare_sets(truth, lsdv1, "Linear", sample) %>%
      mutate(comparison = "Linear"),

    compare_sets(truth, lsdvg, "Giraffe-PVG1", sample) %>%
      mutate(comparison = "Giraffe-PVG1"),

    compare_sets(truth, lsdvg_3_m, "Giraffe-PVG13", sample) %>%
      mutate(comparison = "Giraffe-PVG13"),

    compare_sets(truth, lsdvg_6_m, "Giraffe-PVG16", sample) %>%
      mutate(comparison = "Giraffe-PVG16")

  ) %>%
    mutate(sample = sample,
           depth = depth,
           ref = ref) })

write.csv( results,  "sim_comparison_results.csv", row.names=F)
head(results)
subset(results, comparison=="Giraffe-PVG1")

tp1 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"
                                    & comparison!="Giraffe-PVG16"),
  cols = c(Both), names_to = "metric", values_to = "value"),
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2) + #facet_grid(~ comparison) +
  scale_x_continuous(limits = c(10, 90),  breaks = seq(10, 90, 10)) +
  scale_y_continuous(limits = c(0.5, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Fraction of true positive SNPs") +
  scale_linetype_manual( values = c("Both" = "solid", "Missed" = "dashed"),
  labels = c("Both" = "True Positive", "Missed" = "False Negative"),
  name = "Metric") + guides(linetype = "none") + guides(colour = "none") +
  scale_color_manual(values = c("Linear" = "red", "Giraffe-PVG13" = "blue",
                                "Giraffe-PVG16" = "cyan")) + 
  geom_smooth(se = FALSE, linewidth =2, method = "loess", alpha=0.4)

tp2 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"
                                    & comparison!="Giraffe-PVG16"),
  cols = c(Both), names_to = "metric", values_to = "value"),
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2) + #facet_grid(~ comparison) +
  scale_x_continuous(limits = c(80, 200),  breaks = seq(80, 200, 10)) +
  scale_y_continuous(limits = c(0.94, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Fraction of true positive SNPs") +
  scale_linetype_manual( values = c("Both" = "solid", "Missed" = "dashed"),
  labels = c("Both" = "True Positive", "Missed" = "False Negative"),
  name = "Metric") + guides(linetype = "none") + guides(colour = "none") +
    geom_smooth(se = FALSE, linewidth =2, method = "loess", alpha=0.4) +
      scale_color_manual(values = c("Linear" = "red", "Giraffe-PVG13" = "blue",
                                "Giraffe-PVG16" = "cyan"))

fn1 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1" 
                                    & comparison!="Giraffe-PVG16"),
  cols = c(Missed), names_to = "metric", values_to = "value"),
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2) + #facet_grid(~ comparison) +
  scale_x_continuous(limits = c(10, 90),  breaks = seq(10, 90, 10)) +
  scale_y_continuous(limits = c(0, 0.5)) + theme_bw() +
  labs(x = "Read depth", y = "Fraction of false negative SNPs") +
  scale_linetype_manual( values = c("Both" = "solid", "Missed" = "dashed"),
  labels = c("Both" = "True Positive", "Missed" = "False Negative"),
  name = "Metric") + guides(linetype = "none") + guides(colour = "none") +
  scale_color_manual(values = c("Linear" = "red", "Giraffe-PVG13" = "blue",
                                "Giraffe-PVG16" = "cyan")) + 
  geom_smooth(se = FALSE, linewidth = 2, method = "loess", alpha=0.4)

fn2 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"
                                    & comparison!="Giraffe-PVG16"),
  cols = c(Missed), names_to = "metric", values_to = "value"),
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2) + #facet_grid(~ comparison) +
  scale_x_continuous(limits = c(80, 200),  breaks = seq(80, 200, 10)) +
  scale_y_continuous(limits = c(0, 0.06)) + theme_bw() +
  labs(x = "Read depth", y = "Fraction of false negative SNPs") +
  scale_linetype_manual( values = c("Both" = "solid", "Missed" = "dashed"),
  labels = c("Both" = "True Positive", "Missed" = "False Negative"),
  name = "Metric") + guides(linetype = "none") + guides(colour = "none") +
  scale_color_manual(values = c("Linear" = "red", "Giraffe-PVG13" = "blue",
                                "Giraffe-PVG16" = "cyan")) + 
  geom_smooth(se = FALSE, linewidth = 2, method = "loess", alpha=0.4)

fp1 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(FP), names_to = "metric", values_to = "value"), 
  aes(depth, 100*value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2) + 
  geom_smooth(se = FALSE, linewidth = 2, method = "loess", alpha=0.4) +
   scale_x_continuous(limits = c(0, 200),  breaks = seq(10, 200, 10)) +
   scale_y_continuous(limits = c(0, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Percentage of false positive SNPs") +
  guides(linetype = "none") +
  scale_color_manual(values = c("Linear" = "red", "Giraffe-PVG13" = "blue",
                                "Giraffe-PVG16" = "cyan")) 
pdf("FP_rate.pdf", width =7, height =4)
print(fp1)
dev.off()

tpfp <- ggarrange(tp1, tp2, fn1, fn2, nrow = 2,  ncol = 2, widths = c(1, 1.5), heights = c(1, 1),
          labels = c("A", "B", "C", "D"),  align = "v")         # Optional: vertically aligns axes
ggsave("rates.pdf", tpfp, width =6, height =6)

# plot above for recall and precision and F1 score

recall1 <- ggplot(pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(recall), names_to = "metric", values_to = "value"), 
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2,  alpha=0.3) + 
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.3) +
   scale_x_continuous(limits = c(0, 90),  breaks = seq(10, 90, 20)) +
   scale_y_continuous(limits = c(0.45, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Recall") +
  guides(linetype = "none") +guides(colour = "none") +
  scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

recall2 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(recall), names_to = "metric", values_to = "value"), 
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2,  alpha=0.3) + 
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.4) +
   scale_x_continuous(limits = c(80, 200),  breaks = seq(80, 200, 20)) +
   scale_y_continuous(limits = c(0.925, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Recall") +
  guides(linetype = "none") +guides(colour = "none") +
  scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

precision1 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(precision), names_to = "metric", values_to = "value"), 
  aes(depth, value,  colour = comparison, linetype = metric)) + 
  geom_point(size = 2,  alpha=0.3) +
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.3) +
   scale_x_continuous(limits = c(0, 90),  breaks = seq(10, 90, 20)) +
   scale_y_continuous(limits = c(0.97, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Precision") +
  guides(linetype = "none") +guides(colour = "none") +
  scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

precision2 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(precision), names_to = "metric", values_to = "value"), 
  aes(depth, value,  colour = comparison, linetype = metric)) +       
  geom_point(size = 2,  alpha=0.3) +
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.3) +
   scale_x_continuous(limits = c(80, 200),  breaks = seq(80, 200, 20)) +
   scale_y_continuous(limits = c(0.982, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Precision") +
  guides(linetype = "none") +guides(colour = "none") +
   scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

prec_recall <-  ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(precision), names_to = "metric", values_to = "value"),
  aes(recall, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2,  alpha=0.3) +
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.3) +
#   scale_x_continuous(limits = c(0, 90),  breaks = seq(10, 90, 20)) +
#   scale_y_continuous(limits = c(0.97, 1)) +
  theme_bw() + labs(x = "Recall", y = "Precision") +
  guides(linetype = "none") +guides(colour = "none") +
  scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

F1_1 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(F1), names_to = "metric", values_to = "value"), 
  aes(depth, value,  colour = comparison, linetype = metric)) + 
  geom_point(size = 2,  alpha=0.3) +
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.3) +
   scale_x_continuous(limits = c(0, 90),  breaks = seq(10, 90, 20)) +
   scale_y_continuous(limits = c(0.6, 1)) + theme_bw() +
  labs(x = "Read depth", y = "F1 score") +
  guides(linetype = "none") +guides(colour = "none") +
 scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

F1_2 <- ggplot( pivot_longer( subset(results, comparison!="Giraffe-PVG1"),
  cols = c(F1), names_to = "metric", values_to = "value"),
  aes(depth, value,  colour = comparison, linetype = metric)) +
  geom_point(size = 2,  alpha=0.3) +
  geom_smooth(se = FALSE, linewidth = 1.5, method = "loess", alpha=0.3) +
   scale_x_continuous(limits = c(80, 200),  breaks = seq(80, 200, 20)) +
   scale_y_continuous(limits = c(0.96, 1)) + theme_bw() +
  labs(x = "Read depth", y = "F1 score") +
  guides(linetype = "none") +guides(colour = "none") +
  scale_color_manual(values = c( "Linear" = scales::alpha("red", 0.3),
  "Giraffe-PVG13" = scales::alpha("blue", 0.3), "Giraffe-PVG16" = scales::alpha("cyan", 0.3) ))

scores <- ggarrange(recall1, recall2, precision1, precision2, F1_1, F1_2,
           nrow = 3,  ncol = 2, widths = c(1, 1.3), heights = c(1, 1, 1),
          labels = c("A", "B", "C", "D", "E", "F"),  align = "v")
ggsave("rates_scores.pdf", scores, width =6, height =6)

ggsave("prec_recall.pdf", prec_recall, width =6, height =6)

#results[93,6] = .98 # fix weird value
pdf("F1_depth_PVG1.pdf", width =8, height =4)
ggplot( pivot_longer(
  subset(results, comparison=="Linear" | comparison=="Giraffe-PVG1"),
  cols = c(F1), names_to = "metric", values_to = "value" ), 
  aes(depth, value,  colour = comparison, alpha=0.6, linetype = metric)) +
  geom_line(linewidth = 1) + geom_point(size = 2) + facet_grid(~ref) +
  scale_x_continuous(limits = c(10, 200),  breaks = seq(10, 200, 10)) +
  scale_y_continuous(limits = c(0.6, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Fraction of valid SNPs")
dev.off()

pdf("F1_depth_PVG_zoom2.pdf", width =8, height =4)
ggplot( pivot_longer(
  subset(results, comparison=="Linear" | comparison=="Giraffe-PVG1"),
  cols = c(F1), names_to = "metric", values_to = "value" ), 
  aes(depth, value,  colour = comparison, linetype = metric)) +
      geom_smooth(se = FALSE, linewidth = 1.2, alpha=0.6, method = "loess") +
  geom_point(size = 2) + facet_grid(~ref) +
  scale_x_continuous(limits = c(50, 200),  breaks = seq(10, 200, 20)) +
  scale_y_continuous(limits = c(0.961, 1)) + theme_bw() +
  labs(x = "Read depth", y = "Fraction of valid SNPs")
dev.off()

# ---------- SNP detection heatmap ----------

make_detection_plot <- function(ref_name, giraffe_label, giraffe_path) {
  depths_rev <- sort(depths, decreasing = TRUE)
  truth_sample <- paste0("sim_", ref_name, "_200")
  truth <- bind_rows(
    read_snps(vcf_path("../LSDV1", truth_sample)),
    read_snps(vcf_path("../LSDVG", truth_sample)),
    read_snps(vcf_path("../LSDVG_3", truth_sample)),
    read_snps(vcf_path("../LSDVG_6", truth_sample))
  ) %>% distinct()
  truth$site <- paste(  truth$pos)

  plot_df <- map_dfr(depths_rev, function(d) {
    sample <- paste0("sim_", ref_name, "_", d)
    linear  <- read_snps(vcf_path("../LSDV1", sample))
    giraffe <- read_snps(vcf_path(giraffe_path, sample))
    linear_sites <- linear$pos
    if (giraffe_label == "Giraffe-PVG1") {
      giraffe_sites <- giraffe$pos     } else {
    giraffe_base <- read_snps(vcf_path("../LSDVG", sample))
    giraffe_sites <- unique(c(giraffe_base$pos, giraffe$pos)) }

    tibble( depth = d,site  = truth$site,
      status = case_when(
        site %in% linear_sites & site %in% giraffe_sites ~ "Both",
        site %in% linear_sites & !site %in% giraffe_sites ~ "Linear",
        !site %in% linear_sites & site %in% giraffe_sites ~ "Giraffe",
        TRUE ~ "Neither" ) )})
  plot_df$depth <- factor(plot_df$depth, levels = depths_rev)
  plot_df$comparison <- paste(ref_name, "Linear vs", giraffe_label)
  plot_df}

heatmap_df <- bind_rows(
  make_detection_plot("KX764645", "Giraffe-PVG1", "../LSDVG"),
  make_detection_plot("KX764645", "Giraffe-PVG13", "../LSDVG_3"),
  make_detection_plot("KX764645", "Giraffe-PVG16", "../LSDVG_6"))
heatmap_df$comparison <- factor(   heatmap_df$comparison,
  levels = unique(heatmap_df$comparison) )
heatmap_df <- heatmap_df %>%
  mutate(site_num = as.numeric(as.character(site)))
ordered_sites <- heatmap_df %>%
  distinct(site, site_num) %>% arrange(site_num) %>% pull(site)
heatmap_df$site <- factor(heatmap_df$site, levels = ordered_sites)

heatmap_df2 <- bind_rows( 
  make_detection_plot("OQ511520", "Giraffe-PVG1", "../LSDVG"),
  make_detection_plot("OQ511520", "Giraffe-PVG13", "../LSDVG_3"),
  make_detection_plot("OQ511520", "Giraffe-PVG16", "../LSDVG_6") )
heatmap_df2$comparison <- factor(   heatmap_df2$comparison,
  levels = unique(heatmap_df2$comparison) )
heatmap_df2 <- heatmap_df2 %>%
  mutate(site_num = as.numeric(as.character(site)))
ordered_sites <- heatmap_df2 %>%
  distinct(site, site_num) %>% arrange(site_num) %>% pull(site)
heatmap_df2$site <- factor(heatmap_df2$site, levels = ordered_sites)

p <- ggplot(heatmap_df, aes(x = site, y = depth, fill = status)) +
  geom_tile() + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" )) +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by = 30)]) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")
ggsave("SNP_detection_heatmap_KX764645.pdf", p, width = 16, height = 8)

p2 <- ggplot(heatmap_df2, aes(x = site, y = depth, fill = status)) +
  geom_tile() + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" )) +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by = 18)]) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")
ggsave("SNP_detection_heatmap_OQ511520.pdf", p2, width = 16, height = 8)

pboth <- ggarrange(p, p2, nrow = 2,  ncol = 1,
          labels = c("A", "B"),  align = "v")         # Optional: vertically aligns axes
ggsave("SNP_detection_heatmap_both.pdf", pboth, width = 16, height = 16)

## zoom in

heatmap_df_z <- bind_rows(
   make_detection_plot("KX764645", "Giraffe-PVG13", "../LSDVG_3"))
heatmap_df_z$comparison <- factor(   heatmap_df_z$comparison,
  levels = unique(heatmap_df_z$comparison) )
heatmap_df_z <- heatmap_df_z %>%
  mutate(site_num = as.numeric(as.character(site)))
ordered_sites <- heatmap_df_z %>%
  distinct(site, site_num) %>% arrange(site_num) %>% pull(site)
heatmap_df_z$site <- factor(heatmap_df_z$site, levels = ordered_sites)

heatmap_df_z2 <- bind_rows(
   make_detection_plot("OQ511520", "Giraffe-PVG13", "../LSDVG_3"))
heatmap_df_z2$comparison <- factor(   heatmap_df_z2$comparison,
  levels = unique(heatmap_df_z2$comparison) )
heatmap_df_z2 <- heatmap_df_z2 %>%
  mutate(site_num = as.numeric(as.character(site)))
ordered_sites <- heatmap_df_z2 %>%
  distinct(site, site_num) %>% arrange(site_num) %>% pull(site)
heatmap_df_z2$site <- factor(heatmap_df_z2$site, levels = ordered_sites)

plot_z <- subset(heatmap_df_z %>% arrange((site_num)), depth !=10 & site_num < 2490)
p_z <- ggplot( plot_z, aes(x = site, y = depth, fill = status)) +
  geom_tile(show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

plot_z2 <- subset(heatmap_df_z2 %>% arrange((site_num)), depth !=10 & site_num < 2490)
p_z2 <- ggplot( plot_z2, aes(x = site, y = depth, fill = status)) +
  geom_tile( show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

pboth_z <- ggarrange(p_z, p_z2, nrow =1,  ncol = 2, widths = c(5, 4),
          labels = c("A", "B"),  align = "v")         # Optional: vertically aligns axes
ggsave("SNP_detection_heatmap_both_zoom_2500.pdf", pboth_z, width =4.5, height =4)

# check variation at depth 200 for giraffe-PVG13
plot_y <- subset(heatmap_df_z %>% arrange((site_num)), depth==200 & status=="Giraffe")
str(plot_y)
data.frame(plot_y)

# plot around 5kb region with lots of variation
plot_z5 <- subset(heatmap_df_z %>% arrange((site_num)), 
          depth !=10 & site_num >5500 & site_num < 5600)
p_z5 <- ggplot( plot_z5, aes(x = site, y = depth, fill = status)) +
  geom_tile(show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

# plot around 5kb region with lots of variation
plot_z2_5 <- subset(heatmap_df_z2 %>% arrange((site_num)), 
           depth !=10 &  site_num >5500 & site_num < 5600)
p_z2_5 <- ggplot( plot_z2_5, aes(x = site, y = depth, fill = status)) +
  geom_tile( show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

# combine the two zoomed plots
pboth_z_5 <- ggarrange(p_z5, p_z2_5, nrow =1,  ncol = 2, widths = c(5, 5),
          labels = c("A", "B"),  align = "v")         # Optional: vertically aligns axes
ggsave("SNP_detection_heatmap_both_zoom_5kb.pdf", pboth_z_5, width =5, height =5)

# plot around 58kb region with lots of variation
plot_z58 <- subset(heatmap_df_z %>% arrange((site_num)),
       depth != 10 & site_num >57800 & site_num < 57999)
p_z58 <- ggplot( plot_z58, aes(x = site, y = depth, fill = status)) +
  geom_tile(show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

# plot around 58kb region with lots of variation
plot_z2_58 <- subset(heatmap_df_z2 %>% arrange((site_num)), 
           depth != 10 &  site_num >57800 & site_num < 57999)
p_z2_58 <- ggplot( plot_z2_58, aes(x = site, y = depth, fill = status)) +
  geom_tile( show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

# combine the two zoomed plots
pboth_z58 <- ggarrange(p_z58, p_z2_58, nrow =1,  ncol = 2, widths = c(5,4),
          labels = c("A", "B"),  align = "v")         # Optional: vertically aligns axes
ggsave("SNP_detection_heatmap_both_zoom_58kb.pdf", pboth_z58, width =5, height =5)

# plot around 148900 region with lots of variation
plot_z148 <- subset(heatmap_df_z %>% arrange((site_num)),  depth != 10 & site_num >148700)
p_z148 <- ggplot( plot_z148, aes(x = site, y = depth, fill = status)) +
  geom_tile(show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

# plot around 148900 region with lots of variation
plot_z2_148 <- subset(heatmap_df_z2 %>% arrange((site_num)), depth != 10 & site_num >148700)
p_z2_148 <- ggplot( plot_z2_148, aes(x = site, y = depth, fill = status)) +
  geom_tile( show.legend = F) + facet_grid(comparison ~ .) +
  scale_fill_manual(values = c(
    "Linear" = "red",  "Giraffe" = "cyan",
    "Both" = "black",   "Neither" = "white" ), guide = "none") +
  scale_x_discrete(  breaks = unique(heatmap_df$site)[
      seq(1, length(unique(heatmap_df$site)), by =1)]) +
  theme_bw() + theme(axis.text.x = element_text(angle =90, hjust = 1),
   legend.position = "none") +
  labs(x = "SNP position", y = "Read depth", fill = "Detection")

# combine the two zoomed plots
pboth_z148 <- ggarrange(p_z148, p_z2_148, nrow =1,  ncol = 2, widths = c(5, 3),
          labels = c("A", "B"),  align = "v")         # Optional: vertically aligns axes
ggsave("SNP_detection_heatmap_both_zoom_148kb.pdf", pboth_z148, width =5, height =5) 
