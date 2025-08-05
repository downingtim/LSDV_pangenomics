# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# Read the CSV data
data <- read_csv("get_liftover.csv")
colnames(data) =  gsub("_GB", "GB", colnames(data))

# Convert the data from wide format to long format
long_data <- data %>%
  pivot_longer(cols = starts_with("3_") | starts_with("6_"), #  | starts_with("ALL_"),
               names_to = c("PVG_Type", "Caller"),
               names_sep = "_",
               values_to = "Fraction_Not_Transferred")

long_data$PVG_Type <- factor(long_data$PVG_Type, levels = c("3", "3_GBWT", "6", "6_GBWT", "ALL", "ALL_GBWT"))

# Create the ggplot
plot <- ggplot(long_data, aes(x = Sample, y = Fraction_Not_Transferred)) +
  geom_point(size=0.8, alpha = 0.5) +
  facet_grid(Caller ~ PVG_Type, scales = "free_x") +  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(       x = "Sample", y = "Fraction Not Transferred") +
  theme( axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 10))

# Save the plot to a PDF file
ggsave("get_liftover.pdf", plot, width =21, height =6)
