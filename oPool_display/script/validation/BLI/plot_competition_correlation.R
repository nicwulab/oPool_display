#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)

# Ensure output directory exists
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Load competition data
competition_index_data <- read_tsv("oPool_result/filtered_hits/validated_antibody_competition_index.tsv", show_col_types = FALSE)
competition_percentage_data <- read_tsv("experimental_data/validation/BLI/validated_antibody_competition_percentage.tsv", show_col_types = FALSE)

# Merge tables by antigen and antibody
merged_data <- inner_join(competition_index_data, competition_percentage_data, by = c("Antigen", "Antibody"))

# Define antigen color map (ensuring it matches your previous aesthetics)
color_map <- c(
  "SI06_H1" = "#D18FE7",
  "MI15_H1" = "#4CA998",
  "SP16_H3" = "#E47AB5",
  "QH_H5"   = "#4A9FD8",
  "SH_H7"   = "#5B72D4",
  "Phu_FluB" = "#B08F3E",
  "Lee_FluB" = "#C8645D"
)

# Ensure the output directory exists
ensure_dir("graph")

# Create scatter plot with aesthetics strictly following previous scripts
p <- ggplot(merged_data, aes(x = Competition_Percentage, y = Competition_Index, color = Antigen)) +
  geom_point(size = 2, shape = 16, stroke = 0.8) +  # Point size and shape
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, linetype = "solid", color = "red") +  # Regression line
  scale_color_manual(values = color_map) +  # Custom colors
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_y_continuous(breaks = seq(-2, 1.5, by = 0.5))+
  theme_cowplot(12) +
  theme(
    plot.background = element_rect(fill = "white"),
    axis.text = element_text(size = 7, face = "bold", colour = 'black'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = 'black'),
    axis.title = element_text(size = 7, face = "bold"),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "none",
    plot.margin = margin(5, 5, 10, 5)
  ) +
  labs(
    x = bquote(bold("Competition Percentage (%)")),
    y = bquote(bold("Competition Index"))
    )
  print(paste("Correlation :", cor(merged_data$Competition_Percentage, merged_data$Competition_Index, method = "pearson")))

# Save the plot
output_file <- "graph/validation/BLI/competition_index_vs_percentage.png"
ggsave(output_file, plot = p, width = 2.3, height = 2.3, dpi = 600)

print(paste("Plot saved to:", output_file))