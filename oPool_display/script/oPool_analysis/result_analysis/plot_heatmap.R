library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(gridExtra)

# Heatmap plotting function for a subset of data
plot_enrichment_heatmap <- function(data, legend_title, color_limits, y_labels, textsize = 5) {
  p <- ggplot(data, aes(x = name, y = Enrichment_Score, fill = parameter)) +
    geom_tile(color = "black", size = 0.1) +
    scale_fill_gradientn(
      colours = c("blue","lightblue", "white", "yellow","yellow", "orange","orange", "red", "red"),
      limits = color_limits,
      breaks = c(0, 25, 50),
      labels = c("<0", "25", ">50"),
      guide = "colorbar",
      na.value = "white"  # Ensure values below the cutoff (NA) are white
    ) +
    coord_fixed(ratio = 1) +
    scale_y_discrete(labels = y_labels) +
    theme_cowplot(12) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.text = element_text(size = textsize, face = "bold", colour = 'black'),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = 'black'),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 5, face = "bold"),
      axis.line = element_blank(),  
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      legend.position = "right"
    ) +
    guides(
      fill = guide_colorbar(
        title.theme = element_text(size = 9, face = "bold", colour = 'black', hjust = 0.5),
        label.theme = element_text(size = 9, face = "bold", colour = 'black'),
        frame.colour = "black",
        frame.linewidth = 1,
        ticks = TRUE,
        ticks.colour = "black",
        barwidth = 0.6, barheight = 7,
        title = legend_title
      )
    ) +
    xlab("") +
    ylab("")
  
  return(p)
}

# Load and preprocess datasets
data1 <- read_tsv("oPool_result/enrichment/Full_HA_enrich.tsv")
data2 <- read_tsv("oPool_result/enrichment/HA_stem_enrich.tsv")
data3 <- read_tsv("oPool_result/enrichment/Full_HA_CR9114_competition_enrich.tsv")

data1 <- data1 %>% rename(name = closest_abs)
data2 <- data2 %>% rename(name = closest_abs)
data3 <- data3 %>% rename(name = closest_abs)

# Combine datasets
combined_data <- full_join(data1, data2, by = "name")
combined_data <- full_join(combined_data, data3, by = "name")

# Identify all columns containing "avg_enrich"
enrich_cols <- grep("avg_enrich", names(combined_data), value = TRUE)

table_s3_1 <- combined_data %>%
  select(name, contains("count"))

write_tsv(table_s3_1, "oPool_result/enrichment/table_s3_1.tsv")

table_s3_2 <- combined_data %>%
  select(name, contains("freq"))

write_tsv(table_s3_2, "oPool_result/enrichment/table_s3_2.tsv")

table_s3_3 <- combined_data %>%
  select(name, contains("enrich"))

write_tsv(table_s3_3, "oPool_result/enrichment/table_s3_3.tsv")

# Normalize each column using Median and IQR
data_normalized <- combined_data %>%
  mutate(across(
    all_of(enrich_cols),
    ~ (.-median(., na.rm = TRUE)) / IQR(., na.rm = TRUE),  # Robust normalization
    .names = "{.col}_normalized"
  ))


# Define antigens and their corresponding column patterns
antigens <- c("SI06_H1", "MI15_H1", "SP16_H3", "QH_H5", "SH_H7", "Phu_FluB", "Lee_FluB")

# Define normalized value cutoffs
normalized_cutoffs <- c(
  "Full_HA_H1_stem_avg_enrich_normalized" = -3,
  "Full_HA_H3_stem_avg_enrich_normalized" = -3,
  "Full_HA_SI06_H1_avg_enrich_normalized" = -4,
  "CR9114_compitition_SI06_H1_avg_enrich_normalized" = -4,
  "Full_HA_MI15_H1_avg_enrich_normalized" = -3,
  "CR9114_compitition_MI15_H1_avg_enrich_normalized" = -2,
  "Full_HA_SP16_H3_avg_enrich_normalized" = -5,
  "Full_HA_QH_H5_avg_enrich_normalized" = -5,
  "CR9114_compitition_QH_H5_avg_enrich_normalized" = -3,
  "Full_HA_SH_H7_avg_enrich_normalized" = -3,
  "CR9114_compitition_SH_H7_avg_enrich_normalized" = -4,
  "Full_HA_Phu_FluB_avg_enrich_normalized" = -3,
  "CR9114_compitition_Phu_FluB_avg_enrich_normalized" = -4,
  "Full_HA_Lee_FluB_avg_enrich_normalized" = -3,
  "CR9114_compitition_Lee_FluB_avg_enrich_normalized" = -4
)

# Apply cutoffs: Replace values below the cutoff with 0, and cap values above 50
data_capped_wide <- data_normalized %>%
  mutate(across(
    all_of(names(normalized_cutoffs)),
    ~ ifelse(. < normalized_cutoffs[cur_column()], 0, ifelse(. > 50, 50, .)),
    .names = "{.col}_capped"
  ))


selected_data <- data_capped_wide %>%
  select(name, contains("_avg_enrich_normalized")) %>%
  select(-contains("capped"))
  
write_tsv(selected_data, "result/table_s3_4.tsv")

# Now proceed with pivot_longer for plotting
data_capped <- data_capped_wide %>%
  pivot_longer(
    cols = ends_with("_capped"),
    names_to = "Enrichment_Score",
    values_to = "parameter"
  ) %>%
  drop_na(all_of(enrich_cols))  # Remove rows where parameter is NA

# Replace Y-axis column names with custom labels
y_labels <- c(
  "Full_HA_H1_stem_avg_enrich_normalized_capped" = "H1 stem", 
  "Full_HA_H3_stem_avg_enrich_normalized_capped" = "H3 stem",
  "Full_HA_SI06_H1_avg_enrich_normalized_capped" = "H1/SI06",
  "CR9114_compitition_SI06_H1_avg_enrich_normalized_capped" = "H1/SI06+CR9114",
  "Full_HA_MI15_H1_avg_enrich_normalized_capped" = "H1/MI15",
  "CR9114_compitition_MI15_H1_avg_enrich_normalized_capped" = "H1/MI15+CR9114",
  "Full_HA_SP16_H3_avg_enrich_normalized_capped" = "H3/SP16",
  "Full_HA_QH_H5_avg_enrich_normalized_capped" = "H5/QH05",
  "CR9114_compitition_QH_H5_avg_enrich_normalized_capped" = "H5/QH05+CR9114",
  "Full_HA_SH_H7_avg_enrich_normalized_capped" = "H7/SH13",
  "CR9114_compitition_SH_H7_avg_enrich_normalized_capped" = "H7/SH13+CR9114",
  "Full_HA_Phu_FluB_avg_enrich_normalized_capped" = "B/Phu13",
  "CR9114_compitition_Phu_FluB_avg_enrich_normalized_capped" = "B/Phu13+CR9114",
  "Full_HA_Lee_FluB_avg_enrich_normalized_capped" = "B/Lee40",
  "CR9114_compitition_Lee_FluB_avg_enrich_normalized_capped" = "B/Lee40+CR9114"
)

data_capped$Enrichment_Score <- factor(
  data_capped$Enrichment_Score,
  levels = names(y_labels),
  labels = y_labels
)

# Define color limits and legend title
color_limits <- c(-3, 50)
legend_title <- "Normalized\nScore"

# Split data into subsets for visualization
split_data <- split(data_capped, cut(as.numeric(factor(data_capped$name)), 4))

# Plot the heatmaps
p1 <- plot_enrichment_heatmap(split_data[[1]], legend_title, color_limits, y_labels)
p2 <- plot_enrichment_heatmap(split_data[[2]], legend_title, color_limits, y_labels)
p3 <- plot_enrichment_heatmap(split_data[[3]], legend_title, color_limits, y_labels)
p4 <- plot_enrichment_heatmap(split_data[[4]], legend_title, color_limits, y_labels)
#p5 <- plot_enrichment_heatmap(split_data[[5]], legend_title, color_limits, y_labels)

# Combine the plots
final_plot <- grid.arrange(p1, p2, p3, p4, nrow = 4)

# Save the final plot
ggsave("graph/oPool_analysis/oPool_heatmap.png", final_plot, width = 7, height = 10, dpi = 600)

# Print the final plot
print(final_plot)