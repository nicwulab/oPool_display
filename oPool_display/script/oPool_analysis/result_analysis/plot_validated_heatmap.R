library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)

# Heatmap plotting function
plot_enrichment_heatmap <- function(data, legend_title, color_limits, y_labels, textsize = 7) {
  p <- ggplot(data, aes(x = name, y = Enrichment_Score, fill = parameter)) +
    geom_tile(color = "black", size = 0.3) +
    scale_fill_gradientn(
      colours = c("blue","lightblue", "white", "yellow","yellow", "orange","orange", "red", "red"),
      limits = color_limits,
      breaks = c(0, 25, 50),
      labels = c("<0", "25", ">50"),
      guide = "colorbar",
      na.value = "white"
    ) +
    coord_fixed(ratio = 1) +
    scale_y_discrete(labels = y_labels) +
    theme_cowplot(12) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.text = element_text(size = 7, face = "bold", colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 7, face = "bold", colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_blank(),  
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
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

# Identify enrich columns
enrich_cols <- grep("avg_enrich", names(combined_data), value = TRUE)

# Normalize enrichment columns using median and IQR
data_normalized <- combined_data %>%
  mutate(across(
    all_of(enrich_cols),
    ~ (.-median(., na.rm = TRUE)) / IQR(., na.rm = TRUE),
    .names = "{.col}_normalized"
  ))

# Define cutoff thresholds
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

# Apply capping
data_capped_wide <- data_normalized %>%
  mutate(across(
    all_of(names(normalized_cutoffs)),
    ~ ifelse(. < normalized_cutoffs[cur_column()], 0, ifelse(. > 50, 50, .)),
    .names = "{.col}_capped"
  ))

# Filter to selected scFvs
selected_abs <- c(
  "31.a.55", "16.ND.92", "045-09-1G05", "C3-2-F02", "240-14-IgA-2F02", "56.j.01", "01.ad.01", 
  "042-100809-2F04", "009-10-2F01", "009-10-1G06", "029-09-2A06", "AG11-2F01", "AG11-3C06", 
  "56.e.01", "051-10-2B05", "C7-3-C02", "011-10069-2G04", "017-10116-3D04", "042-10-2D01", 
  "K77-2A04", "3C11", "034-100809-1F05", "AG13-3G01", "150055-015-1D02", "AG2-G02"
)
data_capped_wide <- data_capped_wide %>%
  filter(name %in% selected_abs)

# Pivot longer
data_capped <- data_capped_wide %>%
  pivot_longer(
    cols = ends_with("_capped"),
    names_to = "Enrichment_Score",
    values_to = "parameter"
  ) %>%
  drop_na() %>%
  filter(!grepl("CR9114", Enrichment_Score))

# Label mapping
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

# Plot single heatmap row
legend_title <- "Normalized\nScore"
color_limits <- c(-3, 50)

heatmap_plot <- plot_enrichment_heatmap(data_capped, legend_title, color_limits, y_labels)

# Save plot
ggsave("graph/oPool_analysis/oPool_heatmap_selected_scfvs.png", heatmap_plot, width = 6.5, height = 4, dpi = 600)

# Display plot
print(heatmap_plot)