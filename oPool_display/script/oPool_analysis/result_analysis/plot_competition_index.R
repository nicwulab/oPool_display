# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggpubr)

# Define antigens and corresponding file names
antigens <- c("SI06_H1", "MI15_H1", "SP16_H3", "QH_H5", "SH_H7", "Phu_FluB", "Lee_FluB")
file_paths <- paste0("oPool_result/filtered_hits/filtered_antigen_", antigens, "_custom_cutoff.tsv")

# Define matching colors based on extracted legend
color_map <- c(
  "SI06_H1" = "#D18FE7",
  "MI15_H1" = "#4CA998",
  "SP16_H3" = "#E47AB5",
  "QH_H5"   = "#4A9FD8",
  "SH_H7"   = "#5B72D4",
  "Phu_FluB" = "#B08F3E",
  "Lee_FluB" = "#C8645D"
)

# Plot dimensions for filtered antibodies
plot_heights_filtered <- list(
  "SI06_H1" = 1.2, "MI15_H1" = 1.2, "SP16_H3" = 1.2, "QH_H5" = 1.1,
  "SH_H7" = 1.2, "Phu_FluB" = 1.3, "Lee_FluB" = 1.3
)

# Y-axis limits
y_axis_limits <- list(
  "SI06_H1" = c(-1.5, 1.5), "MI15_H1" = c(-2.0, 2.0), "QH_H5" = c(-1, 1),
  "SH_H7" = c(-1.5, 1.5), "Phu_FluB" = c(-1.5, 1.5), "Lee_FluB" = c(-1.5, 1.5)
)

# Specific antibodies for filtered plots
antibody_lists <- list(
  "SI06_H1" = c("009-10-1G06", "009-10-2F01", "029-09-2A06", "045-09-1G05", "051-10-2B05", "16.ND.92", "31.a.55", "AG11-2F01"),
  "MI15_H1" = c("009-10-1G06", "009-10-2F01", "029-09-2A06", "042-10-2D01","045-09-1G05", "051-10-2B05", "16.ND.92", "31.a.55", "AG11-2F01"),  
  "QH_H5"   = c("045-09-1G05", "3C11", "56.e.01", "C7-3-C02"),
  "SH_H7"   = c("042-10-2D01","045-09-1G05"),
  "Phu_FluB" = c("011-10069-2G04", "017-10116-3D04", "034-100809-1F05", "150055-015-1D02", "K77-2A04"),
  "Lee_FluB" = c("011-10069-2G04", "017-10116-3D04", "034-100809-1F05", "150055-015-1D02", "K77-2A04")
)

# Initialize an empty dataframe to store filtered competition index results
competition_results <- data.frame(Antigen = character(), Antibody = character(), Competition_Index = numeric(), stringsAsFactors = FALSE)

# Read and process each file
for (i in seq_along(file_paths)) {
  file_path <- file_paths[i]
  antigen <- antigens[i]

  if (!file.exists(file_path)) {
    print(paste("Skipping missing file:", file_path))
    next
  }

  data <- read_tsv(file_path, show_col_types = FALSE)
  
  full_ha_col <- paste0("Full_HA_", antigen, "_avg_enrich")
  cr9114_col <- paste0("CR9114_compitition_", antigen, "_avg_enrich")

  if (!(full_ha_col %in% names(data)) || !(cr9114_col %in% names(data))) {
    print(paste("Skipping", antigen, "- Required columns missing"))
    next
  }

  data <- data %>%
    mutate(competition_index = log10(!!sym(full_ha_col) / !!sym(cr9114_col))) %>%
    select(closest_abs, competition_index) %>%
    tidyr::drop_na() %>%
    rename(Antibody = closest_abs, Competition_Index = competition_index) %>%
    mutate(Antigen = antigen)

  # Filter only selected antibodies per antigen
  if (antigen %in% names(antibody_lists)) {
    filtered_data <- data %>% filter(Antibody %in% antibody_lists[[antigen]])
  }

  # Append to results
  competition_results <- bind_rows(competition_results, data)
}

# Save the filtered competition index table
output_file <- "oPool_result/enrichment/competition_index.tsv"
write.table(competition_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

print(paste("Filtered competition index data saved to:", output_file))

# Function to generate bar plots
generate_plot <- function(data, antigen, output_file, plot_height) {
  ymin <- y_axis_limits[[antigen]][1]
  ymax <- y_axis_limits[[antigen]][2]

  bar_width <- 0.3  
  num_antibodies <- nrow(data)
  min_plot_width <- 0.2
  max_plot_width <- 6.5
  plot_width <- min_plot_width + (num_antibodies * bar_width*0.8)
  plot_width <- min(plot_width, max_plot_width)

  p <- ggplot(data, aes(x = reorder(Antibody, -Competition_Index), y = Competition_Index)) +
    geom_bar(stat = "identity", color = "black", fill = color_map[antigen], width = bar_width, linewidth = 0.5) +  
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "solid") +  
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, length.out = 5), labels = scales::label_number(accuracy = 0.1)) +  
    theme_cowplot(12) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.text = element_text(size = 6, face = "bold", colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position = "none",
      plot.margin = margin(5, 5, 10, 5)
    )

  ggsave(output_file, plot = p, width = plot_width, height = plot_height, dpi = 600)
  print(paste("Saved:", output_file))
}

# Generate plots for filtered antibodies only
for (i in seq_along(antigens)) {
  antigen <- antigens[i]
  data <- competition_results %>% filter(Antigen == antigen)

  if (nrow(data) > 0) {
    output_file <- paste0("graph/oPool_analysis/competition/filtered/competition_index_", antigen, "_filtered_antibodies.png")
    generate_plot(filtered_data, antigen, output_file, plot_heights_filtered[[antigen]])
  } else {
    print(paste("No matching antibodies for antigen:", antigen))
  }

  if (nrow(data) > 0) {
    output_file_all <- paste0("graph/oPool_analysis/competition/all/competition_index_", antigen, "_all_antibodies.png")
    generate_plot(data, antigen, output_file_all, plot_heights_all[[antigen]])
  }
}