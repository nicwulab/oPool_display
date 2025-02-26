library(ggplot2)
library(dplyr)
library(readxl)
library(cowplot)

# Adjust time function: Shift Time to start at 0
adjust_time <- function(data) {
  data %>% mutate(Time2 = Time2 - min(Time2, na.rm = TRUE))
}

# Calculate slope function
calculate_slope <- function(time, data) {
  model <- lm(data ~ time)
  return(coef(model)[2])
}

# Ensure output directory exists
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Color mapping for antigens
color_map <- c(
  "SI06_H1" = "#D18FE7",  # Purple
  "MI15_H1" = "#4CA998",  # Teal
  "SP16_H3" = "#E47AB5",  # Pink
  "QH_H5"   = "#4A9FD8",  # Light Blue
  "SH_H7"   = "#5B72D4",  # Blue
  "B_Phu13" = "#B08F3E", # Brown/Gold
  "B_Lee40" = "#C8645D"   # Red
)

# Specific antibody orders for each antigen
antibody_order <- list(
  "SI06_H1" = c("AG11-2F01", "31.a.55", "009-10-1G06", "045-09-1G05", "009-10-2F01", "051-10-2B05"),
  "MI15_H1" = c("31.a.55", "029-09-2A06","009-10-1G06", "045-09-1G05", "042-10-2D01", "009-10-2F01"),
  "QH_H5" = c("3C11", "045-09-1G05", "56.e.01", "C7-3-C02"),
  "SH_H7" = c("045-09-1G05", "042-10-2D01"),
  "B_Phu13" = c("K77-2A04", "017-10116-3D04", "011-10069-2G04", "034-100809-1F05"),
  "B_Lee40" = c("K77-2A04", "011-10069-2G04", "017-10116-3D04", "150055-015-1D02", "034-100809-1F05")
)

# Sensorgram plotting function
plot_sensorgram <- function(data1, data2, antibody, antigen, output_folder) {
  ensure_dir(output_folder)

  ymax <- max(max(data1$Data2, na.rm = TRUE), max(data2$Data2, na.rm = TRUE)) * 1.1
  if (ymax < 1) ymax <- 1
  ymin <- -0.25

  data1 <- data1 %>% mutate(Label = "+CR9114")
  data2 <- data2 %>% mutate(Label = "-CR9114")

  combined_data <- bind_rows(data1, data2)

  p <- ggplot(combined_data, aes(x = Time2, y = Data2, color = Label)) +
    geom_line(linewidth = 0.8) +
    theme_cowplot(12) +
    theme(
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 6, face = "bold"),
      axis.text = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 0.5),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, 'lines'),
      legend.text = element_text(size = 4, face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    ylim(ymin, ymax) +
    xlab(bquote(bold(Time~'(s)'))) +
    ylab(bquote(bold(Response~'(nm)')))

  ggsave(file.path(output_folder, paste0(antigen, "_", antibody, "_Sensorgram.png")), plot = p, width = 1.2, height = 1.2, dpi = 600)
}

# Bar chart plotting function
plot_bar_chart <- function(data, antigen, output_folder) {
  ensure_dir(output_folder)

  bar_width <- 0.35  # Adjusted bar width
  plot_height <- 1.6  # Default plot height

  # Adjust plot height for specific antigens
  if (antigen %in% c("B_Phu13", "B_Lee40")) {
    plot_height <- 1.75
  }

  num_antibodies <- nrow(data)
  min_plot_width <- 0.2  # Minimum width in inches
  max_plot_width <- 3.5  # Maximum width in inches
  plot_width <- min_plot_width + (num_antibodies * bar_width)
  plot_width <- min(plot_width, max_plot_width)  # Cap max width

  # Apply specific antibody order if available
  if (antigen %in% names(antibody_order)) {
    data$antibody <- factor(data$antibody, levels = antibody_order[[antigen]])
  }

  p <- ggplot(data, aes(x = antibody, y = competition, fill = antigen)) +
    geom_bar(stat = "identity", color = "black", width = bar_width, fill = color_map[antigen]) +  # Apply custom color
    scale_y_continuous(limits = c(0, 100), labels = scales::percent_format(scale = 1)) +  # Add % sign
    theme_cowplot(12) +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.text = element_text(size = 7, face = "bold", colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),  
      axis.line = element_blank(),  
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),  
      legend.position = "none",
      plot.margin = margin(5, 5, 10, 5)  # Reduced margins for compactness
    ) +
    labs()

  ggsave(file.path(output_folder, paste0("competition_", antigen, ".png")), plot = p, width = plot_width, height = plot_height, dpi = 600)
}

# Main function
process_data <- function(metadata_file, output_folder) {
  ensure_dir(output_folder)

  metadata <- read_excel(metadata_file)
  competition_data <- list()

  for (i in 1:nrow(metadata)) {
    antibody <- metadata$antibody[i]
    antigen <- metadata$antigen[i]
    file1 <- metadata$data_file_1[i]
    file2 <- metadata$data_file_2[i]

    df1 <- read.table(file1, header = TRUE, sep = "\t", skip = 4, fill = TRUE)
    df2 <- read.table(file2, header = TRUE, sep = "\t", skip = 4, fill = TRUE)

    df1 <- adjust_time(df1)
    df2 <- adjust_time(df2)

    plot_sensorgram(df1, df2, antibody, antigen, output_folder)

    df1_20s <- df1 %>% filter(Time2 <= 60)
    df2_20s <- df2 %>% filter(Time2 <= 60)

    slope1 <- calculate_slope(df1_20s$Time2, df1_20s$Data2)
    slope2 <- calculate_slope(df2_20s$Time2, df2_20s$Data2)

    # Calculate CR9114 competition % with new rules
    competition_percentage <- ifelse(slope2 != 0, 100 - ((slope1 / slope2) * 100), 100)
    competition_percentage <- pmax(pmin(competition_percentage, 100), 0)  # Cap between 0 and 100

    competition_data[[antigen]] <- rbind(competition_data[[antigen]],
                                         data.frame(antibody = antibody, competition = competition_percentage))
  }

  for (antigen in names(competition_data)) {
    plot_bar_chart(competition_data[[antigen]], antigen, output_folder)
  }
}


process_data("experimental_data/validation/BLI/oPool_competition_validation/oPool_competition_validation_sample_names.xlsx", 
    "graph/validation/BLI/oPool_competition_validation")
