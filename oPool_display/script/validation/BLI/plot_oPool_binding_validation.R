library(ggplot2)
library(dplyr)
library(readxl)
library(cowplot)

# Define color map **globally** so all functions can access it
color_map <- c(
  "H1 stem" = "#E6194B",
  "H3 stem" = "#6ACC65", 
  "H1/SI06" = "#D18FE7",
  "H1/MI15" = "#4CA998",
  "H3/SP16" = "#FFB347",
  "H5/QH05" = "#4A9FD8",
  "H7/SH13" = "#5B72D4",
  "B/Phu13" = "#B08F3E",
  "B/Lee40" = "#C8645D"
)

# Adjust time function: Shift Time to start at 0
adjusttime <- function(data) {
  data %>%
    mutate(Time1 = Time1 - min(Time1, na.rm = TRUE))
}

# Sensorgram plotting function
plot_antibody <- function(data, antibody, output_folder, color_map) {
  textsize <- 6
  
  # Dynamically calculate ymax with 10% buffer
  ymax <- max(data$Data1, na.rm = TRUE) * 1.1
  ymax <- max(ymax, 1)  # Ensure ymax is at least 1
  

  # Trim whitespace and ensure correct factor levels
  data <- data %>%
    mutate(Antigen = trimws(Antigen)) %>%
    mutate(Antigen = factor(Antigen, levels = names(color_map)))

  # Set ymin
  ymin <- -0.25

  # Create the plot
  p <- ggplot(data, aes(x = Time1, y = Data1, color = Antigen, group = Antigen)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = color_map, drop = FALSE) +  
    theme_cowplot(12) +
    theme(
      plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = textsize, face = "bold"),
      axis.title.x = element_text(size = textsize, face = "bold"),
      axis.text.y = element_text(size = textsize, face = "bold"),
      axis.text.x = element_text(size = textsize, face = "bold", angle = 90),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, 'lines'),
      legend.text = element_text(size = textsize - 2, face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylim(ymin, ymax) +  
    xlab(bquote(bold(Time~'(s)'))) +
    ylab(bquote(bold(Response~'(nm)'))) +
    ggtitle(paste(antibody))

  # Save the plot
  output_file <- file.path(output_folder, paste0(antibody, "_Sensorgram.png"))
  ggsave(filename = output_file, plot = p, height = 1.3, width = 1.3, dpi = 600)
}

# Process metadata and plot for each antibody
process_metadata <- function(metadata_file, output_folder, color_map) {
  # Read metadata
  metadata <- read_excel(metadata_file)

  # Create output folder if it does not exist
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  # Process each file in the metadata
  combined_data <- list()
  for (i in seq_len(nrow(metadata))) {
    antibody <- metadata$antibody[i]
    antigen <- metadata$antigen[i]
    data_file <- file.path(metadata$data_file[i])

    # Read and process data
    raw_data <- read.table(
      data_file, sep = "\t", header = TRUE, fill = TRUE, skip = 4, check.names = FALSE
    )

    time_series <- raw_data %>%
      select(Time1, Data1) %>%
      filter(!is.na(Time1) & !is.na(Data1)) %>%
      adjusttime() %>%
      mutate(Antigen = trimws(antigen), Antibody = antibody)

    # Store in list for combining later
    combined_data[[antibody]] <- bind_rows(combined_data[[antibody]], time_series)
  }

  # Generate plots for each antibody
  lapply(names(combined_data), function(antibody) {
    plot_antibody(combined_data[[antibody]], antibody, output_folder, color_map) 
  })
}

# Example usage
metadata_file <- "experimental_data/validation/BLI/oPool_binding_validation/oPool_validation_sample_names.xlsx"
output_folder <- "graph/validation/BLI/oPool_binding_validation"
process_metadata(metadata_file, output_folder, color_map)