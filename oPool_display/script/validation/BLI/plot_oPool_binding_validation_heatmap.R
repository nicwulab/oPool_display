library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(cowplot)

# Function to calculate response differences
calculate_differences <- function(data, antibody, antigen) {
  start_time <- min(data$Time1, na.rm = TRUE)+1
  adjusted_end_time <- start_time + 55
  
  response_start <- data %>% filter(Time1 == start_time) %>% pull(Data1)
  response_end <- data %>% filter(Time1 == adjusted_end_time) %>% pull(Data1)
  
  if (length(response_start) == 1 & length(response_end) == 1) {
    difference <- response_end - response_start
  } else {
    difference <- NA
  }
  
  return(data.frame(Antibody = antibody, Antigen = antigen, Difference = difference))
}

# Function to plot the heatmap
plot_heatmap <- function(differences, output_file, legend_title = "Response\n(nm)", overlay_dots = NULL) {
  differences <- differences %>%
    mutate(Difference = as.numeric(Difference)) %>%
    filter(!is.na(Difference))
  
  heatmap_data <- differences %>%
    pivot_wider(names_from = Antibody, values_from = Difference, values_fill = NA)
  
  heatmap_long <- heatmap_data %>%
    pivot_longer(cols = -Antigen, names_to = "Antibody", values_to = "Difference") %>%
    mutate(
      FillColor = ifelse(Difference > 1, "red", NA),
      Antigen = factor(Antigen, levels = c(
        "H1 stem", "H3 stem", "H1/SI06", "H1/MI15", "H3/SP16", "H5/QH05", 
        "H7/SH13", "B/Phu13", "B/Lee40"
      ))
    )
  
  gradient_colors <- c("white", "red")
  
  p <- ggplot(heatmap_long, aes(x = Antibody, y = Antigen, fill = Difference)) +
    geom_tile(color = "black", size = 0.3) +
    geom_tile(data = heatmap_long %>% filter(Difference > 1), fill = "red", color = "black", size = 0.3) +
    geom_tile(data = heatmap_long %>% filter(Difference < 0.08), fill = "white", color = "black", size = 0.3) +
    scale_fill_gradientn(
      colours = gradient_colors,
      limits = c(-0.05, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", ">1"),
      na.value = "white",
      guide = guide_colorbar(
        title = legend_title,
        title.theme = element_text(size = 7, face = "bold", colour = "black", hjust = 0.5),
        label.theme = element_text(size = 7, face = "bold", colour = "black"),
        barwidth = 0.6, 
        barheight = 5,
        frame.colour = "black",
        frame.linewidth = 0.75,
        ticks = TRUE,
        ticks.colour = "black"
      )
    ) +
    coord_fixed(ratio = 1) +
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
    xlab("") +
    ylab("")
  
  if (!is.null(overlay_dots)) {
    p <- p + 
      geom_point(
        data = overlay_dots,
        aes(x = Antibody, y = Antigen),
        size = 0.8,
        shape = 21,
        fill = "black",
        color = "black",
        inherit.aes = FALSE
      )
  }
  
  ggsave(filename = output_file, plot = p, height = 4, width = 6.5, dpi = 300)
}

# Process metadata and calculate differences
process_metadata_with_heatmap <- function(metadata_file, output_folder) {
  metadata <- read_excel(metadata_file)
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  differences <- list()
  
  for (i in 1:nrow(metadata)) {
    antibody <- metadata$antibody[i]
    antigen <- metadata$antigen[i]
    data_file <- metadata$data_file[i]
    
    raw_data <- read.table(
      data_file,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      skip = 4,
      check.names = FALSE
    )
    time_series <- raw_data %>%
      select(Time1, Data1) %>%
      filter(!is.na(Time1) & !is.na(Data1)) %>%
      mutate(Antigen = antigen, Antibody = antibody)
    
    diff <- calculate_differences(time_series, antibody, antigen)
    differences[[length(differences) + 1]] <- diff
  }
  
  differences_df <- bind_rows(differences)
  
  # List-like input for dot overlay
  dot_tiles <- list(
    "009-10-1G06" = c("H3 stem","H1/SI06", "H1/MI15"),
    "009-10-2F01" = c("H3 stem","H1/SI06", "H1/MI15"),
    "01.ad.01" = c("H3 stem","H3/SP16"),
    "011-10069-2G04" = c("B/Phu13", "B/Lee40"),
    "017-10116-3D04" = c("B/Phu13", "B/Lee40"),
    "029-09-2A06" = c("H3 stem", "H1/MI15"),    
    "042-10-2D01" = c("H1/MI15", "H7/SH13"),
    "042-100809-2F04" = c("H3 stem","H3/SP16"),
    "045-09-1G05" = c("H1/SI06", "H1/MI15","H7/SH13"),
    "051-10-2B05" = c("H3 stem","H1/SI06"),
    "16.ND.92" = c("H1 stem"), 
    "240-14-IgA-2F02" = c("H3 stem","H3/SP16"),  
    "31.a.55" = c("H1 stem", "H3 stem", "H1/SI06", "H1/MI15", "H3/SP16"),
    "AG11-2F01" = c("H1 stem", "H1/SI06"),
    "56.e.01" = c("H3/SP16", "H5/QH05"),
    "56.j.01" = c("H3 stem","H3/SP16"), 
    "AG13-3G01" = c("H5/QH05", "B/Phu13"),
    "AG2-G02" = c("H3 stem"),
    "C7-3-C02" = c("H1/MI15","H5/QH05"),
    "K77-2A04" = c("B/Phu13", "B/Lee40"),
    "AG11-3C06" = c("H1/MI15"),
    "034-100809-1F05" = c("B/Phu13", "B/Lee40"),
    "3C11" = c("H5/QH05"),
    "150055-015-1D02" = c("B/Lee40")
  )
  
  # Convert list to data frame for overlay dots
  overlay_dots <- bind_rows(lapply(names(dot_tiles), function(antibody) {
    data.frame(Antibody = antibody, Antigen = dot_tiles[[antibody]])
  }))
  
  heatmap_file <- paste0(output_folder, "/BLI_binding_response_heatmap.png")
  plot_heatmap(differences_df, heatmap_file, overlay_dots = overlay_dots)
}

# Example usage
metadata_file <- "experimental_data/validation/BLI/oPool_binding_validation/oPool_validation_sample_names.xlsx"
output_folder <- "graph/validation/BLI/oPool_binding_validation"
process_metadata_with_heatmap(metadata_file, output_folder)