library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(cowplot)

plot_elisa_heatmap <- function(excel_file, output_image) {
  # Load the data
  data <- read_excel(excel_file, sheet = "Sheet1")
  print(colnames(data))
  
  # Process data for heat map
  data <- data %>%
    select(Antibody, Antigen, Avg_450) %>%
    filter(!Antibody %in% c("CR9114", "PBS only control", "037-10-5E04")) %>%
    drop_na()
  
  # Convert to long format for ggplot
  data_long <- data %>%
    mutate(Antigen = factor(Antigen, levels = c("H1 stem", "H3 stem", "H1/SI06", "H1/MI15", "H3/SP16", "H5/QH05", "H7/SH13", "B/Phu13", "B/Lee40")))
  
  # Define list for overlay dots
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
  
  # Define gradient colors
  gradient_colors <- c("white", "#2171b5")
  
  # Plot heat map
  p <- ggplot(data_long, aes(x = Antibody, y = Antigen, fill = Avg_450)) +
    geom_tile(color = "black", size = 0.3) +
    geom_tile(data = filter(data_long, Avg_450 < 0.5), fill = "white", color = "black", size = 0.3) +
    geom_tile(data = filter(data_long, Avg_450 > 2), fill = "#2171b5", color = "black", size = 0.3) +
    geom_point(
      data = overlay_dots,
      aes(x = Antibody, y = Antigen),
      size = 0.8,
      shape = 21,
      fill = "black",
      color = "black"
    ) +
    scale_fill_gradientn(
      colours = gradient_colors,
      limits = c(0, 2),
      breaks = c(0, 1, 2),
      labels = c("0", "1", ">2"),
      na.value = "white",
      guide = guide_colorbar(
        title = "Avg_450",
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
  
  # Save the plot
  ggsave(output_image, plot = p, height = 4, width = 6.5, dpi = 300)
}


infile = "experimental_data/validation/ELISA/ELISA_Validation_Results.xlsx"
outfile = "graph/validation/ELISA/oPool_validation_heatmap_ELISA.png"
plot_elisa_heatmap(infile, outfile)
