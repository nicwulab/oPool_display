# Load necessary libraries
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(cowplot)

# Function to create heatmaps
create_heatmap <- function(df, colors, output_file, width, height) {
  # Melt the data frame to long format
  df_long <- melt(df, id.vars = "Antibody", variable.name = "Measurement", value.name = "Value")
  
  # Adjust y-axis labels with subscript
  df_long$Measurement <- factor(df_long$Measurement, levels = c("scFv_KD_nM", "Fab_KD_nM"))
  
  # Generate the heatmap
  heatmap_plot <- ggplot(df_long, aes(x = Antibody, y = Measurement, fill = Value)) +
    geom_tile(color = "black", linewidth = 0.8) + # Use linewidth instead of size
    scale_fill_gradientn(colors = colors,  # Color gradient
                         trans = 'log',  # Apply log transformation
                         limits = c(1, 10000),  # Set the limits from 1 to 10000
                         breaks = c(1, 10, 100, 1000, 10000),  # Set log scale breaks
                         labels = c("<1", "10", "100", "1000", ">10000"),  # Set log scale labels
                         na.value = "grey",  # Fill grey for NA values
                         guide = guide_colorbar(barwidth = 0.6, barheight = 4, ticks.length = unit(1, "cm"),ticks.linewidth = 0.8,
                                                frame.colour = "black", frame.linewidth = 0.8,
                                                title = NULL,
                                                label.theme = element_text(size = 6, face = "bold", colour = 'black'),
                                                ticks.colour = "black")) +
    scale_y_discrete(labels = c(scFv_KD_nM = expression(bold(scFv~K[D])), Fab_KD_nM = expression(bold(Fab~K[D])))) +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text = element_text(size = 7, face = "bold", colour = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
          axis.text.y = element_text(hjust = 1, face = "bold", colour = 'black'),
          axis.title = element_text(size = 7, face = "bold"),
          axis.line = element_blank(), # Remove axis lines
          panel.border = element_blank(), # Remove panel border
          plot.margin = unit(c(0, 0, 0, 0), "cm")) + # Remove plot margins
    xlab("") +
    ylab("") +
    coord_fixed() # Ensures the aspect ratio is 1:1, making tiles square
  
  # Print the plot
  print(heatmap_plot)
  
  # Save the heat map
  ggsave(output_file, plot = heatmap_plot, width = width, height = height)
}

# Data for the blue table
df_blue <- data.frame(
  Antibody = c("31.a.55", "AG11-2F01", "16.ND.92", "C3-2.3F02", "AG11-3C06"),
  scFv_KD_nM = c(32.1, 3.48, 1040, 10000, 10000),
  Fab_KD_nM = c(159, 30.3, 160, 10000, 10000),
  check.names = FALSE
)

# Data for the yellow table
df_yellow <- data.frame(
  Antibody = c("240-14 IgA 2F02", "56.j.01", "01.ad.01", "31.a.55", "042-100809_2F04", "01.h.02", "C3-2.3F02", "AG11-3C06"),
  scFv_KD_nM = c(37.2, 32.4, 3.6, 41.2, 131, 5.14, 10000, 10000),
  Fab_KD_nM = c(37.7, 8.8, 1, 260, 12.4, 1, 10000, 10000),
  check.names = FALSE
)

# Color gradient to be used for both heatmaps
H1_colors <- c("#08306b", "#2171b5", "#bdd7e7", "white")
H3_colors <- c("#ff7400","orange","#ffd700", "white")

# Create heatmaps with adjustable dimensions
create_heatmap(df_blue, H1_colors, "graph/H1_validation_heatmap.png", width = 2.8, height = 2.5)
create_heatmap(df_yellow,H3_colors, "graph/H3_validation_heatmap.png", width = 3.7, height = 2.5)