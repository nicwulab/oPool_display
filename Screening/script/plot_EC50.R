# Load necessary libraries
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(stringr)
library(dplyr)
library(gridExtra)
library(cowplot)

# Create a data frame with the provided values
df_new <- data.frame(
  Antigen = c("H1 USSR/77", "H1 Taiwan/86", "H1 Beijing/95", "H1 SI/06", "H1 CA/09", 
              "H1 Mich/15", "H1 Bris/18", "H5 Wash/14","H5 Missouri/15","H5 Florida/22","H1 stem","H3 stem"),
  `AG11-2F01` = c(0.3843, 0.7035, 11.2, 8.013, 0.3408, 0.5358, 1.1, 3.774, 1.304, 0.1, 1.04, 1000),
  `16.ND.92` = c(0.497, 1.699, 0.7202, 0.5899, 1.646, 1.724, 7.714, 8.623, 2.431, 0.6305, 0.7021,1000),
  check.names = FALSE
)


# Reshape the data to long format
df_new_long <- melt(df_new, id.vars = "Antigen", variable.name = "Sample", value.name = "Concentration")

# Ensure the correct order of factors for the y-axis
df_new_long$Sample <- factor(df_new_long$Sample, levels = c("16.ND.92","AG11-2F01"))

df_new_long$Antigen <- factor(df_new_long$Antigen, levels = c("H1 USSR/77", "H1 Taiwan/86", "H1 Beijing/95", 
              "H1 SI/06", "H1 CA/09", 
              "H1 Mich/15", "H1 Bris/18", "H5 Wash/14","H5 Missouri/15","H5 Florida/22","H1 stem","H3 stem"))

# Generate the heat map with updated theme settings
heatmap_plot_new <- ggplot(df_new_long, aes(x = Antigen, y = Sample, fill = Concentration)) +
  geom_tile(color = "black", linewidth = 0.8) + # Use linewidth instead of size
  scale_fill_gradientn(colors = c("#08306b", "#2171b5", "#6baed6", "#bdd7e7", "white"),
                       trans = 'log',  # Apply log transformation
                       limits = c(0.1, 1000),  # Set the limits from 0.1 to 100
                       breaks = c(0.1, 1, 10, 100, 1000),  # Set log scale breaks
                       labels = c("0.1","1", "10", "100",">1000"),  # Set log scale labels
                       na.value = "grey",  # Fill grey for NA values
                       guide = guide_colorbar(barwidth = 0.6, barheight = 4, ticks.length = unit(1, "cm"),ticks.linewidth = 0.8,
                                              frame.colour = "black", frame.linewidth = 0.8,
                                              title = NULL,
                                              label.theme = element_text(size = 6, face = "bold", colour = 'black'),
                                              ticks.colour = "black")) +
  theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        axis.text = element_text(size = 7, face = "bold", colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        axis.text.y = element_text(hjust = 1, colour = 'black'),
        axis.title = element_text(size = 7, face = "bold"),
        axis.line = element_blank(), # Remove axis lines
        panel.border = element_blank(), # Remove panel border
        plot.margin = unit(c(0, 0, 0, 0), "cm")) + # Remove plot margins
  xlab("") +
  ylab("") +
  coord_fixed() # Ensures the aspect ratio is 1:1, making tiles square

# Print the plot
print(heatmap_plot_new)

# Save the heat map with smaller boxes
ggsave("graph/elisa_ec50_heatmap.png", plot = heatmap_plot_new, width = 4, height = 2.5)
