library(ggplot2)
library(cowplot)
library(readr)  # for reading TSV files

# Function to plot with the given style
plot_binding_rate <- function(df, output_file) {
  
  # Extracting data for Standard and Unknown
  standard_data <- df[df$Type == "Standard", ]
  unknown_data <- df[df$Type == "Unknown", ]
  
  textsize <- 7
  
  p <- ggplot() +
    geom_point(data = standard_data, aes(x = `Known Conc.`, y = `Binding Rate`), size = 1.5, alpha= 1, color = 'blue') +
    geom_smooth(data = standard_data, aes(x = `Known Conc.`, y = `Binding Rate`), method = "lm", se = FALSE, size = 0.5, color = 'black') +
    geom_point(data = unknown_data, aes(x = `Well Conc.`, y = `Binding Rate`), size = 1.5, alpha= 1, color = 'darkgrey') +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = textsize, face = "bold"),
          axis.title.x = element_text(size = textsize, face = "bold"),
          axis.text.y = element_text(size = textsize, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_text(size = textsize, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'lines'),
          legend.text = element_text(size = textsize, face = "bold"),
          legend.position = 'none',
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Concentration (µg/ml)",
         y = "Binding Rate (nm/s)")
  
  # Save the plot to a file
  ggsave(filename = output_file, plot = p, height = 1.8, width = 3, dpi = 600)
}

# Reading the TSV file
df <- read_tsv("data/BLI_data/scFv_quantitation.tsv")

# Specify the output file name
output_file <- "scFv_quantitation_curve.png"

# Generate the plot
plot_binding_rate(df, output_file)