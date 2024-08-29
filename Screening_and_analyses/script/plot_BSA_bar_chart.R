library(ggplot2)
library(cowplot)
library(readr)
library(scales)

# Define the function to create and save the stacked bar chart
create_stacked_bar_chart <- function(df, colors, output_file, width, height) {
  df$Sample <- factor(df$Sample, levels = c("AG11-2F01", "16.ND.92", "MEDI-8852", "56.a.09", "54-1G05", "PN-SIA28", "39.29", "429 B01"))
  df$Category <- factor(df$Category, levels = c( "VL","Non-CDR H3 VH","CDR H3"))
  # Create the stacked bar chart
  custom_labels <- c("CDR H3" = "CDRH3", 
                     "Non-CDR H3 VH" = expression(bold("Non-CDRH3 V"[H])), 
                     "VL" = expression(bold("V"[L])))

  stacked_bar_plot <- ggplot(df, aes(x = Sample, y = Percentage, fill = Category)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = colors, labels = custom_labels) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text = element_text(size = 7, face = "bold", colour = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
          axis.text.y = element_text(hjust = 1, face = "bold", colour = 'black'),
          axis.title = element_text(size = 7, face = "bold"),
          axis.line = element_line(linewidth = 0.5, color = "black"), # Bold axis lines
          panel.border = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face = "bold", size = 6, hjust = 0), # Bold legend text
          legend.key.size = unit(0.3, "cm"), # Smaller legend key size
          legend.spacing.x = unit(0, "cm"), # Reduce spacing between legend items and text
          legend.spacing.y = unit(0.2, "cm"),
          legend.position=c(0.2, -0.8),
          legend.direction = "horizontal",
          legend.box.margin = margin(t = -10), # Adjust margin to bring legend closer to plot
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) + # Bold legend text
    labs(x = "", y = "Paratope BSA")

  # Print the plot
  print(stacked_bar_plot)
  
  # Save the stacked bar chart
  ggsave(output_file, plot = stacked_bar_plot, width = width, height = height)
}

# Read the data from the TSV file
input_file <- "data/bsa_percentage.tsv"
data <- read_tsv(input_file)

# Define colors for the categories
colors <- c("CDR H3" = "orange", "Non-CDR H3 VH" = "#FFD700", "VL" = "pink")

# Call the function to create and save the stacked bar chart
create_stacked_bar_chart(data, colors, "graph/BSA_stacked_bar_chart.png", 3, 1.5)