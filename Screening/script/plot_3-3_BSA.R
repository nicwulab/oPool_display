library(ggplot2)
library(cowplot)
library(readr)

# Define the function to create and save the bar chart
create_bar_chart <- function(df, color, output_file, width, height) {
  df$Abs <- factor(df$Abs, levels = c("AG11-2F01", "16.ND.92", "MEDI-8852", "56.a.09", "54-1G05", "PN-SIA28", "39.29", "429 B01"))
  
  # Create the bar chart
  bar_plot <- ggplot(df, aes(x = Abs, y = `IGHD3-3%`)) +
    geom_bar(stat = "identity", fill = color, color = "black", linewidth = 0.5) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text = element_text(size = 7, face = "bold", colour = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
          axis.text.y = element_text(hjust = 1, face = "bold", colour = 'black'),
          axis.title = element_text(size = 7, face = "bold"),
          axis.line = element_line(linewidth = 0.5, color = "black"), # Bold axis lines
          panel.border = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    labs(x = "", y = expression(bold("V"[H]~" BSA")))

  # Print the plot
  print(bar_plot)
  
  # Save the bar chart
  ggsave(output_file, plot = bar_plot, width = width, height = height)
}

# Read the data from the TSV file
input_file <- "data/3-3_BSA_percentage.tsv"
data <- read_tsv(input_file)

# Define the color for CDR H3 from the previous script
cdr_h3_color <- "orange"

# Call the function to create and save the bar chart
create_bar_chart(data, cdr_h3_color, "graph/BSA_IGHD3-3_bar_plot.png", 3, 1.4)