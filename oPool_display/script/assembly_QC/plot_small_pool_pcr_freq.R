# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

# Function to generate plots for a single file
plot_frequencies <- function(infile, outfile, textsize = 10, filter_or_not = "No", values_to_filter = NA) {
  # Read the data
  df <- readr::read_tsv(infile, show_col_types = FALSE)
  
  # Dynamically identify columns for Rep1 and Rep2 frequencies
  x_col <- grep("Rep1.*_freq", colnames(df), value = TRUE)
  y_col <- grep("Rep2.*_freq", colnames(df), value = TRUE)
  
  # Ensure both columns are found
  if (length(x_col) != 1 || length(y_col) != 1) {
    stop("Could not identify unique columns for Rep1_freq and Rep2_freq.")
  }
  
  # Optionally filter the data
  if (filter_or_not == "Yes") {
    df <- df %>%
      filter(.data[['name']] %in% values_to_filter)
  }
  
  # Create the dot plot
  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(size = 0.5, pch = 16, alpha = 0.3) +
    geom_smooth(method = "lm", color = "red", size = 0.5, se = FALSE, linetype = "dashed") + # Fit line
    scale_x_log10(breaks = 10^seq(-5, 0, by = 1), 
                  labels = c(-5, -4, -3, -2, -1, 0), 
                  limits = c(0.00001, 0.3)) +
    scale_y_log10(breaks = 10^seq(-5, 0, by = 1), 
                  labels = c(-5,-4, -3, -2, -1, 0), 
                  limits = c(0.00001, 0.3)) +
    theme_cowplot(12) +
    theme(
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white"),
      axis.title = element_blank(),
      axis.text = element_text(size = textsize, face = "bold", family = "Arial"),
      legend.key.size = unit(0.1, 'in'),
      legend.spacing.x = unit(0.03, 'in'),
      legend.title = element_blank(),
      legend.text = element_text(size = textsize - 2, face = "bold", family = "Arial"),
      legend.position = 'right'
    ) +
    guides(legend = guide_legend(nrow = 2, byrow = TRUE))  + # Wrap legend text by rows
    labs(
     x = bquote(bold(Log[10]~"scFv frequency (Rep 1)")),
     y = bquote(bold(Log[10]~"scFv frequency (Rep 2)"))
    )
  
  # Save the plot
  ggsave(outfile, p, height = 0.8, width = 0.8, dpi = 600)
  
  # Transform values to log10 scale
  log_x <- log10(df[[x_col]])
  log_y <- log10(df[[y_col]])

  # Remove -Inf values that can result from log10(0) or invalid data
  valid_indices <- is.finite(log_x) & is.finite(log_y)
  log_x <- log_x[valid_indices]
  log_y <- log_y[valid_indices]

  # Calculate the Pearson correlation
  pearson_corr <- cor(log_x, log_y, method = "pearson")

  # Print the correlation
  print(paste("Pearson correlation between replicates (log scale):", pearson_corr))
  }

# Directory containing Python output files
python_output_dir <- "oPool_result/assembly_QC/freq" # Replace with your directory path

# Directory to save the plots
output_plot_dir <- "graph/assembly_QC" # Replace with your directory path
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir)

# List all .tsv files in the Python output directory
python_output_files <- list.files(python_output_dir, pattern = "*.tsv", full.names = TRUE)

# Generate plots for each Python output file
for (infile in python_output_files) {
  # Define the output file name for the plot
  plot_file <- file.path(output_plot_dir, paste0(tools::file_path_sans_ext(basename(infile)), "_plot.png"))
  
  # Call the plot function
  plot_frequencies(
    infile = infile,
    outfile = plot_file,
    textsize = 5,
    filter_or_not = "No"
  )
  
  message("Plot generated for: ", basename(infile))
}