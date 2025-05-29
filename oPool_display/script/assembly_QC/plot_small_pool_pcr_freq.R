library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

plot_frequencies <- function(infile, outfile, data_outfile = NULL, textsize = 10, filter_or_not = "No", values_to_filter = NA) {

  df <- readr::read_tsv(infile, show_col_types = FALSE)

  x_col <- grep("Rep1.*_freq", colnames(df), value = TRUE)
  y_col <- grep("Rep2.*_freq", colnames(df), value = TRUE)

  if (length(x_col) != 1 || length(y_col) != 1) {
    stop("Could not identify unique columns for Rep1_freq and Rep2_freq.")
  }

  if (filter_or_not == "Yes") {
    df <- df %>%
      filter(.data[['name']] %in% values_to_filter)
  }

  # Save filtered/plotted data
  #if (!is.null(data_outfile)) {
    #write_tsv(df, data_outfile)
  #}

  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(size = 0.5, pch = 16, alpha = 0.3) +
    geom_smooth(method = "lm", color = "red", size = 0.5, se = FALSE, linetype = "dashed") +
    scale_x_log10(breaks = 10^seq(-5, 0, by = 1),
                  labels = c(-5, -4, -3, -2, -1, 0),
                  limits = c(0.00001, 0.3)) +
    scale_y_log10(breaks = 10^seq(-5, 0, by = 1),
                  labels = c(-5, -4, -3, -2, -1, 0),
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
    guides(legend = guide_legend(nrow = 2, byrow = TRUE))  +
    labs(
      x = bquote(bold(Log[10]~"scFv frequency (Rep 1)")),
      y = bquote(bold(Log[10]~"scFv frequency (Rep 2)"))
    )

  ggsave(outfile, p, height = 0.8, width = 0.8, dpi = 600)

  log_x <- log10(df[[x_col]])
  log_y <- log10(df[[y_col]])
  valid_indices <- is.finite(log_x) & is.finite(log_y)
  log_x <- log_x[valid_indices]
  log_y <- log_y[valid_indices]
  pearson_corr <- cor(log_x, log_y, method = "pearson")

  print(paste("Pearson correlation between replicates (log scale):", pearson_corr))
}



python_output_dir <- "oPool_result/assembly_QC/freq"
output_plot_dir <- "graph/assembly_QC"
output_data_dir <- "graph/assembly_QC/data"
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)
if (!dir.exists(output_data_dir)) dir.create(output_data_dir, recursive = TRUE)

python_output_files <- list.files(python_output_dir, pattern = "*.tsv", full.names = TRUE)

for (infile in python_output_files) {
  base_name <- tools::file_path_sans_ext(basename(infile))
  plot_file <- file.path(output_plot_dir, paste0(base_name, "_plot.png"))
  data_file <- file.path(output_data_dir, paste0(base_name, "_plotted_data.tsv"))

  plot_frequencies(
    infile = infile,
    outfile = plot_file,
    data_outfile = data_file,
    textsize = 5,
    filter_or_not = "No"
  )

  message("Plot and data table generated for: ", basename(infile))
}