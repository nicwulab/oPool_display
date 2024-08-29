library(ggplot2)
library(cowplot)
library(rlang)

plot_filtered_dot <- function(infile, x_col, y_col, outfile, textsize = 8){
  # Read the data
  df <- readr::read_tsv(infile, show_col_types = FALSE)
  
  # Create the dot plot
  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(size=1.2, pch=16, alpha=0.5) +
    scale_x_log10(breaks = 10^seq(-5, 0, by = 1), labels = c(-5,-4,-3,-2,-1,0), limits = c(0.000005, 0.3)) + 
    scale_y_log10(breaks = 10^seq(-5, 0, by = 1), labels = c(-5,-4,-3,-2,-1,0), limits = c(0.000005, 0.3)) + 
    theme_cowplot(12) +
    theme(
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white"),
      axis.title = element_text(size = textsize, face = "bold"),
      axis.text = element_text(size = textsize, face = "bold"),
      legend.key.size = unit(0.1, 'in'),
      legend.spacing.x = unit(0.03, 'in'),
      legend.title = element_blank(),
      legend.text = element_text(size = textsize - 1, face = "bold"),
      legend.position = 'right'
    ) +
    labs(
      x = bquote(bold(paste('replicate 1'))),
      y = bquote(bold(paste('replicate 2')))
    )
  
  ggsave(outfile, p, height = 2, width = 2, dpi = 3000)
  print(paste("Correlation between replicates:", cor(df[[x_col]], df[[y_col]], method = 'pearson')))
}



# Example usage
infile <- 'result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv'
outfile <- 'graph/input_QC.png'

plot_filtered_dot(infile, 'Rep1_prescreen_lib_freq', 'Rep2_prescreen_lib_freq', outfile)
