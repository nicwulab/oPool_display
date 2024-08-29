library(ggplot2)
require(cowplot)
library(rlang)
library(dplyr)


plot_filtered_dot <- function(infile, x_col, y_col, outfile, textsize = 10, filter_or_not = "No", values_to_filter = NA){
  # Read the data
  df <- readr::read_tsv(infile, show_col_types = FALSE)

  if (filter_or_not == "Yes"){
    df <- df %>%
      filter(.data[['closest_abs']] %in% values_to_filter)
  }
  
  # Create the dot plot
  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(size=1.5, pch=16, alpha=0.5) +
    scale_x_log10(breaks = 10^seq(-4, 0, by = 1), labels = c(-4,-3,-2,-1,0), limits = c(0.00005, 0.3)) + 
    scale_y_log10(breaks = 10^seq(-4, 0, by = 1), labels = c(-4,-3,-2,-1,0), limits = c(0.00005, 0.3)) + 
    theme_cowplot(12) +
    theme(
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white"),
      axis.title = element_text(size = textsize, family = "Arial", face = "bold" ),
      axis.text = element_text(size = textsize, face = "bold", family = "Arial"),
      legend.key.size = unit(0.1, 'in'),
      legend.spacing.x = unit(0.03, 'in'),
      legend.title = element_blank(),
      legend.text = element_text(size = textsize - 1, face = "bold", family = "Arial"),
      legend.position = 'right'
    ) +
    labs(
      x = bquote(bold(paste('Replicate 1'))),
      y = bquote(bold(paste('Replicate 2')))
    )
  
  ggsave(outfile, p, height = 2.7, width = 2.7, dpi = 3000)
  print(paste("Correlation between replicates:", cor(df[[x_col]], df[[y_col]], method = 'pearson')))
}

# Example usage
infile <- 'data/50Abs_assembly_PCR_strategy_count.tsv'
outfile_50 <- 'graph/pcr_freq_50.png'
plot_filtered_dot(infile, '50Abs_Pool_1+2_N+N_freq', '50Abs_Pool_1+2_NC+N_freq', outfile_50)

outfile_25 <- 'graph/pcr_freq_25.png'
values_to_filter <- c("007_11_4A04", "FISW-n.01", "150055-023_5C05", "007_11_4E01", "019_10_1E06",
                       "429B01", "LAH3V", "315-02-1A05", "315-04-1B12", "014_10_1D03", "014_10_3C04",
                       "AG9-4B03", "K77-2H01", "R95-1B07", "R95-1D11", "009_10_2G06", "013_10_3F02",
                       "039_10_4D04", "051_09_4D03", "220-14-IgG_1E01", "236-14-IgG_1B03",
                       "042-100809_2F04", "T3-1A", "01.a.30_Heavy", "01.o.02_Heavy")

plot_filtered_dot(infile, '50Abs_Pool_1_N+N_freq', '50Abs_Pool_1_NC+N_freq', outfile_25, textsize = 10, filter_or_not = "Yes", values_to_filter)