library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
require(cowplot)
library(dplyr)

plot_mean_exp_score_by_rep <- function(infile, x_col, y_col, graphname, breaks, labels, limits, list_of_samples, neg_samples, new_group_samples) {
  df <- read_tsv(infile)
  
  # Filter data to include only specific samples
  df <- df %>%
    filter(name %in% list_of_samples)

  df <- df %>% 
  filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col))) 
  
  # Separate negative samples and calculate their mean values
  neg_df <- df %>%
    filter(name %in% neg_samples)
  print(neg_df)
  average_rep1_value <- mean(neg_df[[x_col]], na.rm = TRUE)
  average_rep2_value <- mean(neg_df[[y_col]], na.rm = TRUE)
  print(average_rep1_value)
  print(average_rep2_value)
  
  # Separate new group samples
  new_group_df <- df %>%
    filter(name %in% new_group_samples)
  print(new_group_df)
  
# Color coding based on sample group
df$group <- 'Unknown'  # Default group
df$group[df$name %in% neg_samples] <-  'Head Abs Control' # Negative control group
df$group[df$name %in% new_group_samples] <-  'Stem Abs Control' # New group

# Convert to factor to maintain the order of levels
df$group <- factor(df$group, levels = c("Unknown", "Head Abs Control", "Stem Abs Control"))

textsize <- 7
  # Create the plot
  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col),color = group)) +
    geom_point(size = 1, pch = 16, alpha = 0.8) +  
    geom_smooth(method = "lm", color = "red", size = 0.5, se = FALSE, linetype = "dashed") + 
    scale_color_manual(values = c("grey50", "#2171b5", "green")) +# Fit line
    scale_x_log10(breaks = breaks, labels = labels, limits = limits) +  
    scale_y_log10(breaks = breaks, labels = labels, limits = limits) +  
    theme_cowplot(12) +
    theme(
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white"),
      axis.title = element_blank(),
      axis.text = element_text(size = 9, face = "bold"),
      legend.key.size = unit(0.1, 'in'),
      legend.spacing.x = unit(0.03, 'in'),
      legend.title = element_blank(),
      legend.text = element_text(size = textsize - 1, face = "bold"),
      legend.position = 'none'
    ) +
    labs(
      x = bquote(bold(Log[10]~"Enrichment (replicate 1)")),
      y = bquote(bold(Log[10]~"Enrichment (replicate 2)"))
    )

  # Save the plot
  ggsave(graphname, p, height = 1.5, width = 1.5, dpi = 300)

  # Apply log transformation (adding a small constant to avoid log(0) if needed)
  log_x <- log10(df[[x_col]])  # Add 1e-6 to handle zeros if present
  log_y <- log10(df[[y_col]])

  # Calculate the Pearson correlation on log-transformed data
  correlation <- cor(log_x, log_y, method = 'pearson', use = 'complete.obs')
  print(paste("Correlation between replicates:", correlation))
}

# Read negative sample list
neg_sample_df <- read_tsv("ref_files/neg_abs_list.tsv")
neg_samples <- neg_sample_df$name
print(neg_samples)

# Read library sample list
lib1_ref <- read_csv('ref_files/lib_ref.csv')
lib1_list_of_samples <- lib1_ref$Name

# Define new group sample list
new_group_samples <- c('31.a.55', '042-100809-2F04', 'AG11-2F01')
print(new_group_samples)

# Define parameters for multiple plots
plot_params <- list(
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_SI06_H1_Rep1_enrich',
    y_col = 'Full_HA_SI06_H1_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_SI06_H1_screen_correlation.png',
    breaks = 10^seq(-2, 3, by = 1),
    labels = c(-2, -1, 0, 1, 2, 3),
    limits = c(0.002, 300)
  ),
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_MI15_H1_Rep1_enrich',
    y_col = 'Full_HA_MI15_H1_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_MI15_H1_screen_correlation.png',
    breaks = 10^seq(-2, 2, by = 1),
    labels = c(-2,-1, 0, 1, 2),
    limits = c(0.01, 150)
  ),
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_SP16_H3_Rep1_enrich',
    y_col = 'Full_HA_SP16_H3_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_SP16_H3_screen_correlation.png',
    breaks = 10^seq(-3, 3, by = 1),
    labels = c(-3, -2, -1, 0, 1, 2, 3),
    limits = c(0.001, 3000)
  ),
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_QH_H5_Rep1_enrich',
    y_col = 'Full_HA_QH_H5_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_QH_H5_screen_correlation.png',
    breaks = 10^seq(-2, 2, by = 1),
    labels = c(-2,-1, 0, 1, 2),
    limits = c(0.005, 100)
  ),
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_SH_H7_Rep1_enrich',
    y_col = 'Full_HA_SH_H7_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_SH_H7_screen_correlation.png',
    breaks = 10^seq(-2, 2, by = 1),
    labels = c(-2,-1, 0, 1, 2),
    limits = c(0.01, 150)
    ),
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_Phu_FluB_Rep1_enrich',
    y_col = 'Full_HA_Phu_FluB_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_Phu_FluB_screen_correlation.png',
    breaks = 10^seq(-3, 2, by = 1),
    labels = c(-3, -2,-1, 0, 1, 2),
    limits = c(0.001, 150)
    ),
  list(
    infile = 'oPool_result/enrichment/table_s3_3.tsv',
    x_col = 'Full_HA_Lee_FluB_Rep1_enrich',
    y_col = 'Full_HA_Lee_FluB_Rep2_enrich',
    graphname = 'graph/oPool_analysis/QC/full_HA/HA_Lee_FluB_screen_correlation.png',
    breaks = 10^seq(-3, 2, by = 1),
    labels = c(-3,-2,-1, 0, 1, 2),
    limits = c(0.001, 150)
    )
)

# Loop through parameters to generate plots
for (params in plot_params) {
  print(params$graphname)
  plot_mean_exp_score_by_rep(
    infile = params$infile,
    x_col = params$x_col,
    y_col = params$y_col,
    graphname = params$graphname,
    breaks = params$breaks,
    labels = params$labels,
    limits = params$limits,
    lib1_list_of_samples,
    neg_samples,
    new_group_samples
  )
}

