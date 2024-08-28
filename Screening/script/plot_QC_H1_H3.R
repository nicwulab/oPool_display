library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
require(cowplot)
library(dplyr)

plot_mean_exp_score_by_rep <- function(infile, x_col, y_col, graphname, list_of_samples, neg_samples, new_group_samples) {
  df <- read_tsv(infile)
  
  # Filter data to include only specific samples
  df <- df %>%
    filter(Name %in% list_of_samples)
  print(dim(df))
  
  # Separate negative samples and calculate their mean values
  neg_df <- df %>%
    filter(Name %in% neg_samples)
  print(neg_df)
  average_rep1_value <- mean(neg_df[[x_col]], na.rm = TRUE)
  average_rep2_value <- mean(neg_df[[y_col]], na.rm = TRUE)
  print(average_rep1_value)
  print(average_rep2_value)
  
  # Separate new group samples
  new_group_df <- df %>%
    filter(Name %in% new_group_samples)
  print(new_group_df)
  
  # Color coding based on sample group
df$group <- 'Unknown'  # Default group
df$group[df$Name %in% neg_samples] <-  'Negative Control' # Negative control group
df$group[df$Name %in% new_group_samples] <-  'Positive Control' # New group

# Convert to factor to maintain the order of levels
df$group <- factor(df$group, levels = c("Unknown", "Negative Control", "Positive Control"))

textsize <- 7

# Create the plot
p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col), color = group)) +
  geom_point(size = 1.2, pch = 16, alpha = 0.8) +
  geom_vline(xintercept = 5, linetype = "dotted") +
  geom_hline(yintercept = 2.5, linetype = "dotted") +
  scale_x_log10(breaks = 10^seq(-2, 3, by = 1), labels = c(-2, -1, 0, 1, 2, 3), limits = c(0.05, 300)) +  
  scale_y_log10(breaks = 10^seq(-2, 2, by = 1), labels = c(-2, -1, 0, 1, 2), limits = c(0.01, 100)) + 
  scale_color_manual(values = c("grey50", "coral", "green")) +
  theme_cowplot(12) +
  theme(plot.title = element_blank(),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(size = textsize, face = "bold"),
        axis.text = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.1, 'in'),
        legend.spacing.x = unit(0.03, 'in'),
        legend.title = element_blank(),
        legend.text = element_text(size = textsize - 1, face = "bold"),
        legend.position = 'right') +
  labs(x = bquote(bold(paste('Enrichment (replicate 1)'))), y = bquote(bold(paste('Enrichment (replicate 2)'))))

  # Save the plot
  ggsave(graphname, p, height = 2, width = 3, dpi = 3000)

  print(paste("Correlation between replicates:", cor(df[[x_col]], df[[y_col]], method='pearson')))
}


neg_sample_df <- read_tsv("ref_files/neg_abs_list.tsv")

neg_samples <- neg_sample_df$name
print(neg_samples)

new_group_samples <- c('31.a.55_Heavy','AG11-2F01', '042-100809_2F04')

lib1_ref <- read_tsv('ref_files/300lib.tsv')
lib1_list_of_samples <- lib1_ref$Name


infile <- 'result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv'
column_x1 <- 'H1_stem_avg_enrich'
column_y1 <- 'H3_stem_avg_enrich'
outfile <- 'graph/H1_vs_H3.png'

plot_mean_exp_score_by_rep(infile, column_x1, column_y1, outfile, lib1_list_of_samples, neg_samples, new_group_samples)
