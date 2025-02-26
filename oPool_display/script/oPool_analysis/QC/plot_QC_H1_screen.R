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
p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col), color = group)) +
  geom_point(size = 1.2, pch = 16, alpha = 0.8) +
  geom_vline(xintercept = average_rep1_value, linetype = "dotted") +
  geom_hline(yintercept = average_rep2_value, linetype = "dotted") +
  scale_x_log10(breaks = 10^seq(-2, 3, by = 1), labels = c(-2, -1, 0, 1, 2, 3), limits = c(0.01, 300)) +  
  scale_y_log10(breaks = 10^seq(-2, 3, by = 1), labels = c(-2, -1, 0, 1, 2, 3), limits = c(0.01, 300)) +  
  scale_color_manual(values = c("grey50", "#2171b5", "green")) +
  geom_smooth(method = "lm", color = "red", size = 0.5, se = FALSE, linetype = "dashed") +
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
log_x <- log10(df[[x_col]]) 
log_y <- log10(df[[y_col]])
print(paste("Correlation between replicates:", cor(log_x, log_y, method = 'pearson')))
}

# Read negative sample list
neg_sample_df <- read_tsv("ref_files/neg_abs_list.tsv")
neg_samples <- neg_sample_df$name
print(neg_samples)

# Read library sample list
lib1_ref <- read_csv('ref_files/300lib_clean.csv')
lib1_list_of_samples <- lib1_ref$Name

# Define new group sample list
new_group_samples <- c('31.a.55', 'AG11-2F01')
print(new_group_samples)

# Define file and column names
infile <- 'oPool_result/enrichment/table_s3_3.tsv'
column_x1 <- 'Full_HA_H1_stem_Rep1_enrich'
column_y1 <- 'Full_HA_H1_stem_Rep2_enrich'
outfile <- 'graph/oPool_analysis/QC/full_HA/H1_screen_correlation.png'

# Call the plotting function
plot_mean_exp_score_by_rep(infile, column_x1, column_y1, outfile, lib1_list_of_samples, neg_samples, new_group_samples)
