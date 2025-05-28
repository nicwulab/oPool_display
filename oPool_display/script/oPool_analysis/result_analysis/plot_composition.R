library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(cowplot)
library(scales)

# List all relevant files
files <- list.files(path = "oPool_result/assembly_QC/scFv_count_files/", pattern = "_assembly_freq.tsv$", full.names = TRUE)

# Helper to process each file
process_file <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  
  # Extract the first number in the filename (e.g., 25 from "H3_25_assembly_freq.tsv")
  dataset_id <- as.numeric(str_extract(basename(file), "[0-9]+"))
  
  # If it fails or it's 50, skip
  if (is.na(dataset_id) || dataset_id == 50) return(NULL)
  
  # Identify replicate columns
  rep1_col <- grep("_Rep1_count$", colnames(df), value = TRUE)
  rep2_col <- grep("_Rep2_count$", colnames(df), value = TRUE)
  
  if (length(rep1_col) == 0 || length(rep2_col) == 0) {
    message(paste("Skipping", file, "- missing Rep1 or Rep2 columns"))
    return(NULL)
  }
  
  # Standardize column names
  df <- df %>% rename(Rep1 = all_of(rep1_col), Rep2 = all_of(rep2_col))
  
  # Pivot and classify
  df_long <- df %>%
    select(lev_dist, Rep1, Rep2) %>%
    pivot_longer(cols = c("Rep1", "Rep2"), names_to = "Replicate", values_to = "Count") %>%
    mutate(
      Type = case_when(
        lev_dist == 0 ~ "Native pairing scFv w/o mutations",
        lev_dist > 0 & lev_dist <= 50 ~ "Native pairing scFv w mutations",
        lev_dist > 50 ~ "Non-native pairing scFv"
      ),
      Dataset = dataset_id,
      Sample = paste0(Dataset, " scFvs, ", Replicate)
    )
  return(df_long)
}

# Combine and summarize
all_data <- bind_rows(lapply(files, process_file))

summary_df <- all_data %>%
  group_by(Sample, Type, Dataset) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percentage = Total / sum(Total)) %>%
  ungroup()

# Factor sample levels ordered by increasing Dataset (25 → 200)
summary_df$Sample <- factor(summary_df$Sample, levels = summary_df %>%
                              distinct(Sample, Dataset) %>%
                              arrange(Dataset) %>%
                              pull(Sample))

summary_df$Type <- factor(summary_df$Type, levels = c("Non-native pairing scFv", "Native pairing scFv w mutations", "Native pairing scFv w/o mutations"))

# Plot
p <- ggplot(summary_df, aes(x = Sample, y = Percentage, fill = Type)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black', size = 1) +
  scale_fill_manual(values = c("#E31A1C", "#1F78B4", "#B2DF8A")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.01), expand = c(0, 0)) +
  theme_bw() +
  theme_cowplot(12) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, face = "bold"),
    axis.text.y = element_text(size = 7, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 7),
    plot.title = element_blank(),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    legend.key.size = unit(0.1, 'in'),
    legend.spacing.x = unit(0.03, 'in'),
    legend.title = element_blank(),
    legend.text = element_text(size = 7, face = "bold"),
    legend.position = 'bottom',
    legend.justification = "center"
  ) +
  labs(x = NULL, y = "% of reads")

ggsave("graph/oPool_analysis/QC/assembly_composition.png", p, width = 6, height = 2)