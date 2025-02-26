library(tidyverse)

# Load the input file
data <- read_tsv("oPool_result/enrichment/combined_enrichment.tsv")

# Define antigens and their corresponding custom cutoffs
antigens <- c("H1_stem","H3_stem","SI06_H1", "MI15_H1", "SP16_H3", "QH_H5", "SH_H7", "Phu_FluB", "Lee_FluB")
custom_cutoffs <- list(
  "H1_stem" = c(28),
  "H3_stem" = c(28), 
  "SI06_H1" = c(5, 9),
  "MI15_H1" = c(3, 7),
  "SP16_H3" = c(5),
  "QH_H5"   = c(19, 10),
  "SH_H7"   = c(2.85, 10),
  "Phu_FluB" = c(5, 10),
  "Lee_FluB" = c(5, 11.5)
)

# Initialize an empty list to store results for each antigen
results <- list()

# Process and filter data for each antigen
for (antigen in antigens) {
  full_ha_col_raw <- paste0("Full_HA_", antigen, "_avg_enrich")
  cr9114_col_raw <- paste0("CR9114_compitition_", antigen, "_avg_enrich")
  full_ha_col_norm <- paste0("Full_HA_", antigen, "_avg_enrich_normalized")
  cr9114_col_norm <- paste0("CR9114_compitition_", antigen, "_avg_enrich_normalized")

  # Select relevant columns if they exist
  selected_cols <- c("closest_abs", full_ha_col_raw, cr9114_col_raw, full_ha_col_norm, cr9114_col_norm)
  selected_cols <- selected_cols[selected_cols %in% names(data)]

  # Apply custom cutoff
  cutoffs <- custom_cutoffs[[antigen]]

  # Conditional filtering based on existing columns
  if (full_ha_col_norm %in% names(data) & cr9114_col_norm %in% names(data)) {
    filtered_data <- data %>%
      select(all_of(selected_cols)) %>%
      filter(!!sym(full_ha_col_norm) > cutoffs[1] | !!sym(cr9114_col_norm) > cutoffs[2])
  } else if (full_ha_col_norm %in% names(data)) {
    filtered_data <- data %>%
      select(all_of(selected_cols)) %>%
      filter(!!sym(full_ha_col_norm) > cutoffs[1])
  } else {
    next  # Skip if no relevant columns are available
  }

  # Calculate competition index if both columns exist
  if (full_ha_col_raw %in% names(data) & cr9114_col_raw %in% names(data)) {
    filtered_data <- filtered_data %>%
      mutate(competition_index = log10(!!sym(full_ha_col_raw) / !!sym(cr9114_col_raw))) %>%
      drop_na(competition_index)  # Remove rows with NA values
  }

  # Save the filtered data if there are rows to save
  if (nrow(filtered_data) > 0) {
    write_tsv(filtered_data, paste0("oPool_result/filtered_hits/filtered_antigen_", antigen, "_custom_cutoff.tsv"))
    print(paste("Filtered table saved for antigen:", antigen))
  }
}

print("Cutoff-based filtering completed.")
