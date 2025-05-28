library(ggplot2)
library(readxl)
library(cowplot)
library(dplyr)
library(scales) 

output_dir <- "graph/validation"
if (!dir.exists(output_dir)) dir.create(output_dir)

df <- read_excel("oPool_result/validated_abs_freq.xlsx")

df <- df %>%
  filter(Full_HA_input_Rep1_freq > 0, Full_HA_input_Rep2_freq > 0)

df <- df %>%
  filter(Full_HA_input_Rep1_freq > 0, Full_HA_input_Rep2_freq > 0) %>%
  mutate(
    log10_Rep1 = log10(Full_HA_input_Rep1_freq),
    log10_Rep2 = log10(Full_HA_input_Rep2_freq)
  )

p <- ggplot(df, aes(x = log10_Rep1, y = log10_Rep2)) +
  geom_point(shape = 21, size = 1.8, stroke = 0.2, fill = "#1f77b4", color = "black", alpha = 0.8) +
  scale_x_continuous(
    limits = c(-5, -1),
    breaks = seq(-5, -1, 1),
    labels = seq(-5, -1, 1)
  ) +
  scale_y_continuous(
    limits = c(-5, -1),
    breaks = seq(-5, -1, 1),
    labels = seq(-5, -1, 1)
  ) +
  labs(
    x = expression(log[10]~scFv~frequency~"(Rep 1)"),
    y = expression(log[10]~scFv~frequency~"(Rep 2)")
  ) +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "Full_HA_input_rep_correlation.png"),
  plot = p,
  width = 2, height = 2, dpi = 600, bg = "white"
)