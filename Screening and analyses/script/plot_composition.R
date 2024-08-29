library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)
library(cowplot)
library(scales)

# Loading the data
Anno <- read.table('ref_files/sample_name.tsv', header = 1)
TB <- read.csv('result/PacBio/count_result_prescreen.csv')

# Reshaping and transforming the data
TBm <- melt(TB)
TBm$Type = 'Perfect'
TBm$Type[TBm$X %in% c('TooLong', 'Incomplete', 'NoneNative', 'Mutate', 'FailAlign')] <- TBm$X[TBm$X %in% c('TooLong', 'Incomplete', 'NoneNative', 'Mutate', 'FailAlign')]
TBm$Type <- factor(TBm$Type, levels = levels(as.factor(TBm$Type))[c(1:4,6,5)])
TBm$value[TBm$value <0] = 0

TBm$Sample <- Anno$sample_name[match(TBm$variable, paste('Sample', Anno$sample_ID, sep = '_'))]

# Filtering and adjusting the Type levels for plotting
tmp <- TBm[grep("prescreen", TBm$Sample),]
tmp$group <- data.frame(str_split_fixed(tmp$Sample, '_', 2))[[1]]

tmp$Type <- as.character(tmp$Type)
tmp$Type[!tmp$Type %in% c('Perfect', 'Mutate')] <- 'Other'
# Reorder the Type levels to ensure the correct order in the plot
tmp$Type <- factor(tmp$Type, levels = c('Other','Mutate', 'Perfect'))

# Calculate the proportions and percentages again to ensure the ordering is correct
tmp_summary <- tmp %>%
  group_by(group, Type) %>%
  summarise(value_sum = sum(value)) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(total = sum(value_sum),
         proportion = value_sum / total,
         percentage = proportion * 100)

# Plot with reordered bars
ggplot(tmp_summary, aes(x = value_sum, y = group, fill = Type)) + 
  geom_bar(stat = 'identity', position = 'fill', color = 'black', size = 1) + 
  scale_fill_manual(values= c("#E31A1C", "#1F78B4", "#B2DF8A")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = .5, hjust = 1)) +
  theme_cowplot(12) +
  theme(plot.title=element_blank(),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title=element_text(face="bold"),
        axis.text=element_text(face="bold"),
        legend.key.size=unit(0.1,'in'),
        legend.spacing.x=unit(0.03, 'in'),
        legend.title=element_blank(),
        legend.text=element_text(size=6,face="bold"),
        legend.position='right') + 
  coord_cartesian(xlim = c(0,1.01), expand = F, ylim = c(0.4, 2.6)) + 
  scale_x_continuous(labels = scales::percent)

# Save the updated plot
ggsave('graph/opool_composition.png', w = 6, h = 2)