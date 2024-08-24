library(ggplot2)
library(dplyr)
library(reshape)
library(scales)
library(cowplot)
library(gridExtra)
library(grid)
pdf(NULL)

plotsensorgram <- function(data_sign, data_fit, sample, ymax, ymin){
  textsize <- 7
  p <- ggplot() +
    geom_line(data=data_sign, aes(x=Time, y=Signal, group=SampleID), size=0.8, color='blue') +
    geom_line(data=data_fit, aes(x=Time, y=Signal, group=SampleID), size=0.5, color='#CC6666') +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize, face="bold", hjust=0.5),
          axis.title.y=element_text(size=textsize, face="bold"),
          axis.title.x=element_text(size=textsize, face="bold"),
          axis.text.y=element_text(size=textsize, face="bold", angle=0, hjust=0.5, vjust=0.5),
          axis.text.x=element_text(size=textsize, face="bold", angle=0, hjust=0.5, vjust=0.5), # Horizontal x-axis labels
          legend.title=element_blank(),
          legend.key.size = unit(0.5, 'lines'),
          legend.text=element_text(size=textsize, face="bold"),
          legend.position='none',
          panel.background = element_rect(fill = "white"),        # Set the panel background to white
          plot.background = element_rect(fill = "white"),         # Set the plot background to white
          panel.grid.major = element_blank(),                    # Remove major grid lines
          panel.grid.minor = element_blank()) +                  # Remove minor grid lines
    ylim(ymin, ymax) +
    xlab(bquote(bold(Time~'(s)'))) +
    ylab(bquote(bold(Response~'(nm)'))) +
    ggtitle(sample)
  return(p)
}

adjusttime <- function(t){
  t$Time <- t$Time - min(t$Time)
  return(t)
}

adjustsignal <- function(t, factor){
  t$Signal <- t$Signal * factor
  return(t)
}

plot_wrapper <- function(sample, ymax = 1, ymin = -0.1, output_file){
  folder <- "result"
  fit <- adjusttime(read.table(paste(folder, '/', sample, '_fit_All.compile', sep=''), header=1))
  sign <- adjusttime(read.table(paste(folder, '/', sample, '_sign_All.compile', sep=''), header=1))
  p <- plotsensorgram(sign, fit, sample, ymax, ymin)
  
  # Save the plot to a file
  ggsave(filename=output_file, plot=p, height=2, width=2, dpi=600)
}

# Define your targets and binders
targets <- c('H3-stem','H1-stem')
H3_stem_binders <- c('01.h.02-scFv', '2F02-scFv', '2F04-scFv', '31.a.55-scFv', '56.j.01-scFv', 'AG2-G02-scFv', '01.ad.01-scFv',
  '3C06-scFv','3F02-scFv','01.h.02-Fab', '2F02-Fab', '2F04-Fab', '31.a.55-Fab', '56.j.01-Fab', '01.ad.01-Fab', '3C06-Fab', '3F02-Fab')
H1_stem_binders <- c('16.ND.92-scFv', '2F01-scFv', '31.a.55-scFv','3C06-scFv','3F02-scFv', 
  '16.ND.92-Fab', '2F01-Fab', '31.a.55-Fab','3C06-Fab', '3F02-Fab')

# Store binders and corresponding ymax values in lists
binders_list <- list(H3_stem_binders, H1_stem_binders)
ymax_values_list <- list(c(1, 0.5, 0.5, 0.5, 1.5, 0.5, 1.5, 0.5, 0.5, 1, 2, 2.5, 0.5, 2.5, 2.5, 0.5, 0.5), 
  c(0.5, 1.5, 1.5, 0.5, 0.5, 2, 2, 3, 0.5, 0.5))  # H3 and H1 ymax values

# Iterate over each target and its corresponding binders and ymax values
for (i in 1:length(targets)) {
  target <- targets[i]
  binders <- binders_list[[i]]
  ymax_values <- ymax_values_list[[i]]
  
  for (j in 1:length(binders)) {
    binder <- binders[j]
    ymax <- ymax_values[j]
    
    sample <- paste(binder, '_', target, sep='')
    output_file <- paste('graph/BLI_sensorgram/', binder, '_', target, '_Sensorgram.png', sep='')
    
    # Generate and save the plot
    plot_wrapper(sample, ymax = ymax, ymin = -0.1, output_file = output_file)
  }
}