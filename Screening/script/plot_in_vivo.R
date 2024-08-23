#R code
options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(ggplot2)
library(scales)
library(RColorBrewer)
library(qualpalr)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
#install.packages('data.table')
library(data.table)
library(gridExtra)
library(stringr)
require(cowplot)
#install.packages('ggrepel')
library(ggrepel)

shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
  }

plot_ratio_to_Vero <- function(df, graphname, h, w, xlab, ylab){
  textsize <- 7
  mut_levels <- c('WT','S813V','S813K')
  mean <- rowMeans(select(df,rep1,rep2,rep3,rep4,rep5,rep6,rep7), na.rm=TRUE)
  SD   <- apply(select(df,rep1,rep2,rep3,rep4,rep5,rep6,rep7),1,sd,na.rm = FALSE)
  df <- df %>%
          mutate(mean=mean) %>%
          mutate(SD=SD) %>%
          mutate(mut=factor(mut, levels=mut_levels))
  reshape <- melt(data.table(df),id=c('mut')) %>%
               filter(variable %in% c('rep1', 'rep2', 'rep3', 'rep4','rep5','rep6','rep7'))
  palette <- c('black',"#AD07E3", "#0000FF")
  p <- ggplot() +
    geom_bar(data=df,aes(x=mut,y=mean, fill=mut),stat="identity", position=position_dodge(), width=0.7, alpha=0.5) +
    geom_point(data=reshape,aes(x=mut,y=value, color=mut), size=1, alpha=0.8,pch=16, position=position_jitter(w = 0.2, h = 0)) +
    theme_cowplot(12) +
    scale_fill_manual(values=palette) +
    scale_color_manual(values=palette) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          axis.text.x=element_text(angle = 90, hjust = 1,size=textsize, vjust=0.5,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.05, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=textsize,face="bold"),
          legend.justification = "center",
          legend.position='none') +
    labs(x=xlab,y=ylab)
  if (graphname == 'graph/plaque_ratio_VeroTA.png'){
    p <- p + 
           scale_y_continuous(trans=log10_trans(), limits=c(1,100), breaks=c(1,10,100),
                              labels=c(expression(bold('10'^'0')), expression(bold('10'^'1')), expression(bold('10'^'2'))))
    }
  ggsave(graphname, p, height=h, width=w, dpi=600)
  }

plot_invivo_weight_loss <- function(df, graphname, h, w, xlab, ylab){
  textsize <- 7
  conditions <- c('Protection', 'Treatment', 'Control')
  df <- df %>%
          mutate(exp=factor(exp, levels=conditions))
  palette <- c('orange',"purple", "gray30")
  p <- ggplot() +
    geom_errorbar(data=df,aes(x=dpi,ymax=mean+SD,ymin=mean-SD,color=exp,group=exp),linewidth=0.3, width=0.3, alpha=1) +
    geom_line(data=df,aes(x=dpi,y=mean, color=exp),linewidth=0.3, alpha=0.8) +
    geom_point(data=df,aes(x=dpi,y=mean, color=exp), size=1, alpha=1, pch=16, position=position_jitter(w = 0, h = 0)) +
    theme_cowplot(12) +
    scale_fill_manual(values=palette) +
    scale_color_manual(values=palette) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          axis.text.x=element_text(angle = 0, hjust = 0.5,size=textsize, vjust=1,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.05, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=textsize,face="bold"),
          legend.justification = "center",
          legend.position=c(1.15, 0.5),
          plot.margin = unit(c(1, 5, 1, 1), "lines")) +
    labs(x=xlab,y=ylab) +
    scale_y_continuous(limits=c(63,110), breaks=seq(70,110,10),labels=seq(70,110,10)) +
    scale_x_continuous(limits=c(0,14.5), breaks=seq(0,14,2),labels=seq(0,14,2))
  ggsave(graphname, p, height=h, width=w, dpi=600)
  }

library(ggplot2)
library(cowplot)
library(dplyr)

plot_survival <- function(df, graphname, h, w, xlab, ylab){
  textsize <- 7
  conditions <- c('Protection', 'Treatment', 'Control')
  df <- df %>%
          mutate(exp=factor(exp, levels=conditions))
  palette <- c('orange',"purple", "gray30")
  p <- ggplot() +
    geom_step(data=df,aes(x=dpi,y=survive,color=exp),linewidth=0.3, alpha=1) +
    theme_cowplot(12) +
    scale_fill_manual(values=palette) +
    scale_color_manual(values=palette) +
    theme(
      plot.title=element_blank(),
      plot.background = element_rect(fill = "white"),
      axis.title=element_text(size=textsize,face="bold"),
      axis.text=element_text(size=textsize,face="bold"),
      axis.text.x=element_text(angle = 0, hjust = 0.5,size=textsize, vjust=1,face="bold"),
      legend.key.size=unit(0.1,'in'),
      legend.spacing.x=unit(0.05, 'in'),
      legend.title=element_blank(),
      legend.text=element_text(size=textsize,face="bold"),
      legend.justification = "center",
      legend.position=c(1.05, 0.5),
      plot.margin = unit(c(1, 3, 1, 1), "lines")
    ) +
    labs(x=xlab,y=ylab) +
    scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20),labels=seq(0,100,20)) +
    scale_x_continuous(limits=c(0,14.5), breaks=seq(0,14,2),labels=seq(0,14,2))
  
  ggsave(graphname, p, height=h, width=w, dpi=600)
}

# Example usage
# Assuming `df` is your dataframe and it has the columns `dpi`, `survive`, and `exp`
# plot_survival(df, "survival_plot.png", 6, 8, "dpi", "survival (%)")

plot_lung_titer <- function(df, graphname, h, w, xlab, ylab){
  textsize <- 7
  conditions <- c('Protection', 'Treatment', 'Control')
  df <- df %>%
    mutate(exp = factor(exp, levels = conditions)) %>%
    mutate(TCID50 = log10(TCID50))
  
  df_summary <- df %>%
    group_by(exp) %>%
    summarize(mean = mean(TCID50), .groups = 'drop')
  
  dodge_value <- 0.9
  palette <- c('orange', "purple", "gray30")
  
  # Perform t-tests
  t_tests <- df %>%
    group_by(ab) %>%
    do({
      control_protection <- t.test(TCID50 ~ exp, data = filter(., exp %in% c("Control", "Protection")))
      control_treatment <- t.test(TCID50 ~ exp, data = filter(., exp %in% c("Control", "Treatment")))
      protection_treatment <- t.test(TCID50 ~ exp, data = filter(., exp %in% c("Protection", "Treatment")))
      data.frame(
        ab = unique(.$ab),
        control_protection_p = control_protection$p.value,
        control_treatment_p = control_treatment$p.value,
        protection_treatment_p = protection_treatment$p.value
      )
    })
  print(t_tests)
  
  p <- ggplot() +
    geom_bar(data = df_summary, aes(x = exp, y = mean, fill = exp), stat = "identity",
             position = position_dodge(dodge_value), width = dodge_value, alpha = 0.5) +
    geom_point(data = df, aes(x = exp, y = TCID50, color = exp), size = 1, alpha = 0.8, pch = 16,
               position = position_jitterdodge(jitter.width = 0.7, jitter.height = 0, dodge.width = dodge_value)) +
    theme_cowplot(12) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    theme(plot.title = element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title = element_text(size = textsize, face = "bold"),
          axis.text = element_text(size = textsize, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = textsize, vjust = 1, face = "bold"),
          legend.key.size = unit(0.1, 'in'),
          legend.spacing.x = unit(0.05, 'in'),
          legend.title = element_blank(),
          legend.text = element_text(size = textsize, face = "bold"),
          legend.justification = "center",
          legend.position = c(1.4, 0.5),
          plot.margin = unit(c(1, 5, 1, 1), "lines")) +
    labs(x = xlab, y = ylab) +
    scale_y_continuous(trans = shift_trans(2), limits = c(2, 8.5), breaks = seq(2, 8, 2),
                       labels = c(expression(bold('10'^'2')), expression(bold('10'^'4')),
                                  expression(bold('10'^'6')), expression(bold('10'^'8'))))
  
  ggsave(graphname, p, height = h, width = w, dpi = 600)
}

perform_treatment_ttest <- function(df) {
  # Filter the treatment data
  treatment_data <- df %>% filter(exp == "Treatment", ab %in% c("2F01", "16ND92"))
  print(treatment_data)
  
  # Perform the t-test if there are exactly two levels in Ab
  if (length(unique(treatment_data$ab)) == 2) {
    treatment_comparison <- t.test(TCID50 ~ ab, data = treatment_data)
    cat("\nP-value for 16ND92 Treatment vs 2F01 Treatment: ", treatment_comparison$p.value, "\n")
  } else {
    cat("\nError: Treatment data does not contain exactly two levels of 'Ab'.\n")
  }
}

set.seed(5)
df <- read_tsv('data/invivo_weight_loss.tsv')
plot_invivo_weight_loss(filter(df, ab=='2F01'), 'graph/invivo_weight_loss_2F01.png',1.8,3,'dpi',bquote(bold("% of initial weight")))
plot_invivo_weight_loss(filter(df, ab=='16ND92'), 'graph/invivo_weight_loss_16ND92.png',1.8,3,'dpi',bquote(bold("% of initial weight")))

df <- read_tsv('data/invivo_survival.tsv')
plot_survival(filter(df, ab=='2F01'), 'graph/invivo_survival_2F01.png',1.8,2.7,'dpi',bquote(bold("survival (%)")))
plot_survival(filter(df, ab=='16ND92'), 'graph/invivo_survival_16ND92.png',1.8,2.7,'dpi',bquote(bold("survival (%)")))

df <- read_tsv('data/invivo_lung_titer.tsv')
plot_lung_titer(filter(df, ab=='2F01'), 'graph/invitro_lung_titer_Vero_2F01.png',2.2,2.4,'',expression(bold(log['10']~TCID['50']~"/mL")))
plot_lung_titer(filter(df, ab=='16ND92'), 'graph/invitro_lung_titer_Vero_16ND92.png',2.2,2.4,'',expression(bold(log['10']~TCID['50']~"/mL")))
perform_treatment_ttest(df)
