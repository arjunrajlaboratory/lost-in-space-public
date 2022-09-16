#This plots the player pianos for the methods section of figure 2A

library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
library('R.matlab')
input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files'
save_directory = './../../../plots/plots_for_figures/figure2'

data = readMat(file.path(input_directory, 'S_outpar_20210524_vary_off_add_k_for_paper-rep1_5_1_1.mat'), col_names = F)

plot_five_gene_allele_player_piano = function(data,time1,time2,param,gene_to_plot){
  tdata = data[[1]][[param]][[1]]
  tdata = t(tdata)[c(time1:time2),c(1:10)]
  tdata = tdata %>% as_tibble(rownames = 'time')
  tdata$time = as.integer(tdata$time)
  colnames(tdata) = c('time', 'gene_1_1', 'gene_2_1', 'gene_3_1', 'gene_4_1', 'gene_5_1',
                      'gene_1_2', 'gene_2_2', 'gene_3_2', 'gene_4_2', 'gene_5_2')
  
  tdata = tdata %>% pivot_longer(-time, names_to = 'species')
  
  tdata = tdata %>% mutate(binarized_value = value>3)
  tdata = tdata %>% separate(species, into = c(NA, 'gene', 'allele'), sep = '_')
  tdata$gene = as.integer(tdata$gene)
  tdata$allele = as.integer(tdata$allele)
  
  p1=ggplot(tdata %>% filter(gene == gene_to_plot)) + 
    aes(y = factor(allele), x = time, color = factor(allele), alpha = as.integer(binarized_value)) +
    geom_line(size = 5) +
    scale_alpha(range=c(0,1), guide = 'none')+
    theme_classic()
  
  p2 = ggplot(tdata %>% filter(gene == gene_to_plot)) + 
    aes(y = value, x = time, color = factor(allele)) +
    facet_grid(allele~.) +
    geom_line()+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  plots = list(p2,p1)
  return(plots)
  
}


low_OR_param = 301

high_OR_param = 304

gene_of_interest = 2

t1 = 1900
t2 = t1+500




piano_plot_low_OR = plot_five_gene_allele_player_piano(data, t1, t2, low_OR_param, gene_of_interest)
piano_plot_low_OR[[2]]= piano_plot_low_OR[[2]] + scale_y_discrete(limits = rev)

low_OR_final = plot_grid(plotlist=piano_plot_low_OR, nrow=2, align='v', axis = 'l',rel_heights = c(2,1))


ggsave(file.path(save_directory, 'piano_plot_low_OR_example.eps'), low_OR_final, width = 6, height = 4)

piano_plot_high_OR = plot_five_gene_allele_player_piano(data, t1, t2, high_OR_param, gene_of_interest)
piano_plot_high_OR[[2]]= piano_plot_high_OR[[2]]+ scale_y_discrete(limits = rev)
high_OR_final = plot_grid(plotlist=piano_plot_high_OR, nrow=2, align='v', axis = 'l',rel_heights = c(2,1))

ggsave(file.path(save_directory, 'piano_plot_high_OR_example.eps'), high_OR_final, width = 6, height = 4)
