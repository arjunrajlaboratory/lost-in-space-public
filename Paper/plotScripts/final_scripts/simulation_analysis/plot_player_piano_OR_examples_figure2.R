#This script plots 4 examples the OR metric for the example 5 node, degree 1 network of
#the 700 parameter set varying add off and k
# Plots were re-colored in Illustrator

library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
library('R.matlab')
input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files'
save_directory = './../../../plots/plots_for_figures/'

data = readMat(file.path(input_directory, 'S_outpar_20210604_vary_off_add_k_for_paper_long-rep1_5_1_1.mat'), col_names = F)

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

save_fun = function(plot_to_save, plot_name){
  plot_to_save[[2]] = plot_to_save[[2]] + scale_y_discrete(limits = rev)
  final_plot = plot_grid(plotlist=plot_to_save, nrow=2, align='v', axis = 'l',rel_heights = c(2,1))
  
  ggsave(file.path(save_directory, plot_name), final_plot, width = 6, height = 4)
}
values_to_plot= tibble::tribble(
  ~param_id, ~start_time, ~end_time, ~gene,
       332L,       6000L,     6500L,    1L,
       334L,       6000L,     6500L,    1L,
       335L,      11000L,    11500L,    1L,
       336L,       6000L,     6500L,    1L
  )

plots=list()

for(irow in 1:nrow(values_to_plot)) {
  plots[[irow]] = plot_five_gene_allele_player_piano(data,
                                                   values_to_plot[[irow,2]],
                                                   values_to_plot[[irow,3]],
                                                   values_to_plot[[irow,1]],
                                                   values_to_plot[[irow,4]])
}


lapply(1:4, function(i)
       save_fun(plots[[i]], paste0('OR_example', i, '.eps'))
       )

