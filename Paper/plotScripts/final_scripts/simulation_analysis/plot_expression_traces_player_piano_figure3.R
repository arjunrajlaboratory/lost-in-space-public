library('R.matlab')
library(tidyverse)
library(cowplot)

rawMat = readMat.default('./../../../extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files/S_outpar_20210524_vary_off_add_k_for_paper-rep1_3_1_1_summed.mat')
save_directory = './../../../plots/plots_for_figures/'

plot_three_gene_player_piano = function(data,time1,time2,param){
  tdata = data[[1]][[param]][[1]]
  tdata = t(tdata)[c(time1:time2),c(1:3)]
  tdata = tdata %>% as_tibble(rownames = 'time')
  tdata$time = as.integer(tdata$time)
  colnames(tdata) = c('time', 'gene_1', 'gene_2', 'gene_3')
  
  tdata = tdata %>% pivot_longer(-time, names_to = 'species')
  
  tdata = tdata %>% mutate(binarized_value = value>3)
  tdata = tdata %>% separate(species, into = c(NA, 'gene'), sep = '_')
  tdata$gene = as.integer(tdata$gene)
  if(length(unique(tdata$binarized_value)) == 2){
  p1=ggplot(tdata) + 
    aes(y = factor(gene), x = time, color = factor(gene), alpha = as.integer(binarized_value)) +
    geom_line(size = 5) +
    scale_alpha(range=c(0,1), guide = 'none')+
    theme_classic()
  } else{
    p1=ggplot(tdata) + 
      aes(y = factor(gene), x = time, color = factor(gene)) +
      geom_line(size = 5) +
      theme_classic()
  }
  p2 = ggplot(tdata) + 
    aes(y = value, x = time, color = factor(gene)) +
    geom_step()+
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
idx = c(341, 302, 304)
t1=1280
t2=1396

plots= lapply(idx, function(iter) plot_three_gene_player_piano(rawMat, t1, t2, iter))

lapply(1:3, function(i)
  save_fun(plots[[i]], paste0('3dplot_example_', i, '.svg'))
)
