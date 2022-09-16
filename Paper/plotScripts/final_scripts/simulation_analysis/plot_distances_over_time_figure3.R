library('R.matlab')
library(tidyverse)
library(cowplot)
squared_L2_from_median= readMat.default('./../../../extractedData/20210528_for_paper/vary_off_add_k_700params/mat_files/squared_L2_from_median_20210524_vary_off_add_k_for_paper-rep1_3_1.mat')
save_directory = './../../../plots/plots_for_figures/'


plot_three_distances = function(input_data,time1,time2,param){
  
  raw_data = input_data[[1]][[param]][[1]]
  
  raw_data = t(raw_data)[c(time1:time2),c(1:4)]
  raw_data = raw_data %>% as_tibble(rownames = 'time')
  raw_data$time = as.integer(raw_data$time)
  colnames(raw_data) = c('time', 'gene_1', 'gene_2', 'gene_3', 'total')
  
  temp = raw_data %>% select(-c(total)) %>% pivot_longer(cols=gene_1:gene_3) %>% group_by(time) %>% filter(value == max(value)) %>% arrange(time)
  temp = temp %>% add_count(time) %>% mutate(max_gene = ifelse(n==1, name, 'none'))
  
  full_data = raw_data %>% select(time, total)
  
  full_data = temp %>% select(time, max_gene) %>% distinct(.) %>% right_join(full_data)
  full_data$max_gene = factor(full_data$max_gene, levels = c('gene_1', 'gene_2', 'gene_3', 'none'))
  
  p1=ggplot(full_data) +
    aes(x=time, y= sqrt(total)) +
    geom_step(aes(group=1, color = time)) +
    labs(y= 'euclidian distance from median') +
    scale_color_viridis(limits = c(0,120), option = 'magma')+
    theme_classic()+
    ylim(c(0,200))+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  return(p1)
}

idx = c(341, 302, 304)
t1=1280
t2=1396


plots = lapply(idx, function(x) plot_three_distances(squared_L2_from_median, t1, t2, x))
plots[[3]] + ylim(0,200)
plot_grid(plotlist =  plots, nrow=1)
for (i in 1:length(idx)) {
  ggsave(paste0(save_directory, 'distances_param', idx[i], '.svg'), plots[[i]], width = 6, height = 4)
}
