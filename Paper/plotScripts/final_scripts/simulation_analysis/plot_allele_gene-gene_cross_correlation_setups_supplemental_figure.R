library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
gray_palette = c('#D8D5D5', '#898082', '#353131')

input_directory = './../../../extractedData/20210528_for_paper/cross_correlation/'
save_directory = './../../../plots/plots_for_figures/supplemental_figures/'

max_allele_xcorr = tibble()
max_allele_xcorr = read_csv(paste0(input_directory, 'max_xcorrs_allele_species5-rep1.csv'), col_names = F)
colnames(max_allele_xcorr) = c("n_species", 'subnet', 'param_id', "rep", 'ispecies', "allele_cor", 'lag')

mean_allele_xcorr = max_allele_xcorr %>% group_by(n_species, subnet, param_id, rep) %>% summarize(mean_allele_cor = mean(allele_cor))

max_burst_xcorr = tibble()
max_burst_xcorr = read_csv(paste0(input_directory, 'max_xcorrs_burst_species5-rep1.csv'), col_names = F)
colnames(max_burst_xcorr) = c("n_species", 'subnet', 'param_id', "rep", 'ispecies_1', 'ispecies_2', "gene_gene_cor", 'lag')

single_burst_xcorr = max_burst_xcorr %>% filter(ispecies_1 == 1 & ispecies_2 == 2)

metadata = read_csv(paste0(input_directory, 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")

sim_class = read_csv('./../../../extractedData/20210528_for_paper/vary_off_add_k_700params/long_simulation_class-rep1.csv', col_names = F)
colnames(sim_class) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

metadata = metadata %>% add_column(param_id = 1:nrow(metadata))
max_allele_xcorr = max_allele_xcorr %>% left_join(metadata)
mean_allele_xcorr = mean_allele_xcorr %>% left_join(metadata)
single_burst_xcorr = single_burst_xcorr %>% left_join(metadata)
max_allele_xcorr = max_allele_xcorr %>% left_join(sim_class)
mean_allele_xcorr = mean_allele_xcorr %>% left_join(sim_class)
single_burst_xcorr = single_burst_xcorr %>% left_join(sim_class)

max_burst_xcorr = max_burst_xcorr %>% left_join(metadata) %>% left_join(sim_class)

all_data = mean_allele_xcorr %>% left_join(single_burst_xcorr)

ggplot(max_allele_xcorr %>% filter(k==110, is_constant_high == 0, subnet %in% c(1,5,6,10))) +
  aes(x=lag, y = allele_cor) +
  theme_bw()+
  facet_grid(.~subnet)+
  geom_point(shape = 16)
ggsave(paste0(save_directory, 'allele_cor_lags.svg'), width = 10, height = 4)

ggplot(mean_allele_xcorr %>% filter(k==110, subnet %in% c(1,5,6,10))) +
  aes(x = r_add, y = r_off) +
  geom_tile(data=mean_allele_xcorr %>% filter(k==110, is_constant_high == 0, subnet %in% c(1,5,6,10)), aes(fill = mean_allele_cor))+
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis(name= 'mean_allele_cor', direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_grid(.~subnet)+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, 'allele_cor_heatmaps.svg'), width = 10, height = 4)




ggplot(single_burst_xcorr %>% filter(k==110, is_constant_high == 0, subnet %in% c(1,5,6,10))) +
  aes(x=lag, y = gene_gene_cor) +
  theme_bw()+
  facet_grid(.~subnet)+
  geom_point(shape = 16)
ggsave(paste0(save_directory, 'gene-gene_cor_lags.svg'), width = 10, height = 4)
  

ggplot(single_burst_xcorr %>% filter(k==110, subnet %in% c(1,5,6,10))) +
  aes(x = r_add, y = r_off, fill = gene_gene_cor) +
  geom_tile(data=single_burst_xcorr %>% filter(k==110, is_constant_high == 0, subnet %in% c(1,5,6,10)), aes(fill = gene_gene_cor))+
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis(name= 'direct_gene_gene_xcor', direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_grid(.~subnet)+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, 'gene-gene_heatmaps.svg'), width = 10, height = 4)

direct_connection_summary = max_burst_xcorr %>% filter(k==110, ispecies_1==1, is_constant_high==0, subnet==1)

direct_connection_summary = direct_connection_summary %>% mutate(direct_connection = ifelse(ispecies_2 == 2 | ispecies_2 == 5, T, F))

ggplot(direct_connection_summary)+
  aes(x=direct_connection, y = gene_gene_cor, fill = direct_connection) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic()
ggsave(paste0(save_directory, 'gene-gene_direct_connection.svg'), width = 5, height = 4)
