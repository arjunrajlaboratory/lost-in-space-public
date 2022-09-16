library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
library(datapasta)

#####
#inputs
input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/pre_post_signal/'
input_directory2 = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/'

save_directory = './../../../plots/plots_for_figures/supplemental_figures/'

sim_class = read_csv(paste0(input_directory2, 'long_simulation_class-rep', 1, '.csv'), col_names = F)
colnames(sim_class) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

input_directory_ccr = './../../../extractedData/20210528_for_paper/cross_correlation/'

max_allele_xcorr = tibble()

max_allele_xcorr = read_csv(paste0(input_directory_ccr, 'max_xcorrs_allele_species5-rep1.csv'), col_names = F)
colnames(max_allele_xcorr) = c("n_species", 'subnet', 'param_id', "rep", 'ispecies', "cor", 'lag')

mean_allele_xcorr = max_allele_xcorr %>% group_by(n_species, subnet, param_id, rep) %>% summarize(mean_cor = mean(cor))

max_burst_xcorr = tibble()

max_burst_xcorr = read_csv(paste0(input_directory_ccr, 'max_xcorrs_burst_species5-rep1.csv'), col_names = F)
colnames(max_burst_xcorr) = c("n_species", 'subnet', 'param_id', "rep", 'ispecies_1', 'ispecies_2', "cor", 'lag')

single_burst_xcorr = max_burst_xcorr %>% filter(ispecies_1 == 1 & ispecies_2 == 2)

input_directory_burst_window = './../../../extractedData/20210528_for_paper/burst_window_test/'
burst_window = read_csv(paste0(input_directory_burst_window, 'long_bursts-rep1.csv'), col_names = F)
colnames(burst_window) = c("n_species", 'subnet', 'param_id', "rep", 'burst_window', "number_bursts", "mean_burst_length", "mean_total_window_nodes", 'mean_concurrent_nodes', "mean_fraction_on")

full_data=readRDS('./../../../extractedData/20210528_for_paper/signal_data.rds')
full_data = full_data %>% left_join(mean_allele_xcorr)

full_data = single_burst_xcorr %>% select(c(n_species:rep, cor)) %>% right_join(full_data)
full_data = full_data %>% left_join(sim_class)
full_data = burst_window %>% select(c("n_species", 'subnet', 'param_id', "rep", 'burst_window', 'mean_total_window_nodes')) %>% right_join(full_data)

plot_data = full_data %>% filter(species ==2, mean_dynamic_range > 2, is_constant_high == 0)
ggplot(plot_data) +
  aes(x = cor, y = mean_time_nan_zero, color = factor(k),label = param_id) +
  geom_point()+
  geom_line() +
  theme_classic()+
  scale_color_viridis_d(option='plasma', end = 0.95)+
  labs(x = 'gene-gene correlation',
       y = 'time constant',
       title= 'All dynamic range at least 4 fold')
ggsave(paste0(save_directory, 'gene-gene_cor_v_time_const.svg'),width = 8, height = 4)

ggplot(plot_data) +
  aes(x = mean_total_window_nodes, y = mean_time_nan_zero, color = factor(k),label = param_id) +
  geom_point()+
  geom_line() +
  theme_classic()+
  facet_grid(burst_window~.)+
  scale_color_viridis_d(option='plasma', end = 0.95)+
  labs(x = 'gene-gene correlation',
       y = 'time constant',
       title= 'All dynamic range at least 4 fold')
ggsave(paste0(save_directory, 'burst_window_v_time_const.svg'),width = 8, height = 16)
