#This script visualizes the bursting and OR metrics for the all nodes/subnets of
#the 700 param set varying add off and k

library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
gray_palette = c('#D8D5D5', '#898082', '#353131')

input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/'

save_directory = './../../../plots/plots_for_figures/supplemental_figures/'
#Read in data and add columns (no columns in MATLAB)
sim_classification = read_csv(paste0(input_directory, 'long_simulation_class-rep', 1, '.csv'), col_names = F)
colnames(sim_classification) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

bursts = tibble()
files_to_read = paste0(input_directory, 'long_bursts-rep', 1, '.csv')
bursts = bind_rows(bursts, files_to_read %>% 
                     lapply(read_csv, col_names = FALSE))
OR = tibble()
files_to_read = paste0(input_directory, 'long_OR-rep', 1, '.csv')
OR = bind_rows(OR, files_to_read %>% 
                 lapply(read_csv, col_names = FALSE))

colnames(OR) = c("n_species", 'subnet', 'param_id', "rep", 'mean_OR')

colnames(bursts) = c("n_species", 'subnet', 'param_id', "rep", "number_bursts", "mean_burst_length", "mean_total_nodes", 'mean_concurrent_nodes', "mean_fraction_on")

metadata = read_csv(paste0(input_directory, 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")

metadata = metadata %>% add_column(param_id = 1:nrow(metadata))


all_data = left_join(OR, metadata) %>% left_join(bursts) %>% left_join(sim_classification) %>% mutate(log_mean_OR = log(mean_OR))


all_k_values_OR=ggplot(all_data %>% filter(n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(n_species ==5, subnet == 1, is_constant_high==0), aes(fill = log_mean_OR)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_wrap(.~k, labeller = 'label_both')+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, "all_k_values_OR.svg"), all_k_values_OR, height = 10, width = 10)

all_k_values_mean_total_nodes=ggplot(all_data %>% filter(n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(n_species ==5, subnet == 1, is_constant_high==0), aes(fill = mean_total_nodes)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  theme_bw()+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  facet_wrap(.~k, labeller = 'label_both')+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, "all_k_values_mean_total_nodes.svg"), all_k_values_mean_total_nodes, height = 10, width = 10)

all_k_values_number_bursts=ggplot(all_data %>% filter(n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(n_species ==5, subnet == 1, is_constant_high==0), aes(fill = number_bursts)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_wrap(.~k)+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))

all_k_values_mean_burst_length=ggplot(all_data %>% filter(n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(n_species ==5, subnet == 1, is_constant_high==0), aes(fill = log(mean_burst_length))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_wrap(.~k)+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))