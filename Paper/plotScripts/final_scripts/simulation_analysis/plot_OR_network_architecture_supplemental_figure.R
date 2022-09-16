#This script visualizes the bursting and OR metrics for the all nodes/subnets of
#the 700 param set varying add off and k

library(tidyverse)
library(plotly)
library(cowplot)
library(viridis)
library(ggnewscale)
library(datapasta)
library(svglite)

gray_palette = c('#D8D5D5', '#898082', '#353131')

#####
#In/out directories
input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/'

save_directory = './../../../plots/plots_for_figures/supplemental_figures/'
#####
#Read in data and add columns (no columns in MATLAB)
sim_classification = read_csv(paste0(input_directory, 'long_simulation_class-rep', 1, '.csv'), col_names = F)
colnames(sim_classification) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

network_metadata = read_csv(paste0('./../../../extractedData/20210528_for_paper/', 'topology_metadata.csv'), col_names = F)
colnames(network_metadata) = c('n_species', 'subnet', 'connectivity', 'auto_regulation', 'characteristic_distance')

network_metadata$auto_regulation = as.logical(network_metadata$auto_regulation)

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
#####
#combine data
all_data = left_join(OR, metadata) %>% left_join(bursts) %>% left_join(network_metadata) %>% left_join(sim_classification) %>% mutate(log_mean_OR = log(mean_OR))


#####
#Netsize plots

facet_netsize=ggplot(all_data %>% filter(subnet == 1, k==110)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(subnet == 1, k==110, is_constant_high == 0), aes(fill = log(mean_OR)))+
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))+
  facet_grid(n_species~.)
ggsave(paste0(save_directory, 'facet_netsize.svg'), facet_netsize,width = 4, height = 6)

facet_netsize_25=ggplot(all_data %>% filter(subnet == 1, k==110, n_species %in% c(2,5))) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(subnet == 1, k==110, is_constant_high == 0, n_species %in% c(2,5)), aes(fill = log(mean_OR)))+
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))+
  facet_grid(.~n_species)
ggsave(paste0(save_directory, 'facet_netsize_25.svg'), facet_netsize_25,width = 6, height = 3)

#####
#Subnet plots

facet_subnet=ggplot(all_data %>% filter(n_species ==5, k ==110, subnet < 6)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(n_species ==5, k ==110, subnet < 6, is_constant_high==0), aes(fill = log(mean_OR))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(limits = c(-.75, 4), direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))+
  facet_grid(subnet~.)

facet_subnet_self_loops=ggplot(all_data %>% filter(n_species ==5, k ==110, subnet > 5)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data =all_data %>% filter(n_species ==5, k ==110, subnet > 5, is_constant_high==0), aes(fill=log(mean_OR))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(limits = c(-.75, 4), direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))+
  facet_grid(subnet~.)

ggsave(paste0(save_directory, 'facet_subnet.svg'), facet_subnet,width = 4, height = 6)
ggsave(paste0(save_directory, 'facet_subnet_self_loops.svg'), facet_subnet_self_loops,width = 4, height = 6)


facet_subnet_select=ggplot(all_data %>% filter(n_species ==5, k ==110, subnet %in% c(1,5,6,10))) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(n_species ==5, k ==110, subnet %in% c(1,5,6,10), is_constant_high==0), aes(fill=log(mean_OR))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))+
  facet_wrap(.~subnet)
ggsave(paste0(save_directory, 'facet_subnet_select.svg'), facet_subnet_select,width = 10, height = 6)


#####
# Single parameter set changing over n_species and connectivity


param355_over_n_species= ggplot(all_data %>% filter(param_id %in% c(355), subnet ==1)) +
  aes(x = n_species, y = log_mean_OR) +
  geom_point()+
  geom_line() +
  theme_classic()+
  ylim(0,2)
ggsave(paste0(save_directory, 'param355_over_n_species.svg'), param355_over_n_species, width = 6, height = 4)

param355_over_connectivity=ggplot(all_data %>% filter(param_id %in% c(355), n_species ==5)) +
  aes(x = connectivity, y = log_mean_OR, color = auto_regulation) +
  geom_point()+
  scale_color_manual(values = c("black", "red"))+
  theme_classic() +
  ylim(0,2) +
  geom_line()
ggsave(paste0(save_directory, 'param355_over_connectivity.svg'), param355_over_connectivity, width = 10, height = 4)
