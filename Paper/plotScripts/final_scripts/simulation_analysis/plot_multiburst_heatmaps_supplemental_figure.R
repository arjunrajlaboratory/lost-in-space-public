#This script visualizes the mutli-burst metrics for supplemental figure on multi-burst
# definition

library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
gray_palette = c('#D8D5D5', '#898082', '#353131')
sim_length = 199900

input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/'

save_directory = './../../../plots/plots_for_figures/supplemental_figures/'
#Read in data and add columns (no columns in MATLAB)

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
sim_class = read_csv(paste0(input_directory, 'long_simulation_class-rep1.csv'), col_names = F)
colnames(sim_class) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')


all_data = left_join(OR, metadata) %>% left_join(bursts) %>% left_join(sim_class) %>% mutate(log_mean_OR = log(mean_OR))
# 
# #Annotate technical outliers for OR
# all_data = all_data %>%
#   mutate(filtered_log_OR = ifelse(log_mean_OR < 0 | is.infinite(log_mean_OR), NA, log_mean_OR)) # Only keep ORs that are not below 1 and are finite
# 
# all_data = all_data %>%
#   mutate(filtered_out_OR = ifelse(log_mean_OR < 0 & is.finite(log_mean_OR), 'sub-sampled', # OR < 1 implies anticorrelation, which is impossible and is therefore undersampled
#                                   ifelse(log_mean_OR < 0 & is.infinite(log_mean_OR), 'zero', # OR = 0 is special case that implies a strong correlation in our data
#                                          ifelse(is.na(log_mean_OR), NA, 'number'))))
# #Annotate technical outliers for burst number
# 
# all_data = all_data %>%
#   mutate(filtered_out_number_bursts = ifelse(number_bursts == 1, 'one', NA))
# 
# all_data = all_data %>%
#   mutate(filtered_number_bursts = ifelse(number_bursts == 1, NA, number_bursts))
# 
# #Annotate technical outliers for burst length
# full_sim_cutoff = max(all_data$mean_burst_length, na.rm=TRUE) - max(all_data$mean_burst_length, na.rm=TRUE) * 0.01
# 
# 
# all_data = all_data %>%
#   mutate(filtered_out_burst_length = ifelse(mean_burst_length >= full_sim_cutoff, 'full-simulation', NA))
# 
# all_data = all_data %>%
#   mutate(filtered_burst_length = ifelse(mean_burst_length >= full_sim_cutoff, NA, mean_burst_length))
# 


example_number_bursts=ggplot(all_data %>% filter(k == 110, n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(k == 110, n_species ==5, subnet == 1, is_constant_high==0), aes(fill = number_bursts)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(name= 'normalized number', direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, "example_number_bursts.svg"), example_number_bursts)

example_burst_length=ggplot(all_data %>% filter(k == 110, n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(k == 110, n_species ==5, subnet == 1, is_constant_high==0), aes(fill = log(mean_burst_length))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, "example_burst_length.svg"), example_burst_length)

example_number_nodes=ggplot(all_data %>% filter(k == 110, n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off, label = paste0(param_id, '_', round(mean_OR,2))) +
  geom_tile(data = all_data %>% filter(k == 110, n_species ==5, subnet == 1, is_constant_high==0), aes(fill = mean_total_nodes)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))

ggsave(paste0(save_directory, "example_number_nodes.svg"), example_number_nodes)
