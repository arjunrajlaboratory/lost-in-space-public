library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
gray_palette = c('#D8D5D5', '#898082', '#353131')

#####
#inputs
input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/pre_post_signal/'
input_directory2 = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/'

save_directory = './../../../plots/plots_for_figures/'

#####
#data import
sim_class = read_csv('./../../../extractedData/20210528_for_paper/vary_off_add_k_700params/long_simulation_class-rep1.csv', col_names = F)
colnames(sim_class) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

time_constant_df = read_csv(paste0(input_directory, 'export_time_constant.csv'), col_names = FALSE)
colnames(time_constant_df) = c("pre_rep", "param_id", "species", "time_constant")

time_constant_df$species = as.character(time_constant_df$species)
pre_means = read_csv(paste0(input_directory, 'pre_signal_means.csv'), col_names = FALSE)
colnames(pre_means) = c("pre_rep", "param_id", "mean_1", "mean_2", "mean_3", "mean_4", "mean_5")

post_means = read_csv(paste0(input_directory, 'post_signal_means.csv'), col_names = FALSE)
colnames(post_means) = c("pre_rep", "param_id", "mean_1", "mean_2", "mean_3", "mean_4", "mean_5")

metadata = read_csv(paste0('./../../../extractedData/20210528_for_paper/vary_off_add_k_700params/', 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")

starting_values = read_csv(paste0(input_directory, 'starting_values.csv'), col_names = FALSE)
colnames(starting_values) = c("pre_rep", "param_id", "starting_1", "starting_2", "starting_3", "starting_4", "starting_5")

metadata = metadata %>% add_column(param_id = 1:nrow(metadata))

pre_means = pre_means %>% pivot_longer(mean_1:mean_5, names_to = "species", values_to = 'pre_mean') %>% separate(species, into = c(NA, 'species'), sep = '_')
post_means = post_means %>% pivot_longer(mean_1:mean_5, names_to = "species", values_to = 'post_mean') %>% separate(species, into = c(NA, 'species'), sep = '_')
starting_values = starting_values %>% pivot_longer(starting_1:starting_5, names_to = "species", values_to = 'starting_value') %>% separate(species, into = c(NA, 'species'), sep = '_')

all_signal_data = time_constant_df %>% left_join(pre_means) %>% left_join(post_means) %>% left_join(starting_values) %>% left_join(metadata)
all_signal_data = all_signal_data %>% mutate(dynamic_range = log2(post_mean/pre_mean))
all_signal_data = all_signal_data %>% mutate(time_constant_no_nan = time_constant)
all_signal_data$time_constant_no_nan[is.na(all_signal_data$time_constant_no_nan)] = 0

all_signal_data = all_signal_data %>% left_join(sim_class)

summary_signal_data = all_signal_data %>% group_by(param_id, species) %>% summarize(mean_time = mean(time_constant, na.rm=T),
                                                                                    mean_time_nan_zero = mean(time_constant_no_nan),
                                                                                    mean_pre_mean = mean(pre_mean),
                                                                                    mean_post_mean= mean(post_mean))

summary_signal_data = summary_signal_data %>% mutate(mean_dynamic_range = log2(mean_post_mean/mean_pre_mean))

summary_signal_data = summary_signal_data %>% left_join(metadata) %>% add_column(n_species = 5,
                                                                                 subnet = 1)

summary_signal_data = summary_signal_data %>% left_join(sim_class)
bursts = tibble()
files_to_read = paste0(input_directory2, 'long_bursts-rep', 1, '.csv')
bursts = bind_rows(bursts, files_to_read %>% 
                     lapply(read_csv, col_names = FALSE))
OR = tibble()
files_to_read = paste0(input_directory2, 'long_OR-rep', 1, '.csv')
OR = bind_rows(OR, files_to_read %>% 
                 lapply(read_csv, col_names = FALSE))
colnames(OR) = c("n_species", 'subnet', 'param_id', "rep", 'mean_OR')

OR$mean_OR = as.double(OR$mean_OR)

colnames(bursts) = c("n_species", 'subnet', 'param_id', "rep", "number_bursts", "mean_burst_length", "mean_total_nodes", 'mean_concurrent_nodes', "mean_fraction_on")

metadata = read_csv(paste0(input_directory2, 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")

metadata = metadata %>% add_column(param_id = 1:nrow(metadata))


all_data = left_join(OR, metadata) %>% left_join(bursts) %>% mutate(log_mean_OR = log(mean_OR))

full_data = summary_signal_data %>% left_join(all_data)
saveRDS(full_data, './../../../extractedData/20210528_for_paper/signal_data.rds')
#####
#plots-heatmaps
dynamic_range_plot = ggplot(summary_signal_data %>% filter(k == 110)) +
  aes(x = r_add, y= r_off, label = param_id) +
  geom_tile(data = summary_signal_data %>% filter(k == 110, is_constant_high == 0), aes(fill = mean_dynamic_range)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1, species = 2') +
  scale_fill_viridis(direction = 1, option = "mako", name ='DNR')+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_grid(.~species)+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, 'mean_DNR.svg'), dynamic_range_plot,width = 5, height = 4)

dynamic_range_plot_node2 = ggplot(summary_signal_data %>% filter(k == 110, species ==2)) +
  aes(x = r_add, y= r_off, label = param_id) +
  geom_tile(data = summary_signal_data %>% filter(k == 110, species ==2, is_constant_high == 0), aes(fill = mean_dynamic_range)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1, species = 2') +
  scale_fill_viridis(direction = 1, option = "mako", name ='DNR')+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, 'DNR_node2.svg'), dynamic_range_plot_node2,width = 5, height = 4)


time_constant=ggplot(summary_signal_data %>% filter(k == 110, species ==2)) +
  aes(x = r_add, y= r_off, label = param_id) +
  geom_tile(data = summary_signal_data %>% filter(k == 110, species ==2, is_constant_high == 0), aes(fill = mean_time)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1, species = 2') +
  scale_fill_viridis(direction = 1, option = "mako", name ='time')+
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  facet_grid(.~species)+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, 'mean_time_constant_node2.svg'), time_constant,width = 5, height = 4)

plot_data = full_data %>% filter(species ==2, mean_dynamic_range > 2, is_constant_high==0)
mean_total_nodes_v_time_constant_filt_DNR= ggplot(plot_data) +
  aes(x = mean_total_nodes, y = mean_time_nan_zero, color = factor(k),label = param_id) +
  geom_point(shape=16) +
  geom_line() +
  theme_classic()+
  scale_color_viridis_d(option='plasma', end = 0.95)+
  labs(x = 'mean_total_nodes',
       y = 'time constant',
       title= 'All dynamic range at least 4 fold')
ggsave(paste0(save_directory, 'corr_v_time_const.svg'), mean_total_nodes_v_time_constant_filt_DNR,width = 8, height = 4)
