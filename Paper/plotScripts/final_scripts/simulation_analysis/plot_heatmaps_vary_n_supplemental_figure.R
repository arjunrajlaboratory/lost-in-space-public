library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
library(RColorBrewer)
gray_palette = c('#D8D5D5', '#898082', '#353131')

input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_n_params/'

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

sim_class = read_csv(paste0(input_directory, 'simulation_class-rep1.csv'), col_names = F)
colnames(sim_class) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

metadata = read_csv(paste0(input_directory, 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")

metadata = metadata %>% add_column(param_id = 1:nrow(metadata))

all_data = OR %>% left_join(bursts) %>% left_join(metadata) %>% left_join(sim_class)
all_data = all_data %>% mutate(log_mean_OR = log(mean_OR))
# all_data = all_data %>%
#   mutate(filtered_log_OR = ifelse(log_mean_OR < 0 | is.infinite(log_mean_OR), NA, log_mean_OR)) # Only keep ORs that are not below 1 and are finite
# 
# all_data = all_data %>%
#   mutate(filtered_out_OR = ifelse(log_mean_OR < 0 & is.finite(log_mean_OR), 'sub-sampled', # OR < 1 implies anticorrelation, which is impossible and is therefore undersampled
#                                   ifelse(log_mean_OR < 0 & is.infinite(log_mean_OR), 'zero', # OR = 0 is special case that implies a strong correlation in our data
#                                          ifelse(is.na(log_mean_OR), NA, 'number'))))


facet_n_same_scale=ggplot(all_data) +
  aes(x = r_add, y = r_off) +
  geom_tile(data=all_data %>% filter(is_constant_high==0), aes(fill = log_mean_OR)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  facet_grid(.~n, labeller = 'label_both') +
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, "facet_n_same_scale_OR.svg"), facet_n_same_scale, height = 4, width = 20)


n_plots = lapply(c(0.2, 0.5, 1, 2, 4), function(x) ggplot(all_data %>% filter(n==x))+
                                                            aes(x = r_add, y = r_off) +
                                                            geom_tile(data=all_data %>% filter(n==x, is_constant_high==0), aes(fill = log_mean_OR)) +
                                                            scale_x_log10() +
                                                            scale_y_log10() +
                                                            labs(title = paste0('N = ', x)) +
                                                            scale_fill_viridis(direction = 1, option = "mako")+
                                                            new_scale_fill()+
                                                            geom_tile(aes(fill=as.factor(is_constant_high)))+
                                                            scale_fill_manual(values = c('transparent', gray_palette[3]))+
                                                            theme_bw()+
                                                            theme(panel.border = element_blank(), 
                                                                  panel.grid.major = element_blank(), 
                                                                  axis.text = element_text(size = 12), 
                                                                  plot.title = element_text(size = 16, face = 'bold'), 
                                                                  axis.title = element_text(face = 'bold')))
plot_grid(plotlist = n_plots, nrow = 1)
ggsave(paste0(save_directory, "facet_n_individual_scale_OR.svg"), facet_n_individual_scale, height = 4, width = 30)
facet_n_same_scale_total_nodes=ggplot(all_data) +
  aes(x = r_add, y = r_off) +
  geom_tile(data=all_data %>% filter(is_constant_high==0), aes(fill = mean_total_nodes)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  facet_wrap(.~n, labeller = 'label_both') +
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(paste0(save_directory, "facet_n_same_scale_total_nodes.svg"), facet_n_same_scale_total_nodes, height = 4, width = 8)
