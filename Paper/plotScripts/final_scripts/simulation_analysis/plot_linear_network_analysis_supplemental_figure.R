library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
gray_palette = c('#D8D5D5', '#898082', '#353131')

input_directory = './../../../extractedData/20210528_for_paper/linear_networks/'
save_directory = './../../../plots/plots_for_figures/supplemental_figures/'
#Read in data and add columns (no columns in MATLAB)
sim_classification = read_csv(paste0(input_directory, 'long_simulation_class-rep', 1, '.csv'), col_names = F)
colnames(sim_classification) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

bursts = tibble()
files_to_read = paste0(input_directory, 'long_bursts-species5.csv')
bursts = bind_rows(bursts, files_to_read %>% 
                     lapply(read_csv, col_names = FALSE))
bursts = bursts[,1:8]

per_gene_above = read_csv(paste0(input_directory, 'per_gene_above-species5.csv'), col_names=F)
OR = tibble()
files_to_read = paste0(input_directory, 'long_OR-species5.csv')
OR = bind_rows(OR, files_to_read %>% 
                 lapply(read_csv, col_names = FALSE))


cross_corrs = read_csv(paste0(input_directory, 'max_xcorrs_burst_species5-rep1.csv'), col_names = F)
colnames(cross_corrs) = c("n_species", 'subnet', 'param_id', "rep", 'ispecies_1', 'ispecies_2', "cor", 'lag')
allele_cross_corrs = read_csv(paste0(input_directory, 'max_xcorrs_allele_species5-rep1.csv'), col_names = F)
colnames(allele_cross_corrs) = c("n_species", 'subnet', 'param_id', "rep", 'allele_node', "allele_cor", 'lag')
allele_cross_corrs$allele_node = paste0('allele_node_', allele_cross_corrs$allele_node)

colnames(OR) = c("n_species", 'subnet', 'param_id', "rep", 'allele_node_1', 'allele_node_2', 'allele_node_3', 'allele_node_4', 'allele_node_5')
colnames(per_gene_above) = c("n_species", 'subnet', 'param_id', "rep", 'node_1', 'node_2', 'node_3', 'node_4', 'node_5')

colnames(bursts) = c("n_species", 'subnet', 'param_id', "rep", "number_bursts", "mean_burst_length", "mean_total_nodes", 'mean_concurrent_nodes')

metadata = read_csv(paste0(input_directory, 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")

metadata = metadata %>% add_column(param_id = 1:nrow(metadata))

all_data = left_join(OR, metadata) %>% left_join(per_gene_above) %>% left_join(bursts) %>% left_join(sim_classification)
cross_corrs = left_join(cross_corrs, metadata)

cross_corrs = cross_corrs %>% unite('combined', ispecies_1:ispecies_2, remove = F)
all_data_long = all_data %>% pivot_longer(cols=allele_node_1:allele_node_5, names_to= 'allele_node', values_to = 'allelic_OR')
all_data_long = all_data_long %>% pivot_longer(cols=node_1:node_5, names_to = 'node_above_thresh', values_to = 'time_above_thresh')

all_data_long$allele_node = factor(all_data_long$allele_node, levels = c('allele_node_5', 'allele_node_4', 'allele_node_3', 'allele_node_2', 'allele_node_1'))

direct_cors = cross_corrs %>% filter(combined %in% c('1_2', '2_3', '3_4', '4_5'))
direct_cors = direct_cors %>% mutate(allele_node = case_when(combined == '1_2' ~ 'allele_node_1',
                                                         combined == '2_3' ~ 'allele_node_2',
                                                         combined == '3_4' ~ 'allele_node_3',
                                                         combined == '4_5' ~ 'allele_node_4'))
combined_OR_direct_cors = all_data_long %>% select(n_species, subnet, param_id, rep, allele_node, allelic_OR) %>% right_join(direct_cors) %>% distinct(.)
combined_OR_direct_cors = allele_cross_corrs %>% select(-lag) %>% right_join(combined_OR_direct_cors)
combined_OR_direct_cors$allele_node = factor(combined_OR_direct_cors$allele_node, levels = c('allele_node_5', 'allele_node_4', 'allele_node_3', 'allele_node_2', 'allele_node_1'))

direct_cors$combined = factor(direct_cors$combined, levels = c('4_5', '3_4', '2_3', '1_2'))

first_minus_last_node = ggplot(all_data %>% filter(k == 110, n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off) +
  geom_tile(aes(fill = node_1-node_5)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(file.path(save_directory, 'linear_first_minus_last_node.eps'), first_minus_last_node, width = 6, height = 4)

facet_OR_node_heatmaps=ggplot(all_data_long %>% filter(k == 110, n_species ==5)) +
  aes(x = r_add, y = r_off) +
  geom_tile(data = all_data_long %>% filter(k == 110, n_species ==5, is_constant_high == 0), aes(fill = log(allelic_OR))) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5') +
  scale_fill_viridis(direction = 1, option = "mako") +
  new_scale_fill()+
  geom_tile(aes(fill=as.factor(is_constant_high)))+
  scale_fill_manual(values = c('transparent', gray_palette[3]))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold')) +
  facet_grid(.~allele_node)
ggsave(file.path(save_directory, 'linear_facet_OR_node_heatmaps.eps'), facet_OR_node_heatmaps, width = 10, height = 6)

mean_total_nodes = ggplot(all_data %>% filter(k == 110, n_species ==5, subnet == 1)) +
  aes(x = r_add, y = r_off) +
  geom_tile(aes(fill = mean_total_nodes)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = 'K = 110, nodes = 5, subnet = 1') +
  scale_fill_viridis(direction = 1, option = "mako")
ggsave(file.path(save_directory, 'mean_total_nodes.eps'), mean_total_nodes, width = 6, height = 4)

direct_gene_correlations_facet = ggplot(direct_cors) +
  aes(x = r_add, y = r_off, fill = cor) +
  geom_tile()+
  scale_x_log10() +
  scale_y_log10() +
  theme_cowplot()+
  scale_fill_viridis(direction = 1, option = "mako") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))+
  facet_grid(.~combined)
ggsave(file.path(save_directory, 'linear_direct_gene_correlations_facet.eps'), direct_gene_correlations_facet, width = 10, height = 6)

log_allelic_OR_v_gene_cor = ggplot(combined_OR_direct_cors) +
  aes(x=log(allelic_OR), y = cor) +
  geom_point(shape = 16) +
  facet_grid(.~allele_node) +
  theme_cowplot()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(file.path(save_directory, 'linear_log_allelic_OR_v_gene_cor.eps'),log_allelic_OR_v_gene_cor, width = 10, height = 6)

allele_cor_v_gene_cor=ggplot(combined_OR_direct_cors) +
  aes(x=allele_cor, y = cor) +
  geom_point(shape=16) +
  facet_grid(.~allele_node)+
  theme_cowplot()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.title = element_text(size = 16, face = 'bold'), 
        axis.title = element_text(face = 'bold'))
ggsave(file.path(save_directory, 'linear_allele_cor_v_gene_cor.eps'),allele_cor_v_gene_cor, width = 10, height = 6)
