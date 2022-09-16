library(tidyverse)
library(plotly)
library(cowplot)
library(ggnewscale)
library(viridis)
library(RColorBrewer)

gray_palette = c('#D8D5D5', '#898082', '#353131')

input_directory = './../../../extractedData/20210528_for_paper/vary_off_add_k_700params/'
save_directory = './../../../plots/plots_for_figures/'
normalize <- function(x) {
  return((x- min(x, na.rm=T)) /(max(x, na.rm=T)-min(x, na.rm=T)))
}
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
OR$mean_OR = as.numeric(OR$mean_OR)

colnames(bursts) = c("n_species", 'subnet', 'param_id', "rep", "number_bursts", "mean_burst_length", "mean_total_nodes", 'mean_concurrent_nodes', "mean_fraction_on")

metadata = read_csv(paste0(input_directory, 'metadata.csv'), col_names = FALSE)
colnames(metadata) = c("r_prod", "r_deg", "r_add", "r_off", "proddiff", "r_on", "k", "n", "p")
metadata = metadata %>% add_column(param_id = 1:nrow(metadata))

max_allele_xcorr = tibble()
max_allele_xcorr = read_csv(paste0('./../../../extractedData/20210528_for_paper/cross_correlation/', 'max_xcorrs_allele_species5-rep1.csv'), col_names = F)
colnames(max_allele_xcorr) = c("n_species", 'subnet', 'param_id', "rep", 'allele_node', "allele_cor", 'lag')
mean_allele_xcorr = max_allele_xcorr %>% group_by(n_species, subnet, param_id, rep) %>% summarize(mean_allele_cor = mean(allele_cor))

max_burst_xcorr = tibble()
max_burst_xcorr = read_csv(paste0('./../../../extractedData/20210528_for_paper/cross_correlation/', 'max_xcorrs_burst_species5-rep1.csv'), col_names = F)
colnames(max_burst_xcorr) = c("n_species", 'subnet', 'param_id', "rep", 'ispecies_1', 'ispecies_2', "gene_gene_cor", 'lag')
single_burst_xcorr = max_burst_xcorr %>% filter(ispecies_1 == 1 & ispecies_2 == 2)

sim_class = read_csv('./../../../extractedData/20210528_for_paper/vary_off_add_k_700params/long_simulation_class-rep1.csv', col_names = F)
colnames(sim_class) = c("n_species", 'subnet', 'param_id', "rep", "fraction_on", 'is_constant_high')

all_data = left_join(OR, metadata) %>% left_join(bursts)
all_data = all_data %>% mutate(log_mean_OR = log(mean_OR))
all_data = all_data %>%
  mutate(filtered_log_OR = ifelse(log_mean_OR < 0 | is.infinite(log_mean_OR), NA, log_mean_OR)) # Only keep ORs that are not below 1 and are finite

all_data = all_data %>%
  mutate(filtered_out_OR = ifelse(log_mean_OR < 0 & is.finite(log_mean_OR), 'sub-sampled', # OR < 1 implies anticorrelation, which is impossible and is therefore undersampled
                                  ifelse(log_mean_OR < 0 & is.infinite(log_mean_OR), 'zero', # OR = 0 is special case that implies a strong correlation in our data
                                         ifelse(is.na(log_mean_OR), NA, 'number'))))

all_data = all_data %>% left_join(sim_class) %>% left_join(mean_allele_xcorr) %>% left_join(single_burst_xcorr)

all_data = all_data %>% group_by(n_species, subnet) %>% filter(is_constant_high == 0) %>% mutate(normalized_log_OR = normalize(filtered_log_OR),
                                                                                                 normalized_mean_total_nodes = normalize(mean_total_nodes))

unique_r_off_to_plot = unique(all_data$r_off)[c(3, 5, 7)]

discrete_pal = viridis(6, option = 'plasma', dir = -1, begin = 0.1, end = 0.7)
purple_pal = c('#4E02C0', '#B686FE', '#E0CBFE')
orange_pal = c('#8d3a02', '#fc7922', '#feb686')

r_off_data = all_data %>% filter(k==110, n_species == 5, subnet ==1, r_off %in% unique_r_off_to_plot)

p1_OR = nls(normalized_log_OR ~ SSlogis(r_add, Asym,xmid,scal), data = r_off_data, subset = (r_off == unique_r_off_to_plot[1]))
p2_OR = update(p1_OR, subset = (r_off == unique_r_off_to_plot[2]))
p3_OR = update(p2_OR, subset = (r_off == unique_r_off_to_plot[3]))

pframe0 = data.frame(r_add = seq(0.001, 1, length.out = 10000))

predictions_OR= rbind(
  data.frame(pframe0, normalized_log_OR=  predict(p1_OR, pframe0),
             r_off=unique_r_off_to_plot[1], km = coef(p1_OR)[2]),
  data.frame(pframe0, normalized_log_OR=  predict(p2_OR, pframe0),
             r_off=unique_r_off_to_plot[2], km = coef(p2_OR)[2]),
  data.frame(pframe0, normalized_log_OR=  predict(p3_OR, pframe0),
             r_off=unique_r_off_to_plot[3], km = coef(p3_OR)[2])
)

p1_burst = nls(normalized_mean_total_nodes ~ SSlogis(r_add, Asym,xmid,scal), data = r_off_data, subset = (r_off == unique_r_off_to_plot[1]))
p2_burst = update(p1_burst, subset = (r_off == unique_r_off_to_plot[2]))
p3_burst = update(p2_burst, subset = (r_off == unique_r_off_to_plot[3]))

predictions_burst= rbind(
  data.frame(pframe0, normalized_mean_total_nodes=  predict(p1_burst, pframe0),
             r_off=unique_r_off_to_plot[1], km = coef(p1_burst)[2]),
  data.frame(pframe0, normalized_mean_total_nodes=  predict(p2_burst, pframe0),
             r_off=unique_r_off_to_plot[2], km = coef(p2_burst)[2]),
  data.frame(pframe0, normalized_mean_total_nodes=  predict(p3_burst, pframe0),
             r_off=unique_r_off_to_plot[3], km = coef(p3_burst)[2])
)


log_OR_unique_off=ggplot(r_off_data) +
  aes(x=r_add, y= normalized_log_OR, color = factor(r_off), label = param_id) +
  geom_point(shape=16)+
  geom_line(data = predictions_OR, aes(label = NULL)) +
  theme_classic() +
  geom_vline(data=predictions_OR, linetype = "dashed", aes(label = NULL, xintercept = km, color = factor(r_off))) +
  scale_x_log10() +
  scale_color_manual(values = purple_pal)

total_nodes_unique_off=ggplot(r_off_data) +
  aes(x=r_add, y= normalized_mean_total_nodes, color = factor(r_off)) +
  geom_point(shape=16)+
  geom_line(data = predictions_burst) +
  theme_classic() +
  geom_vline(data=predictions_burst, linetype = "dashed", aes(xintercept = km, color = factor(r_off))) +
  scale_x_log10()+
  scale_color_manual(values = purple_pal)
both_plots_unique_off = plot_grid(log_OR_unique_off,total_nodes_unique_off, nrow=2, align = 'hv')
ggsave(paste0(save_directory, 'r_off_unique_linegraphs_log_OR_total_nodes.svg'), both_plots_unique_off, width = 6, height = 4)



ggplot(all_data %>% filter(n_species ==5, is_constant_high==0, k==110, subnet==1)) +
  aes(x=normalized_log_OR, y=normalized_mean_total_nodes) +
  geom_point(shape=16)+
  theme_classic()
ggsave(paste0(save_directory, 'supplemental_figures/', 'OR_multi_burst_correspondence.svg'), width = 4.5, height = 3)


ggplot(all_data %>% filter(n_species ==5, is_constant_high==0, k==110, subnet==1)) +
  aes(x=normalized_log_OR, y=normalized_mean_total_nodes, color = as.factor(r_off)) +
  geom_point(shape=16)+
  facet_grid(.~r_off)+
  theme_classic() +
  scale_color_viridis(discrete = T, end = 0.9)
ggsave(paste0(save_directory, 'supplemental_figures/', 'OR_multi_burst_correspondence_faceted.svg'), width = 9, height = 2)

ggplot(all_data %>% filter(n_species ==5, is_constant_high==0, subnet==1)) +
  aes(x=normalized_log_OR, y=normalized_mean_total_nodes, color = as.factor(r_off)) +
  geom_point(shape=16)+
  theme_classic() +
  facet_grid(k~r_off)+
  scale_color_viridis(discrete = T, end = 0.9)
ggsave(paste0(save_directory, 'supplemental_figures/', 'OR_multi_burst_correspondence_over_k_faceted.svg'), width = 6.5, height = 4)


ggplot(all_data %>% filter(n_species ==5, is_constant_high==0, k==110, subnet %in% c(2:5))) +
  aes(x=normalized_log_OR, y=normalized_mean_total_nodes, color = as.factor(r_off)) +
  geom_point(shape=16)+
  theme_classic() +
  facet_grid(subnet~r_off)+
  scale_color_viridis(discrete = T, end = 0.9)
ggsave(paste0(save_directory, 'supplemental_figures/', 'OR_multi_burst_correspondence_over_subnet_faceted.svg'), width = 6.5, height = 2.6)

ggplot(all_data %>% filter(n_species ==5, k==110, subnet==1, is_constant_high==0)) +
  aes(x=mean_allele_cor, y=gene_gene_cor) +
  geom_point(shape=16)+
  labs(x='allele cor', y = 'gene cor')+
  theme_classic()
ggsave(paste0(save_directory, 'supplemental_figures/', 'allele_gene-gene_correlation_correspondence.svg'), width = 4.5, height = 3)

ggplot(all_data %>% filter(n_species ==5, subnet==1, is_constant_high==0)) +
  aes(x=mean_allele_cor, y=gene_gene_cor) +
  geom_point(shape=16)+
  labs(x='allele cor', y = 'gene cor')+
  theme_classic() +
  facet_grid(k~.)
ggsave(paste0(save_directory, 'supplemental_figures/', 'allele_gene-gene_correlation_correspondence_over_k.svg'), width = 1.75, height = 4)

ggplot(all_data %>% filter(n_species ==5, k==110, is_constant_high==0, subnet %in% c(2:5))) +
  aes(x=mean_allele_cor, y=gene_gene_cor) +
  geom_point(shape=16)+
  labs(x='allele cor', y = 'gene cor')+
  theme_classic() +
  facet_grid(subnet~.)
ggsave(paste0(save_directory, 'supplemental_figures/', 'allele_gene-gene_correlation_correspondence_over_subnet.svg'), width = 1.75, height = 2.6)