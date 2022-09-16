library(RColorBrewer)
library(tidyverse)
library(readr)
save_directory = "./../../../../plots/plots_for_figures/supplemental_figures/"

tiled_plots_fib <- read_delim("tiled plots_fib.txt", 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cast_k = read_csv('./../../../../extractedData/20210528_for_paper/analysis_of_larsson_2019/41586_2018_836_MOESM3_ESM_CAST.csv',
                  )
cast_k = cast_k %>% select(c(1,2))

cast_k= cast_k %>% separate(2, into =c('kon_cast', 'koff_cast'), sep = "(?<=\\d) +") %>% separate('kon_cast', into = c(NA, 'kon_cast'), sep = '[[]')
cast_k$kon_cast = as.numeric(cast_k$kon_cast)
cast_k$koff_cast = as.numeric(cast_k$koff_cast)

tiled_plots_fib = tiled_plots_fib %>% left_join(cast_k)
# pdf(file = "~/fib_point.pdf",   # The directory you want to save the file in
#     width = 20, # The width of the plot in inches
#     height = 20) # The height of the plot in inches
# 
# fib_point <- ggplot(tiled_plots_fib, aes(x = log10(kon), y = log10(koff))) +
#   geom_point(data = subset(subset(tiled_plots_fib, OR > 0), !is.na(OR)), 
#              aes(fill = log(OR)), 
#              alpha = 0.8,
#              width=0.03,
#              height=0.04) +
#   geom_tile(data = subset(subset(tiled_plots_fib, OR > 0), !is.na(OR)), 
#             aes(fill = log(OR)), 
#             alpha = 0.8,
#             width=0.03,
#             height=0.04) +
#   geom_tile(data = subset(tiled_plots_fib, OR == 0),
#             linetype = 0,
#             fill = "pink",
#             alpha = 0.5,
#             width=0.03,
#             height=0.04) +
#   geom_tile(data = subset(tiled_plots_fib, Ubiq==1),
#             linetype = 0,
#             fill = "purple", # subset of 0
#             width=0.05,
#             height=0.07) +
#   geom_tile(data = subset(tiled_plots_fib, is.na(OR)),
#             linetype = 0,
#             fill = "grey",
#             width=0.03,
#             height=0.04) +
#   geom_tile(data = subset(tiled_plots_fib, Zero_Zero==1),
#             linetype = 0,
#             fill = "gold",
#             width=0.05,
#             height=0.07) +
#   geom_tile(data = subset(tiled_plots_fib, FullCorr==1),
#             linetype = 0,
#             fill = "green", # subset of NaN
#             width=0.05,
#             height=0.07) +
#   geom_tile(data = subset(tiled_plots_fib, FullCorr_Ubiq==1),
#             linetype = 0,
#             fill = "red", # subset of 0/0
#             width=0.05,
#             height=0.07) +
#   theme_classic() + scale_fill_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9"))
# 
# dev.off()
# 


fibHM <- ggplot(tiled_plots_fib, 
                aes(x = log10(kon), y = log10(koff), z = log10(OR))) + 
  stat_summary_2d(data = tiled_plots_fib, 
                  fun = function(x) mean(x)) +
  theme_classic() + scale_fill_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9"))

ggsave(file.path(save_directory, 'fibHM_c57_no_nan.svg'), fibHM, width = 10, height = 6)


fibHM <- ggplot(tiled_plots_fib, 
                aes(x = log10(kon), y = log10(koff), z = log10(OR))) + 
  stat_summary_2d(data = tiled_plots_fib, 
                  fun = function(x) mean(x)) +
  stat_summary_2d(data = subset(tiled_plots_fib, OR <= 0),
                  fill = "pink",
                  alpha = 1,
                  aes(z=0)) + # low correlation
  # stat_summary_2d(data = subset(tiled_plots_fib, Ubiq==1),
  #                 fill = "purple",
  #                 aes(z=0)) + # subset of 0
  stat_summary_2d(data = subset(tiled_plots_fib, is.na(OR)),
                  alpha = 1,
                  fill = "grey",
                  aes(z=0)) + #NaN
# stat_summary_2d(data = subset(tiled_plots_fib, Zero_Zero==1),
#                 fill = "gold",
#                 aes(z=0)) + # 0/0
# stat_summary_2d(data = subset(tiled_plots_fib, FullCorr==1),
#                 fill = "green",
#                 aes(z=0)) + # subset of NaN
# stat_summary_2d(data = subset(tiled_plots_fib, FullCorr_Ubiq==1),
#                 fill = "red",
#                 aes(z=0)) + # subset of 0/0
theme_classic() + scale_fill_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9"))

ggsave(file.path(save_directory, 'fibHM_c57.eps'), fibHM, width = 10, height = 6)


fibHM_intermediate_check = fibHM + geom_rect(aes(xmin = -.3927, xmax = -.1683, ymin = 0.2596, ymax = 0.6057), fill = 'transparent', color = 'black')
ggsave(file.path(save_directory, 'fibHM_intermediate_check_c57.eps'), fibHM_intermediate_check, width = 10, height = 6)

fib_intermediate_zoom = ggplot(tiled_plots_fib %>% filter(log10(kon) > -.3927, log10(kon) < -.1683, log10(koff) > 0.2596, log10(koff) < 0.6057), 
                               aes(x = log10(kon), y = log10(koff), color = log10(OR))) +
  geom_point(shape=16) +
  geom_text_repel(data = tiled_plots_fib %>% filter(log10(kon) > -.3927, log10(kon) < -.1683, log10(koff) > 0.2596, log10(koff) < 0.6057, log10(OR) > 1),
                  aes(label = Gene), color = 'black')+ scale_color_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9")) +  theme_classic()

ggsave(file.path(save_directory, 'fib_intermediate_zoom_c57.eps'), fib_intermediate_zoom, width = 10, height = 6)


fibHM_cast <- ggplot(tiled_plots_fib, 
                aes(x = log10(kon_cast), y = log10(koff_cast), z = log10(OR))) + 
  stat_summary_2d(data = tiled_plots_fib, 
                  fun = function(x) mean(x)) +
  stat_summary_2d(data = subset(tiled_plots_fib, OR <= 0),
                  fill = "pink",
                  alpha = 1,
                  aes(z=0)) + # low correlation
  # stat_summary_2d(data = subset(tiled_plots_fib, Ubiq==1),
  #                 fill = "purple",
  #                 aes(z=0)) + # subset of 0
  stat_summary_2d(data = subset(tiled_plots_fib, is.na(OR)),
                  alpha = 1,
                  fill = "grey",
                  aes(z=0)) + #NaN
  # stat_summary_2d(data = subset(tiled_plots_fib, Zero_Zero==1),
  #                 fill = "gold",
  #                 aes(z=0)) + # 0/0
  # stat_summary_2d(data = subset(tiled_plots_fib, FullCorr==1),
  #                 fill = "green",
  #                 aes(z=0)) + # subset of NaN
  # stat_summary_2d(data = subset(tiled_plots_fib, FullCorr_Ubiq==1),
  #                 fill = "red",
  #                 aes(z=0)) + # subset of 0/0
  theme_classic() + scale_fill_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9"))
ggsave(file.path(save_directory, 'fibHM_cast.eps'), fibHM_cast, width = 10, height = 6)
