library(ggplot2)
library(RColorBrewer)
library(ggridges)
library(viridis)
library(hrbrthemes)

library(readr)


###ODDS RATIOS###
compiled_OR_ridgelnePlot <- read_delim("~/compiled_OR_ridgelnePlot.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

compiled_OR_ridgelnePlot2 <- compiled_OR_ridgelnePlot
compiled_OR_ridgelnePlot2$Function <- factor(compiled_OR_ridgelnePlot$Function,     # Reorder factor levels
                                                c("mRNA Processing", "Ribosomal", "ETC", "TGFb", "SMAD", "NOTCH", "Stromal", "Alveolar", "Myeloblast", 
                                                  "Germ", "Fibroblast", "ESC"))

box_ESC_OR <- ggplot(compiled_OR_ridgelnePlot2, 
                      aes(x = OR_ESC, 
                          y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, OR_ESC>0), 
                      mapping = aes(x = log10(OR_ESC), 
                                    y = Function), 
                      bandwidth = 0.1, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_ESC_OR

box_MAF_OR <- ggplot(compiled_OR_ridgelnePlot2, 
                     aes(x = OR_MAF, 
                         y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, OR_MAF>0), 
                      mapping = aes(x = log10(OR_MAF), 
                                    y = Function), 
                      bandwidth = 0.1, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_MAF_OR

box_comb_OR <- ggplot(compiled_OR_ridgelnePlot2, 
                      aes(x = OR_MAF, 
                          y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, OR_MAF>0), 
                      mapping = aes(x = log10(OR_MAF), 
                                    y = Function, fill = "blue"), 
                      bandwidth = 0.1, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, point_color = "blue",
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, OR_ESC>0), 
                      mapping = aes(x = log10(OR_ESC), 
                                    y = Function, fill = "red"), 
                      bandwidth = 0.1, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, point_color = "red",
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()

box_comb_OR

###### kon ######
box_ESC_kon <- ggplot(compiled_OR_ridgelnePlot2, 
                     aes(x = kon_ESC, 
                         y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, kon_ESC>0), 
                      mapping = aes(x = log10(kon_ESC), 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_ESC_kon

box_MAF_kon <- ggplot(compiled_OR_ridgelnePlot2, 
                     aes(x = kon_MAF, 
                         y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, kon_MAF>0), 
                      mapping = aes(x = log10(kon_MAF), 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_MAF_kon


###### kon/koff ######
box_ESC_ratio <- ggplot(compiled_OR_ridgelnePlot2, 
                      aes(x = kon_MAF/koff_MAF, 
                          y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, kon_ESC>0), 
                      mapping = aes(x = log10(kon_ESC), 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_ESC_ratio

box_MAF_ratio <- ggplot(compiled_OR_ridgelnePlot2, 
                      aes(x = kon_MAF/koff_MAF, 
                          y = Function)) + 
  geom_density_ridges(subset(compiled_OR_ridgelnePlot2, kon_MAF>0), 
                      mapping = aes(x = log10(kon_MAF/koff_MAF), 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  
  theme_ridges()
box_MAF_ratio


########## Functional_Cat_Tiles
tile_celltype_esc <- ggplot(compiled_OR_ridgelnePlot2, 
                aes(x = log10(kon_ESC), y = log10(koff_ESC), z = log10(OR_ESC))) + 
  stat_summary_2d(data = compiled_OR_ridgelnePlot2, 
                  fun = function(x) mean(x)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "mRNA Processing"),
                  fill = "light green",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "ETC"),
                                             fill = "light green", 
                                             aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Ribosomal"),
                  fill = "light green",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "NOTCH"),
                  fill = "light blue",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "SMAD"),
                  fill = "light blue",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "TGFb"),
                  fill = "light blue",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Alveolar"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Myeloblast"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Germ"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Stromal"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Fibroblast"),
                  fill = "purple",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "ESC"),
                  fill = "black",
                  aes(z=0)) +
  theme_classic() + scale_fill_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9"))

tile_celltype_esc

tile_celltype_MAF <- ggplot(compiled_OR_ridgelnePlot2, 
                            aes(x = log10(kon_MAF), y = log10(koff_MAF), z = log10(OR_MAF))) + 
  stat_summary_2d(data = compiled_OR_ridgelnePlot2, 
                  fun = function(x) mean(x)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "mRNA Processing"),
                  fill = "light green",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "ETC"),
                  fill = "light green", 
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Ribosomal"),
                  fill = "light green",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "NOTCH"),
                  fill = "light blue",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "SMAD"),
                  fill = "light blue",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "TGFb"),
                  fill = "light blue",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Alveolar"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Myeloblast"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Germ"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Stromal"),
                  fill = "pink",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "Fibroblast"),
                  fill = "purple",
                  aes(z=0)) +
  stat_summary_2d(data = subset(compiled_OR_ridgelnePlot2, Function == "ESC"),
                  fill = "black",
                  aes(z=0)) +
  theme_classic() + scale_fill_gradientn(colours = c("#081d58", "#7fcdbb", "#ffffd9"))

tile_celltype_MAF



