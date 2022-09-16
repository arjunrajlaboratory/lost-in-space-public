library(ggplot2)
library(RColorBrewer)
library(ggridges)
library(viridis)
library(hrbrthemes)

library(readr)

entropy_expression_kon_koff <- read_delim("~/entropy_expression_kon_koff.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

entropy_expression_kon_koff2 <- entropy_expression_kon_koff
entropy_expression_kon_koff2$Function <- factor(entropy_expression_kon_koff$Function,     # Reorder factor levels
       c("mRNA Processing", "Ribosomal", "ETC", "TGFb", "SMAD", "NOTCH", "Stromal", "Alveolar", "Myeloblast", 
         "Germ", "Fibroblast", "ESC"))
entropy_expression_kon_koff2$Function2 <- factor(entropy_expression_kon_koff$Function2,     # Reorder factor levels
                                                c("mRNA Processing", "Ribosomal", "ETC", "Signaling", "Stromal", "Alveolar", "Myeloblast", 
                                                  "Germ", "Fibroblast", "ESC"))

box_ESC_C57 <- ggplot(entropy_expression_kon_koff2, 
                      aes(x = Entropy_ESC_C57, 
                          y = Function)) + 
                  geom_density_ridges(entropy_expression_kon_koff2, 
                                      mapping = aes(x = Entropy_ESC_C57, 
                                          y = Function), 
                                      #point_color = entropy_expression_kon_koff$`Mean Expression_ESC_C57`, 
                                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                                       jittered_points = TRUE,
                                       position = position_points_jitter(width = 0.05, height = 0),
                                       point_shape='|', point_size = 1, 
                                       point_alpha = 1, alpha = 0.7, 
                                       quantile_lines = TRUE, quantiles = 2) +
                  #scale_fill_manual(aesthetics = "point_color", 
                   #                     values = entropy_expression_kon_koff$`Mean Expression_ESC_C57`) +
                  theme_ridges()
box_ESC_C57

box_ESC_CAST <- ggplot(entropy_expression_kon_koff2, 
                      aes(x = Entropy_ESC_CAST, 
                          y = Function)) + 
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_ESC_CAST, 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_ESC_CAST

box_MAF_CAST <- ggplot(entropy_expression_kon_koff2, 
                       aes(x = Entropy_MAF_CAST, 
                           y = Function)) + 
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_MAF_CAST, 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_MAF_CAST

box_MAF_C57 <- ggplot(entropy_expression_kon_koff2, 
                       aes(x = Entropy_MAF_C57, 
                           y = Function)) + 
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_MAF_C57, 
                                    y = Function), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, 
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_MAF_C57

box_comb_CAST_Ent <- ggplot(entropy_expression_kon_koff2, 
                         aes(x = Entropy_ESC_CAST, 
                             y = Function)) + 
  
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_ESC_CAST, 
                                    y = Function, fill = "blue"), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, point_color = "blue",
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_MAF_CAST, 
                                    y = Function, fill = "red"), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, point_color = "red",
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_comb_CAST_Ent  

box_comb_C57_Ent <- ggplot(entropy_expression_kon_koff2, 
                            aes(x = Entropy_ESC_C57, 
                                y = Function)) + 
  
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_ESC_C57, 
                                    y = Function, fill = "blue"), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, point_color = "blue",
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  geom_density_ridges(entropy_expression_kon_koff2, 
                      mapping = aes(x = Entropy_MAF_C57, 
                                    y = Function, fill = "red"), 
                      bandwidth = 0.08, rel_min_height = 0.001, scale = 0.8, 
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape='|', point_size = 1, point_color = "red",
                      point_alpha = 1, alpha = 0.7, 
                      quantile_lines = TRUE, quantiles = 2) +
  theme_ridges()
box_comb_C57_Ent
