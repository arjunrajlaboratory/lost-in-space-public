library(ggplot2)
library(RColorBrewer)
library(ggridges)
library(viridis)
library(hrbrthemes)

library(readr)

MAFresponseR <- read_delim("~/MAFresponseR.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)
MAFplot <- ggplot(MAFresponseR, aes(x=as.factor(Category), y=logOR)) + 
  geom_jitter(color="black", size=0.9, alpha=0.9)+ 
  geom_violin(fill="slateblue", alpha=0.4)+theme_bw() 
