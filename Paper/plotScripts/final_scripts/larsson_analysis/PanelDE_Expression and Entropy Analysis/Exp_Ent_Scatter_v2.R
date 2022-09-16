library(ggplot2)
library(readxl)

compiled_exp_ent <- read_excel("~/all_Exp_EntVals.xlsx")

par(mfrow=c(2,2))

E57 <- ggplot(compiled_exp_ent, aes(x = log10(MeanExpression_E57), y = (Entropy_E57))) +
  geom_point() + theme_classic()
E57

ECAST <- ggplot(compiled_exp_ent, aes(x = log10(MeanExpression_ECAST), y = (Entropy_ECAST))) +
  geom_point() + theme_classic()
ECAST

F57 <- ggplot(compiled_exp_ent, aes(x = log10(MeanExpression_FC57), y = (Entropy_FC57))) +
  geom_point() + theme_classic()
F57

FCAST <- ggplot(compiled_exp_ent, aes(x = log10(MeanExpression_FCAST), y = (Entropy_FCAST))) +
  geom_point() + theme_classic()
FCAST


library(gridExtra)
grid.arrange(E57, ECAST, F57, FCAST, ncol = 2)
