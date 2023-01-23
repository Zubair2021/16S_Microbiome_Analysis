library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(ggstatsplot)
library(ggpubr)


# create a dataframe

df <- read.table("pielous_evenness.tsv", header = T)

# Create comparisons for significance

my_comparisons <- list(c("Control", "Pre_inoc"), c("Control", "Low_AL"),
                       c("Control", "High_AL"), c("Control", "Low_S1133"),
                       c("Control", "High_S1133"))


# create the boxplot

ggplot(df, aes(x =`Treatment`, y = `Pielou_Evenness`, fill = `Treatment`)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  theme_ipsum() +
  theme(axis.title.y = element_text(hjust = 0.5, size = 20, vjust = 6),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 15), 
        legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15)) +
  scale_fill_discrete(limits = c("Control", "Pre_inoc", "Low_AL", "High_AL", "Low_S1133", "High_S1133"), 
                      guide = guide_legend(keyheight =3)) +
  labs(y = "Pielou's Evenness") 


