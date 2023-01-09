library(ggplot2, quietly=T)
library(dplyr, quietly = T)
library(cowplot, quietly = T)

## Figure 3a - Lollipop plot for the proportion colocalized

# Reading number of colocalized loci by traits
coloc_counts <- read.csv('../data/Fig3/JLIM_COLOC.fdr.counts.txt', header = T, sep ='\t')

p_3a <- coloc_counts %>% mutate(Category = factor(Category, levels=c("Granulocyte", "Monocyte", "Lymphocyte", "MatureRed", "ImmatureRed","Platelet"))) %>%
  mutate(Trait = factor(Trait, levels=c("plt_dist_width", "plt_crit", "plt_count","MPV", "ret_percentage",  "ret_count", "MRV",
                                        "imm_ret_frac","high_light_scat_ret_percentage", "high_light_scat_ret_count", "rbc_dist_width",  "rbc_count",
                                        "MSCV", "MCV", "MCHC", "MCH", "ht_percentage", "hb_concentration",
                                        "lym_percentage",  "lym_count", "mono_percentage", "mono_count", "wbc_count", "neut_percentage", "neut_count",
                                        "eosino_percentage","eosino_count", "baso_percentage",  "baso_count"))) %>%
  filter(Proportion_colocalized > 0) %>%
  ggplot(aes(x = Trait, y=Proportion_colocalized)) + theme_classic() + background_grid(major = 'x', minor='x') +
  geom_point(size = 3.5, aes(colour = Category)) +
  geom_segment(aes(xend=Trait, y=0, yend=Proportion_colocalized, colour = Category), size=1.1, alpha=0.5) +
  coord_flip() + scale_y_reverse(position = "left", expand = expansion(mult = c(0.1, 0))) + scale_x_discrete(position = "top") +
  xlab("") + ylab("Proportion colocalized") +
  scale_color_manual(values=c('#F98400', '#351C75','#40A2EB', '#A61C00', '#FFB6C1', '#008000'), name="Trait group") +
  theme(legend.position = "left", axis.text.x=element_text(size=16), axis.text.y=element_blank(),
        axis.title.x = element_text(size=18),
        legend.title = element_text(size=17), legend.text = element_text(size=15),
        panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent'))


ggsave('../figures/Fig3a.pdf', plot = p_3a, bg='transparent', device='pdf', width=6.97, height=5.42, units="in")
