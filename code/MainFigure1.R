library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(stringr)
library(patchwork)

## Figure 1b  -- Enrichment of # blood cell trait association tagging variants


## Figure 1b - gkmSVM vs per-allele effect size

pu1_gkmsvm <- read.csv('../data/Fig1/PU1.gkmsvm.nominal.snps.txt', header = T, sep = '\t')


# gkmsvm
p_1b <- ggplot(pu1_gkmsvm, aes(x=ALT_gkmsvm - REF_gkmsvm, y=beta, color=-log10(p))) +
  geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
  geom_point(alpha=0.8) +
  labs(y="PU.1 bQTL effect size", x = expression(Delta~"gkm-SVM") , color=expression(paste(-log[10],"(",italic(p),")"))) +
  scale_colour_gradient(high = "black", low = "#E9E9E9") +
  theme_classic() + ylim(c(-1.9,1.9)) + xlim(c(-9,9)) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text=element_text(size=16),
        legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = "white", color = "black"),
        aspect.ratio=1) +
  geom_point(data=pu1_gkmsvm %>% filter(PU1_peak=='PU1_21900'), aes(x=ALT_gkmsvm - REF_gkmsvm, y=beta), shape=23, size = 3, fill='purple') +
  geom_point(data=pu1_gkmsvm %>% filter(PU1_peak=='PU1_53081'), aes(x=ALT_gkmsvm - REF_gkmsvm, y=beta), shape=23, size = 3, fill='purple') +
  geom_point(data=pu1_gkmsvm %>% filter(PU1_peak=='PU1_70939'), aes(x=ALT_gkmsvm - REF_gkmsvm, y=beta), shape=23, size = 3, fill='purple')



## Figure 1d  --  Number of motif-altering variants at the center of peak
motif_variant_location <- read.csv('../data/Fig1/motif_variant_location.txt', header=T, sep='\t')

p_1c <- motif_variant_location %>% mutate(pu1motif = factor(pu1motif, levels=c("none","motif"))) %>%
ggplot(aes(x=location, y=number)) +
  geom_col(aes(fill = location, alpha=pu1motif, linetype=pu1motif), width = 0.8, color = 'black') +
  theme_classic() +
  ylab(expression(paste("Number of PU.1 bQTLs"))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        aspect.ratio = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values =c('#F8766D', '#B8B4B4')) +
  scale_alpha_manual(values = c(0,1),labels = c('no','yes'), name = "PU.1 motif altered") +
  scale_linetype_manual(values = c('dashed', 'solid'), labels = c('no', 'yes'), name = expression(paste(italic("PU.1")," motif altered"))) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 6)) +
  guides(fill = "none", alpha = "none", linetype="none")


p_1_merged <- p_1b + p_1c + plot_layout(nrow = 2, heights = c(1, 0.7))


ggsave('../figures/Fig1.pdf',
       plot = p_1_merged,
       device='pdf',
       width=200, height=260, units="mm")
