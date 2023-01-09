library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(stringr)
library(patchwork)

## Figure S4b  -- Enrichment of # blood cell trait association tagging variants


## Figure S4b - gkmSVM vs per-allele effect size

pu1_gkmsvm <- read.csv('../data/ExtDataFig2/PU1_motif_score_bqtl.txt', header = T, sep = '\t')


# gkmsvm
p_supp_4b <- ggplot(pu1_gkmsvm, aes(x=ALT_gkmsvm - REF_gkmsvm, y=beta, color=-log10(p))) +
  geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
  geom_point(alpha=0.8, size = 3) +
  labs(y="PU.1 bQTL effect size", x = expression(Delta~"gkm-SVM") , color=expression(paste(-log[10],"(",italic(p),")"))) +
  scale_colour_gradient(high = "black", low = "#E9E9E9", limits = c(0,NA)) +
  theme_classic() + ylim(c(-1.7,1.7)) + xlim(c(-7,7)) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        legend.position = c(0.87, 0.25),
        aspect.ratio=1)

# p value for Pearson correlation
# cor.test(pu1_gkmsvm$ALT_gkmsvm-pu1_gkmsvm$REF_gkmsvm, pu1_gkmsvm$beta)$p.value

ggsave('../figures/FigS4b.pdf',
       plot = p_supp_4b,
       device='pdf',
       width=150, height=150, units="mm")
