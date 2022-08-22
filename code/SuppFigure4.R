library(ggplot2)
library(cowplot)
library(data.table)
library(locuscomparer)
library(patchwork)

source("misc/locuscompare_updated.R")



## Supp Fig. 4a - Merged association plot
chr = "5"
population = "EUR"

### Scatter Plot
pu1_file = "../data/SuppFig4/PU1_27777.pu1bqtl.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')

gwas_file = "../data/SuppFig4/ZNF608.lym_count.metal.txt"
gwas_stat = read_metal(gwas_file, marker_col = 'rsid', pval_col = 'pval')

merged = merge(gwas_stat, pu1_stat, by = "rsid", suffixes = c("1", "2"), all = FALSE)
snp = 'rs12517864'

ld = retrieve_LD(chr, snp, population)
color = assign_color2(merged$rsid, snp, ld)

shape = ifelse(merged$rsid == snp, 23, 21)
names(shape) = merged$rsid

size = ifelse(merged$rsid == snp, 3, 2)
names(size) = merged$rsid

merged$label = ifelse(merged$rsid == snp, merged$rsid, '')



p_supp_4a <- ggplot(merged, aes(x=logp1, y=logp2)) +
  geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  xlab(expression(paste(-log[10],"(",italic(p)["Lym count"],")"))) +
  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=16), axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=14)) +
  ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf) +
  geom_point(aes(x=merged[merged$rsid ==snp,]$logp1, y=merged[merged$rsid ==snp,]$logp2), shape=23, size=3, fill="purple")



## Supp Fig. 4b - Z score plot

rs12517864_effects <- read.csv('../data/SuppFig4/rs12517864_effects.txt', header = T, sep='\t')

p_supp_4b <- ggplot(rs12517864_effects, aes(y = beta / se, x= source)) +
 geom_col(fill="gray", color="black") +
 theme_classic() + background_grid(major = 'y') +
 scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
 theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
       aspect.ratio = 1) +
 geom_vline(xintercept=0, size=1) + ylab("Z score (rs12517864)")


p_supp_4 <- p_supp_4a + p_supp_4b + plot_layout(widths = c(1,0.5))

ggplot2::ggsave('../figures/Supp4.pdf',
               plot = p_supp_4,
               device='pdf',
               width=300, height=150, units="mm")
