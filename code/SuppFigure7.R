library(ggplot2)
library(cowplot)
library(data.table)
library(locuscomparer)
library(patchwork)

source("misc/locuscompare_updated.R")



## Supp - MS GWAS
ms_2013_file = "../data/SuppFig7/ZC2HC1A_MS_2013.metal.txt"
gwas_stat = read_metal(ms_2013_file, marker_col = 'rsid', pval_col = 'pval')
#gwas_stat = get_position(gwas_stat)


pu1_file = "../data/SuppFig7/PU1_40678.pu1bqtl.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')
#pu1_stat = get_position(pu1_stat)

merged = merge(gwas_stat, pu1_stat, by = "rsid", suffixes = c("1", "2"), all = FALSE)
chr = '8'
snp = 'rs3808619'
population = 'EUR'

ld = retrieve_LD(chr, snp, population)
color = assign_color2(merged$rsid, snp, ld)

shape = ifelse(merged$rsid == snp, 23, 21)
names(shape) = merged$rsid

size = ifelse(merged$rsid == snp, 3, 2)
names(size) = merged$rsid

merged$label = ifelse(merged$rsid == snp, merged$rsid, '')



#merged_stat <- merge(gwas_stat, pu1_stat, by=c("rsid", "pos"))

#chr = '8'
#snp = 'rs3808619'
#population = 'EUR'
#ld = retrieve_LD(chr, snp, population)

#color = assign_color2(merged_stat$rsid, snp, ld)

#shape = ifelse(merged_stat$rsid == snp, 23, 21)
#names(shape) = merged_stat$rsid

#size = ifelse(merged_stat$rsid == snp, 3, 2)
#names(size) = merged_stat$rsid

#merged_stat$label = ifelse(merged_stat$rsid == snp, merged_stat$rsid, '')


#p_supp_7a <- make_scatterplot2(merged, title1 = 'GWAS', title2 = 'eQTL', color, shape,
#                     size, legend = T, legend_position = 'bottomright')


#p_supp_7a <- p_supp_7a +
#  xlab(expression(paste(-log[10],"(",italic(p)["MS GWAS"],")"))) +
#  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
#  theme_minimal() +
#  theme(aspect.ratio = 1,
#        axis.title.x = element_text(size=14),
#        axis.title.y = element_text(size=14),
#        axis.text.x = element_text(size=12),
#        axis.text.y = element_text(size=12))


p_supp_7a <- ggplot(merged, aes(x=logp1, y=logp2)) +
  geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  xlab(expression(paste(-log[10],"(",italic(p)["MS GWAS"],")"))) +
  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=16), axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=14)) +
  ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf) +
  geom_point(aes(x=merged[merged$rsid ==snp,]$logp1, y=merged[merged$rsid ==snp,]$logp2), shape=23, size=3, fill="purple")



## Supp Fig. 7b - Z score plot

rs3808619_effects <- read.csv('../data/SuppFig7/rs3808619_effects.txt', header = T, sep='\t')

p_supp_7b <- ggplot(rs3808619_effects, aes(y = beta / se, x= source)) +
 geom_col(fill="gray", color="black") +
 theme_classic() + background_grid(major = 'y') +
 scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
 theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
       aspect.ratio = 1) +
 geom_vline(xintercept=0, size=1) + ylab("Z score (rs3808619)")


p_supp_7 <- p_supp_7a + p_supp_7b + plot_layout(widths = c(1,0.5))

ggplot2::ggsave('../figures/Supp7.pdf',
               plot = p_supp_7,
               device='pdf',
               width=300, height=150, units="mm")
