library(ggplot2)
library(cowplot)
library(data.table)
library(locuscomparer)
library(dplyr)
library(patchwork)

source("misc/locuscompare_updated.R")


## Figure S9a - Merged association plot

## Lym count GWAS
lym_count_file = "../data/Fig6/ZC2HC1A.lym_count.metal.txt"
gwas_stat = read_metal(lym_count_file, marker_col = 'rsid', pval_col = 'pval')

pu1_file = "../data/Fig6/PU1_40678.pu1bqtl.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')

merged_lym = merge(gwas_stat, pu1_stat, by = "rsid", suffixes = c("1", "2"), all = FALSE)
chr = '8'
snp = 'rs3808619'
population = 'EUR'

ld = retrieve_LD(chr, snp, population)
color = assign_color2(merged_lym$rsid, snp, ld)

shape = ifelse(merged_lym$rsid == snp, 23, 21)
names(shape) = merged_lym$rsid

size = ifelse(merged_lym$rsid == snp, 3, 2)
names(size) = merged_lym$rsid

merged_lym$label = ifelse(merged_lym$rsid == snp, merged_lym$rsid, '')

p_supp_9a <- ggplot(merged_lym, aes(x=logp1, y=logp2)) +
  geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  xlab(expression(paste(-log[10],"(",italic(p)["Lym count GWAS"],")"))) +
  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=16), axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=14)) +
  ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf) +
  geom_point(aes(x=merged_lym[merged_lym$rsid ==snp,]$logp1, y=merged_lym[merged_lym$rsid ==snp,]$logp2), shape=23, size=3, fill="purple")



## Figure S9b - Z score plot

blood_rs3808619_effects <- read.csv('../data/FigS9/blood_rs3808619_effects.txt', header = T, sep='\t')

p_supp_9b <- blood_rs3808619_effects %>% mutate(source = factor(source, levels=c("PU.1 binding", "Lym count", "Lym %", "Mono %", "Neut %", "WBC count"))) %>%
  ggplot(aes(y = beta / se, x= source)) +
  geom_col(aes(fill=(beta < 0)), color="black") +
  theme_classic() + background_grid(major = 'y') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_discrete(c("#F8766D", "#00BFC4"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=14, angle = 45, hjust=1), axis.text.y = element_text(size=14),
       legend.position="none", aspect.ratio = 1) +
  geom_vline(xintercept=0, size=1) + ylab("Z score (rs3808619)")

## Figure S9c - Merged association plot with MS GWAS

ms_2013_file = "../data/FigS9/ZC2HC1A_MS_2013.metal.txt"
gwas_stat = read_metal(ms_2013_file, marker_col = 'rsid', pval_col = 'pval')

merged_ms = merge(gwas_stat, pu1_stat, by = "rsid", suffixes = c("1", "2"), all = FALSE)
chr = '8'
snp = 'rs3808619'
population = 'EUR'

ld = retrieve_LD(chr, snp, population)
color = assign_color2(merged_ms$rsid, snp, ld)

shape = ifelse(merged_ms$rsid == snp, 23, 21)
names(shape) = merged_ms$rsid

size = ifelse(merged_ms$rsid == snp, 3, 2)
names(size) = merged_ms$rsid

merged_ms$label = ifelse(merged_ms$rsid == snp, merged_ms$rsid, '')

p_supp_9c <- ggplot(merged_ms, aes(x=logp1, y=logp2)) +
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
  geom_point(aes(x=merged_ms[merged_ms$rsid ==snp,]$logp1, y=merged_ms[merged_ms$rsid ==snp,]$logp2), shape=23, size=3, fill="purple")



## Figure S9d - Z score plot

MS_rs3808619_effects <- read.csv('../data/FigS9/MS_rs3808619_effects.txt', header = T, sep='\t')

p_supp_9d <- ggplot(MS_rs3808619_effects, aes(y = beta / se, x= source)) +
 geom_col(fill="#F8766D", color="black") +
 theme_classic() + background_grid(major = 'y') +
 scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
 theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
       aspect.ratio = 1) +
 geom_vline(xintercept=0, size=1) + ylab("Z score (rs3808619)")


p_supp_9 <- p_supp_9a + p_supp_9b + p_supp_9c + p_supp_9d + plot_layout(nrow = 2, ncol = 2, widths = c(1,0.5))

ggplot2::ggsave('../figures/FigS9.pdf',
               plot = p_supp_9,
               device='pdf',
               width=300, height=300, units="mm")
