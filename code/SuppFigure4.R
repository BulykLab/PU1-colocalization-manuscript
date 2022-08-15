library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(patchwork)
library(dplyr)
library(ggplot2)

source("misc/locuscompare_updated.R")



#### Supp Figures ####


### Scatter Plot
pu1_file = "~/Google Drive/Manuscript/PU1_bloodcelltrait/data/ZNF608/PU1_27777.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')

gwas_file = "~/Google Drive/Manuscript/PU1_bloodcelltrait/data/ZNF608/ZNF608.lym_count.metal.txt"
gwas_stat = read_metal(gwas_file, marker_col = 'rsid', pval_col = 'pval')

merged = merge(gwas_stat, pu1_stat, by = "rsid", suffixes = c("1", "2"), all = FALSE)
snp = 'rs12517864'

ld = retrieve_LD(chr, snp, population)
color = assign_color(merged$rsid, snp, ld)

shape = ifelse(merged$rsid == snp, 23, 21)
names(shape) = merged$rsid

size = ifelse(merged$rsid == snp, 3, 2)
names(size) = merged$rsid

merged$label = ifelse(merged$rsid == snp, merged$rsid, '')

p_supp4_a = make_scatterplot2(merged, title1 = 'GWAS', title2 = 'eQTL', color, shape,
                     size, legend = F, legend_position = 'bottomright')


p_supp4_a <- p_supp4_a +
  xlab(expression(paste(-log[10],"(",italic(p)["Lym #"],")"))) +
  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Supp4_a.pdf',
                plot = p_supp4_a,
                device='pdf',
                width=130, height=130, units="mm")
