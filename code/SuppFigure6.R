library(ggplot2)
library(cowplot)
library(data.table)
library(locuscomparer)
library(patchwork)

source("misc/locuscompare_updated.R")



## Supp - MS GWAS
ms_2013_file = "../data/SuppFig6/ZC2HC1A_MS_2013.metal.txt"
gwas_stat = read_metal(ms_2013_file, marker_col = 'rsid', pval_col = 'pval')
gwas_stat = get_position(gwas_stat)


pu1_file = "../data/SuppFig6/PU1_40678.pu1bqtl.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')
pu1_stat = get_position(pu1_stat)


merged_stat <- merge(gwas_stat, pu1_stat, by=c("rsid", "pos"))

chr = '8'
snp = 'rs3808619'
population = 'EUR'
ld = retrieve_LD(chr, snp, population)

color = assign_color2(merged_stat$rsid, snp, ld)

shape = ifelse(merged_stat$rsid == snp, 23, 21)
names(shape) = merged_stat$rsid

size = ifelse(merged_stat$rsid == snp, 3, 2)
names(size) = merged_stat$rsid



p_supp_6a <- ggplot(merged_stat) + geom_point(aes(x=-log10(pval.x), y=-log10(pval.y), fill = rsid, size = rsid, shape = rsid)) +
  theme_minimal() +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  xlab(expression(paste(-log[10],"(",italic(p)["MS GWAS"],")"))) +
  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=16), axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=16), axis.text.y = element_text(size=14)) +
  geom_point(aes(x=gwas_stat[gwas_stat$pos ==79577726,]$logp, y=pu1_stat[pu1_stat$pos ==79577726,]$logp), shape=23, size=3, fill="purple")


## Effect size plot

### NEED TO WORK ON THIS
#rs3808619_effects <- data.frame(source=as.factor(c("MS", "PU.1 binding")),
#                             beta=c(0.11422114409002286, -0.029094),
#                             se=c(0.018738168792812798, 0.00213564))

rs3808619_effects <- read.csv('../data/SuppFig6/rs3808619_effects.txt', header = T, sep='\t')

p_supp_6b <- ggplot(rs3808619_effects, aes(y = beta / se, x= source)) +
 geom_col(fill="gray", color="black") +
 theme_classic() + background_grid(major = 'y') +
 scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
 theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
       aspect.ratio = 1) +
 geom_vline(xintercept=0, size=1) + ylab("Z score (rs3808619)")


p_supp_6 <- p_supp_6a + p_supp_6b + plot_layout(widths = c(1,0.6))

ggplot2::ggsave('../figures/Supp6.pdf',
               plot = p_supp_6,
               device='pdf',
               width=300, height=150, units="mm")
