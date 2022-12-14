library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(cowplot, quietly = T)
library(data.table, quietly = T)
library(AnnotationDbi, quietly = T)
library(ggbio, quietly = T)
library(GenomicRanges, quietly = T)
library(locuscomparer, quietly = T)
library(patchwork, quietly = T)


#### Figure S6a - GWAS effect sizes of rs5827412 for 5 blood cell traits

#rs5827412_gwas <- data.frame(source=as.factor(c("mono %", "mono #", "neut %", "neut #", "wbc #")),
#                             beta=c(-0.0451654, -0.029094, 0.0304417, 0.0262255, 0.0160937),
#                             se=c(0.00216547, 0.00213564, 0.00220315, 0.00219324, 0.0021686))

rs5827412_gwas <- read.csv('../data/FigS6/rs5827412_gwas.txt', header = T, sep ='\t')


p_supp_6_a <- ggplot(rs5827412_gwas, aes(x = beta, xmin=beta-1.96*se, xmax=beta+1.96*se, y= source)) +
  geom_point(fill="gray", color="black", size=2) +
  theme_minimal_vgrid() +
  geom_errorbar(color = "black", width = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limit=c(-0.06,0.06)) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=14),
        aspect.ratio = .7) +
  geom_vline(xintercept=0, size=1) + xlab("Effect size")

ggplot2::ggsave('../figures/FigS6a.pdf',
                plot = p_supp_6_a,
                device='pdf',
                width=150, height=100, units="mm")






#### Figure S6b - QTL plots for rs5827412

## PU.1 bQTL
PU1_71930 <- read.csv('../data/FigS6/qtl/PU1_71930.cpm.boxplot.txt', header = F, sep ='\t')

p_pu1_qtl_5 <- ggplot(PU1_71930, aes(x=as.factor(V2), y=V1)) +
  geom_boxplot(fill='#46ACC8', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(y="Read per million", title = "PU.1 ChIP") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

## chromatin accessibility QTL
ATAC_146492 <- read.csv('../data/FigS6/qtl/ATAC_146492.cpm.boxplot.txt', header = F, sep ='\t')

p_atac_qtl_5 <- ggplot(ATAC_146492, aes(x=as.factor(V2), y=V1)) +
  geom_boxplot(fill='#5166CC', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "ATAC") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

## H3K4me1 QTL
H3K4me1_19 <- read.csv('../data/FigS6/qtl/H3K4me1_19_18511451_18514338.cpm.boxplot.txt', header = F, sep ='\t')

p_h3k4me1_qtl_5 <- ggplot(H3K4me1_19, aes(x=as.factor(V2), y=V1)) +
  geom_boxplot(fill='#E7B800', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "H3K4me1", x="Allele dosage") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

## H3K27ac QTL
H3K27ac_19 <- read.csv('../data/FigS6/qtl/H3K27ac_19_18511593_18514423.cpm.boxplot.txt', header = F, sep ='\t')

p_h3k27ac_qtl_5 <- ggplot(H3K27ac_19, aes(x=as.factor(V2), y=V1)) +
  geom_boxplot(fill='#009E73', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "H3K27ac") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

## LRRC25 eQTL
LRRC25 <- read.csv('../data/FigS6/qtl/LRRC25.rpkm.boxplot.txt', header = F, sep ='\t')

p_lrrc25_qtl_5 <- ggplot(LRRC25, aes(x=as.factor(V2), y=V1)) +
  geom_boxplot(fill='darkgray', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="RPKM", title = expression(italic("LRRC25"))) +
  ylim(-0.15, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))


## plotting the panels together
p_supp_6_b_horizontal <- p_pu1_qtl_5 +
  p_atac_qtl_5 +
  p_h3k4me1_qtl_5 +
  p_h3k27ac_qtl_5 +
  p_lrrc25_qtl_5 +
  plot_layout(ncol = 5, widths = c(1,1,1,1,1))


ggplot2::ggsave('../figures/FigS6b.pdf',
                plot = p_supp_6_b_horizontal,
                device='pdf',
                width=250, height=60, units="mm")
