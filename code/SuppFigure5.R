library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(patchwork)


#### SUPP

#rs5827412_gwas <- data.frame(source=as.factor(c("mono %", "mono #", "neut %", "neut #", "wbc #")),
#                             beta=c(-0.0451654, -0.029094, 0.0304417, 0.0262255, 0.0160937),
#                             se=c(0.00216547, 0.00213564, 0.00220315, 0.00219324, 0.0021686))

rs5827412_gwas <- read.csv('../data/SuppFig5/rs5827412_gwas.txt', header = T, sep ='\t')


p_supp_5_a <- ggplot(rs5827412_gwas, aes(x = beta, xmin=beta-1.96*se, xmax=beta+1.96*se, y= source)) +
  geom_point(fill="gray", color="black", size=2) +
  theme_minimal_vgrid() +
  geom_errorbar(color = "black", width = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limit=c(-0.06,0.06)) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=14),
        aspect.ratio = .7) +
  geom_vline(xintercept=0, size=1) + xlab("Effect size")

ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Supp5_a.pdf',
                plot = p_supp_5_a,
                device='pdf',
                width=150, height=100, units="mm")


## QTLs

## PU.1 peak

PU1_71930 <- read.csv('~/Projects/PU1_gwas/plots/PU1_71930/PU1_71930.cpm.boxplot.txt', header = F, sep ='\t')

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

## ATAC peak over PU.1 peak

ATAC_146492 <- read.csv('~/Projects/PU1_gwas/plots/PU1_71930/ATAC_146492.cpm.boxplot.txt', header = F, sep ='\t')

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

## H3K4me1

H3K4me1_19 <- read.csv('~/Projects/PU1_gwas/plots/PU1_71930/H3K4me1_19_18511451_18514338.cpm.boxplot.txt', header = F, sep ='\t')

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

## H3K27ac

H3K27ac_19 <- read.csv('~/Projects/PU1_gwas/plots/PU1_71930/H3K27ac_19_18511593_18514423.cpm.boxplot.txt', header = F, sep ='\t')

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

## LRRC25

LRRC25 <- read.csv('~/Projects/PU1_gwas/plots/PU1_71930/LRRC25.rpkm.boxplot.txt', header = F, sep ='\t')

p_lrrc25_qtl_5 <- ggplot(LRRC25, aes(x=as.factor(V2), y=V1)) +
  geom_boxplot(fill='darkgray', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="RPKM", title = expression(italic("LRRC25"))) +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))



p_supp_5_c_horizontal <- p_pu1_qtl_5 +
  p_atac_qtl_5 +
  p_h3k4me1_qtl_5 +
  p_h3k27ac_qtl_5 +
  p_lrrc25_qtl_5 +
  plot_layout(ncol = 5, widths = c(1,1,1,1,1))

p_supp_5_c_horizontal

ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Supp5_c.pdf',
                plot = p_supp_5_c_horizontal,
                device='pdf',
                width=250, height=60, units="mm")


## Not used  - LRRC25 eQTL plot
LRRC25_eqtl <- read.csv('~/Google Drive/Manuscript/PU1_bloodcelltrait/data/LRRC.eqtl.txt', header = T, sep ='\t')



ggplot(LRRC25_eqtl, aes(x= interaction(Study, CellType), y = mlog10p, colour = CellType, fill = CellType,shape=Direction)) + theme_classic() +
  geom_point(size = 3.5) + scale_shape_manual(name = "eQTL effect direction", values=c(25,24))

ggplot(LRRC25_eqtl, aes(x= interaction(Study, CellType), y = Beta, ymin=Beta - Beta_se*1.96, ymax=Beta + Beta_se*1.96, colour = CellType)) +
  theme_minimal_hgrid() +
  geom_point(size = 3.5) + ylim(c(-1,0.02)) + xlab("") +
  scale_x_discrete(position='top') + geom_errorbar(width = 0.2) + geom_hline(yintercept = 0, size=1.5) +
  theme(axis.text.x=element_text(size=10, angle = 90))



ggplot(LRRC25_eqtl[LRRC25_eqtl$CellType =='monocyte',], aes(x= interaction(Study, CellType_2), y = Beta, ymin=Beta - Beta_se*1.96, ymax=Beta + Beta_se*1.96, colour = CellType)) +
  theme_minimal_hgrid() +
  geom_point(size = 3.5) + geom_errorbar(width = 0.1) +
  ylim(c(-0.15,0.02)) + xlab("") + ylab("Effect size") +
  scale_x_discrete(position='top') +  geom_hline(yintercept = 0, size=1.5) +
  theme(axis.text.x=element_text(size=10, angle = 90))

ggplot(LRRC25_eqtl[LRRC25_eqtl$CellType_2 =='classical_monocyte',], aes(x= Study, y = Beta, ymin=Beta - Beta_se*1.96, ymax=Beta + Beta_se*1.96, colour = CellType)) +
  theme_minimal_hgrid() +
  geom_point(size = 3.5) + geom_errorbar(width = 0.1) +
  ylim(c(-0.15,0.02)) + xlab("") + ylab("Effect size") +
  scale_x_discrete(position='top') +  geom_hline(yintercept = 0, size=1.5) +
  theme(axis.text.x=element_text(size=10, angle = 90))

##






## Supp Figure - LRRC25 expression across all blood cell types
heme_RNA <- read.csv('../data/misc/heme_RNA.txt', header = T, sep ='\t', row.names = 1)
LRRC25_RNA <- transpose((heme_RNA / colSums(heme_RNA) * 10^6)["LRRC25",1:49])
celltype <- c("HSC", "HSC", "HSC", "HSC", "MPP", "MPP", "MPP", "MPP", "LMPP", "LMPP", "LMPP", "CMP", "CMP", "CMP", "CMP", "GMP", "GMP", "GMP", "GMP", "MEP", "MEP", "MEP", "MEP", "Mono", "Mono", "Mono", "Mono", "CD4T", "CD4T", "CD4T", "CD4T", "CD8T", "CD8T", "CD8T", "CD8T", "NK", "NK", "NK", "NK", "B", "B", "B", "B", "CLP", "CLP", "CLP", "Ery", "Ery", "Ery")
LRRC25_RNA["celltype"] <- celltype

p_lrrc25_rna_all <- ggplot(LRRC25_RNA, aes(x=reorder(as.factor(celltype), -V1), y=(V1+0.1), fill=(celltype=='Mono'))) +
  geom_boxplot(color="black") + geom_jitter(shape=16, position=position_jitter(0.1)) +
  theme_classic() + scale_y_continuous(trans='log10') +
  scale_fill_manual(values = c("#B8B4B4", "#FF7373")) +
  labs(y=expression(paste(italic('LRRC25'), " count per million")), x = "Cell type") +
  theme(axis.title.x = element_text(size=18, face="bold"), axis.title.y = element_text(size=18, face="bold"),
        axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        legend.position="none")
p_lrrc25_rna_all

ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Supp5_b.pdf',
                plot = p_lrrc25_rna_all,
                device='pdf',
                width=250, height=150, units="mm")
