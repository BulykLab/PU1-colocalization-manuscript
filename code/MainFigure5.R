library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(patchwork)

source("misc/locuscompare_updated.R")


## Figure 5a - Colocalization plot

gwas_file="../data/Fig5/coloc/mono_p.LRRC25.withbeta.filtered.txt"
pu1_file="../data/Fig5/coloc/PU1_71930.bqtl.withbeta.filtered.txt"
ld_file="../data/Fig5/coloc/rs5827412.ld.txt"
snp = 'rs5827412'

gwas_stat <- read.table(gwas_file, header=F)
colnames(gwas_stat) <- c("rsid", "chr", "pos", "pval", "beta", "se")
gwas_stat$z <- gwas_stat$beta / gwas_stat$se

pu1_stat <- read.table(pu1_file, header=F)
colnames(pu1_stat) <- c("rsid", "chr", "pos", "pval", "beta")
pu1_stat$z <- qnorm((pu1_stat$pval)/2, lower.tail = F) * (pu1_stat$beta / abs(pu1_stat$beta))

merged_stat <- merge(gwas_stat, pu1_stat, by=c("rsid", "pos"))


ld = read.table(ld_file, header=T)

color = assign_color2(merged_stat$rsid, snp, ld)

shape = ifelse(merged_stat$rsid == snp, 23, 21)
names(shape) = merged_stat$rsid

size = ifelse(merged_stat$rsid == snp, 3, 2)
names(size) = merged_stat$rsid


p_5_a <- ggplot(merged_stat) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(aes(x=z.x, y=z.y, fill = rsid, size = rsid, shape = rsid)) +
  theme_minimal() +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  xlab(expression(paste(italic(Z)["Mono %"]))) +
  ylab(expression(paste(italic(Z)["PU.1 bQTL"]))) +
  xlim(c(-21,12)) + ylim(c(-7,4)) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))



ggplot2::ggsave('../figures/Fig5_a.pdf',
                plot = p_5_a,
                device='pdf',
                width=110, height=110, units="mm")





## Figure 5c - rs5827412 MPRA
rs5827412_MPRA <- data.frame(source=as.factor(c("Tewhey2016", "Abell2022")), logSkew=c(-0.390574802, -0.342679482),logSkewSE=c(0.09816138,0.307302851))

p_rs5827412_mpra <- ggplot(rs5827412_MPRA, aes(x= source, y = logSkew, ymin=logSkew-1.96*logSkewSE, ymax=logSkew+1.96*logSkewSE)) +
  geom_col(fill="gray", color="black") +
  theme_classic() + scale_x_discrete(position='top') +
  geom_errorbar(color = "black", width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(-1,0.3)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=14),
        aspect.ratio = 1.1) +
  geom_hline(yintercept=0, size=1) + ylab("rs5827412 \n allelic skew")

p_rs5827412_mpra


### Figure 5d - PU1KO
pu1ko_atac_lrrc25 <- read.csv('../data/Fig5/pu1_ko_atac_lrrc25.txt', header = T, sep ='\t')

p_pu1ko <- pu1ko_atac_lrrc25 %>% mutate(condition = factor(condition, levels=c("SPI1+/+", "SPI-/-"))) %>%
  ggplot(aes(x=as.factor(condition), y=tmm)) +
  geom_boxplot(fill='white', color="black") + theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1), size=3) +
  labs(y="Accessibility (CPM)") + ylim(0, NA) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    aspect.ratio = 1.1)


## Figure 5e - LRRC25 expression across blood cell types
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


p_lrrc25_rna <- LRRC25_RNA %>% filter(celltype %in% c("HSC", "MPP", "CMP", "GMP", "Mono")) %>%
  mutate(celltype = factor(celltype, levels=c("HSC", "MPP", "CMP", "GMP", "Mono"))) %>%
  ggplot(aes(x=celltype, y=(V1+0.1), fill=celltype)) +
    geom_boxplot(color="black") + geom_jitter(shape=16, position=position_jitter(0.1)) +
    theme_classic() + scale_y_continuous(trans='log10') +
    labs(y=expression(paste(italic('LRRC25'), " CPM"))) +
    scale_fill_manual(values = c("#333333", "#7F3939", "#B25050", "#E56767", "#FF7373")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14),
          axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
          aspect.ratio = 0.8,
          legend.position = "none")
p_lrrc25_rna

p_5_cde <- p_rs5827412_mpra + p_pu1ko + p_lrrc25_rna + plot_layout(ncol=3, widths = c(1,1,1.4))


ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Fig5_cde.pdf',
                plot = p_5_cde,
                device='pdf',
                width=250, height=75, units="mm")


#### Figure 5f,g - LRRC25 eQTL (monocyte) and gene plot
### Gene plot

chr <- "19"
#start <- 18450000
#end <- 18550000

start <- 18495845
end <- 18529789


txdb <- AnnotationDbi::loadDb("../data/misc/txdb_v19_hg19.sqlite")
gr = GenomicRanges::GRanges(seqnames = "chr19", ranges = IRanges(start, end))

p_lrrc25 <- ggplot() + theme_classic() +
  geom_alignment(
    txdb,
    which = gr,
    cds.rect.h = 0.1,
    color = "black",
    fill = "black",
    label.size = 4,
    arrow.rate = 0,
    length = unit(0.2, "cm"),
    gap.geom = 'segment'
  ) +
  ylim(c(0.75,1.25)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_vline(xintercept = 18512817, linetype=2) +
  geom_point(aes(x=18512817, y=1), shape=23, size=3, fill="purple")


##
bp_lrrc25_file="~/Google Drive/Manuscript/PU1_bloodcelltrait/data/LRRC25/blueprint.LRRC25.eQTL.txt"
ld_file="~/Google Drive/Manuscript/PU1_bloodcelltrait/data/LRRC25/rs5827412.ld.txt"
bp_lrrc25_stat <- read.table(bp_lrrc25_file, header=F)
colnames(bp_lrrc25_stat) <- c("chr", "pos", "rsid", "pval", "beta")

bp_lrrc25_stat <- bp_lrrc25_stat %>% filter(pos >= start & pos < end)

shape = ifelse(bp_lrrc25_stat$pval > 1e-3, 21, ifelse(bp_lrrc25_stat$beta > 0, 24, 25))
names(shape) = bp_lrrc25_stat$rsid

snp = 'rs5827412'
size = ifelse(bp_lrrc25_stat$rsid == snp, 3, 2)
names(size) = bp_lrrc25_stat$rsid

#color = ifelse(bp_lrrc25_stat$rsid == snp, 'purple', "#FF7373")
#names(color) = bp_lrrc25_stat$rsid

ld = read.table(ld_file, header=T)

color = assign_color2(bp_lrrc25_stat$rsid, snp, ld)
names(color) = bp_lrrc25_stat$rsid


title <- "LRRC25"
p_bp_lrrc25 <- ggplot(bp_lrrc25_stat) +
  geom_point(aes(x=pos, y=-log10(pval)), shape=shape, fill=color,size =size) +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_vline(xintercept = 18512817, linetype=2) +
  ylab(expression(paste(-log[10],"(",italic(p),")"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))


p_bp_lrrc25 <- ggplot(bp_lrrc25_stat) +
  geom_point(aes(x=pos, y=-log10(pval), shape=rsid, fill=rsid, size =rsid)) +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_vline(xintercept = 18512817, linetype=2) +
  ylab(expression(paste(-log[10],"(",italic(p),")"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1))))
p_bp_lrrc25


title <- expression(paste("Monocyte ", italic('LRRC25'), " eQTL (BLUEPRINT)"))
p_bp_lrrc25 <- p_bp_lrrc25 + labs(title = title)


p_5_f <- p_bp_lrrc25 + p_lrrc25 + plot_layout(nrow=2, heights = c(2,1))

ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Fig5_f.pdf',
                plot = p_5_c,
                device='pdf',
                width=200, height=72, units="mm")


#### SUPP

rs5827412_gwas <- data.frame(source=as.factor(c("mono %", "mono #", "neut %", "neut #", "wbc #")),
                             beta=c(-0.0451654, -0.029094, 0.0304417, 0.0262255, 0.0160937),
                             se=c(0.00216547, 0.00213564, 0.00220315, 0.00219324, 0.0021686))


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


## ALTERNATE - colocalization plot with - log 10 p values

gwas_file="~/Google Drive/Manuscript/PU1_bloodcelltrait/data/LRRC25/mono_p.LRRC25.filtered.txt"
pu1_file="~/Google Drive/Manuscript/PU1_bloodcelltrait/data/LRRC25/PU1_71930.bqtl.filtered.txt"
ld_file="~/Google Drive/Manuscript/PU1_bloodcelltrait/data/LRRC25/rs5827412.ld.txt"
snp = 'rs5827412'

gwas_stat <- read.table(gwas_file, header=F)
colnames(gwas_stat) <- c("rsid", "chr", "pos", "pval")
pu1_stat <- read.table(pu1_file, header=F)
colnames(pu1_stat) <- c("rsid", "chr", "pos", "pval")

merged_stat <- merge(gwas_stat, pu1_stat, by=c("rsid", "pos"))


ld = read.table(ld_file, header=T)

color = assign_color2(merged_stat$rsid, snp, ld)

shape = ifelse(merged_stat$rsid == snp, 23, 21)
names(shape) = merged_stat$rsid

size = ifelse(merged_stat$rsid == snp, 3, 2)
names(size) = merged_stat$rsid




p_5_b <- ggplot(merged_stat) + geom_point(aes(x=-log10(pval.x), y=-log10(pval.y), fill = rsid, size = rsid, shape = rsid)) +
  theme_minimal() +
  scale_fill_manual(values = color, guide = "none") +
  scale_shape_manual(values = shape, guide = "none") +
  scale_size_manual(values = size, guide = "none") +
  xlab(expression(paste(-log[10],"(",italic(p)["Mono %"],")"))) +
  ylab(expression(paste(-log[10],"(",italic(p)["PU.1 bQTL"],")"))) +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

p_5_b

geom_point(data=merged_stat %>% filter(rsid == snp), aes(x=-log10(pval.x), y=-log10(pval.y)), shape=23, size=3, fill="purple")

ggplot2::ggsave('~/Google Drive/Manuscript/PU1_bloodcelltrait/figures/Fig5_b.pdf',
                plot = p_5_b,
                device='pdf',
                width=110, height=100, units="mm")

ggrepel::geom_text_repel(data=merged_stat %>% filter(rsid == snp), aes(label=label, x=-log10(pval.x), y=-log10(pval.y)), max.overlaps = Inf)
