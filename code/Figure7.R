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


### Figure 7b-d association plots
## Gene plot
chr='8'
population = 'EUR'

start <- 79420000
end <- 79730000

txdb <- AnnotationDbi::loadDb("../data/misc/txdb_v19_hg19.sqlite")
gr = GenomicRanges::GRanges(seqnames = "chr8", ranges = IRanges(start, end))

p_zc2hc1a <- ggplot() + theme_classic() +
  geom_alignment(
    txdb,
    which = gr,
    cds.rect.h = 0.1,
    color = "black",
    fill = "black",
    label.size = 5,
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
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_vline(xintercept = 79577726, linetype=2) +
  geom_point(aes(x=79577726, y=1), shape=23, size=3, fill="purple")



## PU1_40678
pu1_file = "../data/Fig7/PU1_40678.pu1bqtl.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')
pu1_stat = get_position(pu1_stat)
snp = 'rs3808619'
ld = retrieve_LD(chr, snp, population)
#color = assign_color(pu1_stat$rsid, snp, ld)
color = assign_color2(pu1_stat$rsid, snp, ld)

shape = ifelse(pu1_stat$rsid == snp, 23, 21)
names(shape) = pu1_stat$rsid

size = ifelse(pu1_stat$rsid == snp, 3, 2)
names(size) = pu1_stat$rsid

pu1_stat$label = ifelse(pu1_stat$rsid == snp, pu1_stat$rsid, '')

metal = pu1_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
#colnames(metal)[which(colnames(metal) == 'logp1')] = 'logp'
title = "PU.1 bQTL"
p_pu1 <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)
p_pu1 <- p_pu1 + geom_vline(xintercept = 79577726, linetype=2) +
  geom_point(aes(x=79577726, y=pu1_stat[pu1_stat$pos ==79577726,]$logp), shape=23, size=3, fill="purple") +
  labs(title = title)

## lym count GWAS
gwas_file = "../data/Fig7/ZC2HC1A.lym_count.metal.txt"
gwas_stat = read_metal(gwas_file, marker_col = 'rsid', pval_col = 'pval')
gwas_stat = get_position(gwas_stat)
snp = 'rs3808619'
ld = retrieve_LD(chr, snp, population)
#color = assign_color(gwas_stat$rsid, snp, ld)
color = assign_color2(gwas_stat$rsid, snp, ld)


shape = ifelse(gwas_stat$rsid == snp, 23, 21)
names(shape) = gwas_stat$rsid

size = ifelse(gwas_stat$rsid == snp, 3, 2)
names(size) = gwas_stat$rsid

gwas_stat$label = ifelse(gwas_stat$rsid == snp, gwas_stat$rsid, '')

metal = gwas_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
#colnames(metal)[which(colnames(metal) == 'logp1')] = 'logp'
title = "Lymphocyte count GWAS"
p_gwas <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)
p_gwas <- p_gwas + geom_vline(xintercept = 79577726, linetype=2) +
  geom_point(aes(x=79577726, y=gwas_stat[gwas_stat$pos ==79577726,]$logp), shape=23, size=3, fill="purple")


## SUSIE fine-mapping - lym_count
lymcount_susie <- read.csv('../data/Fig7/ZC2HC1A.Lym.CS.bed', header = T, sep ='\t')

title = "Lymphocyte count credible set (SuSiE)"

p_fm <- ggplot(lymcount_susie, aes(x=end, y=pip)) +
  geom_point(alpha=1, shape=21, size=2, fill='#00AFBB') + theme_classic() +
  labs(y="PIP", x = "chr8", fill="Credible set ID", title = title) +
  scale_x_continuous(lim = c(start, end), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(lim = c(0,1), breaks = seq(0,1, 0.5)) +
  theme(legend.position = c(0.9, 0.55),
        legend.title=element_text(size=11), legend.text = element_text(size=11),
        legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1)))
  ) +
  geom_vline(xintercept = 79577726, linetype=2) +
  geom_point(aes(x=79577726, y=lymcount_susie[lymcount_susie$end ==79577726,]$pip), shape=23, size=3, fill="purple")


title <- expression(paste("Lymphocyte count 95% credible set (SuSiE) (", italic('n'), "=44)"))
p_fm <- p_fm + labs(title = title)

### TOGETHER

p_7_bcd <- p_pu1 + scale_y_continuous(expand = expansion(mult = c(0, 0.08))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  p_gwas + scale_y_continuous(expand = expansion(mult = c(0, 0.08))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  p_fm + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  p_zc2hc1a +
  plot_layout(nrow = 4, heights = c(1, 1, 0.4, 0.4))


ggplot2::ggsave('../figures/Fig7bcd.pdf',
                plot = p_7_bcd,
                device='pdf',
                width=250, height=125, units="mm")



### Figure 7e QTLs
## PU.1 peak

PU1_40678 <- read.csv('../data/Fig7/qtl/PU1_40678.cpm.boxplot.txt', header = T, sep ='\t')

p_pu1_qtl_7 <- ggplot(PU1_40678, aes(x=as.factor(genotype), y=cpm)) +
  geom_boxplot(fill='#46ACC8', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(y="Read per million", title = "PU.1 ChIP") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))

## ATAC peak over PU.1 peak

ATAC_87956 <- read.csv('../data/Fig7/qtl/ATAC_87956.cpm.boxplot.txt', header = T, sep ='\t')

p_atac_qtl_7 <- ggplot(ATAC_87956, aes(x=as.factor(genotype), y=cpm)) +
  geom_boxplot(fill='#5166CC', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "ATAC") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))

## H3K4me1

H3K4me1_8_7 <- read.csv('../data/Fig7/qtl/H3K4me1_8_79576449_79578176.cpm.boxplot.txt', header = T, sep ='\t')

p_h3k4me1_qtl_7 <- ggplot(H3K4me1_8_7, aes(x=as.factor(genotype), y=cpm)) +
  geom_boxplot(fill='#E7B800', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "H3K4me1") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))

## H3K4me3

H3K4me3_8_7 <- read.csv('../data/Fig7/qtl/H3K4me3_8_79576923_79580028.cpm.boxplot.txt', header = T, sep ='\t')

p_h3k4me3_qtl_7 <- ggplot(H3K4me3_8_7, aes(x=as.factor(genotype), y=cpm)) +
  geom_boxplot(fill='#FC4E07', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "H3K4me3") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))

## H3K27ac

H3K27ac_8_7 <- read.csv('../data/Fig7/qtl/H3K27ac_8_79576467_79579592.cpm.boxplot.txt', header = T, sep ='\t')

p_h3k27ac_qtl_7 <- ggplot(H3K27ac_8_7, aes(x=as.factor(genotype), y=cpm)) +
  geom_boxplot(fill='#009E73', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="Read per million", title = "H3K27ac") +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))


## ZC2HC1A

ZC2HC1A <- read.csv('../data/Fig7/qtl/ZC2HC1A.rpkm.boxplot.txt', header = T, sep ='\t')

p_zc2hc1a_qtl_7 <- ggplot(ZC2HC1A, aes(x=as.factor(genotype), y=rpkm)) +
  geom_boxplot(fill='darkgray', color="black") +
  theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(y="RPKM", title = expression(italic("ZC2HC1A"))) +
  ylim(0, NA) +
  theme(plot.title = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))


p_7_e_vertical <- p_pu1_qtl_7 +
  p_atac_qtl_7 +
  p_h3k4me1_qtl_7 +
  p_h3k4me3_qtl_7 +
  p_h3k27ac_qtl_7 +
  p_zc2hc1a_qtl_7 +
  plot_layout(nrow = 6)


ggplot2::ggsave('../figures/Fig7e.pdf',
                plot = p_7_e_vertical,
                device='pdf',
                width=60, height=220, units="mm")


### Figure 7f,g - MPRA

ZC2HC1A_MPRA <- read.table(file = "../data/Fig7/Abell2021.ZC2HC1A.MPRA.txt", header = TRUE)
ZC2HC1A_MPRA["status"] <- ifelse(ZC2HC1A_MPRA$padj_allele >= 0.05, "NotSignificant", ifelse(ZC2HC1A_MPRA$eQTL_beta * ZC2HC1A_MPRA$log2FoldChange_allele >0, "Significant_concordant", "Significant_discordant"))

p_zc2hc1a_mpra <- ggplot(ZC2HC1A_MPRA, aes(x= log2FoldChange_allele, y = logpadj_allele, xmin=log2FoldChange_allele-lfcSE_allele, xmax=log2FoldChange_allele+lfcSE_allele, shape=direction,fill=status, color=status)) +
  geom_point(size=3) +
  theme_classic() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') +
  scale_shape_manual(name = "eQTL effect direction", values=c(25,24)) +
  scale_color_manual(values=c('gray','red', 'black'), name="Significance",
                     labels=c("Not Significant", "Significant & Concordant", "Significant & Discordant")) +
  scale_fill_manual(values=c('gray','red','black'), name="Significance",
                    labels=c("Not Significant", "Significant & Concordant", "Significant & Discordant")) +
  xlab("MPRA Allelic Effect") + ylab(expression(paste(-log[10],"(",italic(p)[adjusted],")"))) +
  xlim(c(-5.2,5.2)) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title=element_text(size=14), legend.text=element_text(size=14),
        aspect.ratio = 1)



## Figure 7f - inset

rs3808619_MPRA <- read.csv('../data/Fig7/rs3808619_mpra.txt', header = T, sep = '\t')

#data.frame(source=as.factor(c("Tewhey2016", "Abell2022")), logSkew=c(0.441827591, 1.366886623),logSkewSE=c(0.1612375,0.33971268))


p_rs3808619_mpra <- ggplot(rs3808619_MPRA, aes(x= source, y = logSkew, ymin=logSkew-1.96*logSkewSE, ymax=logSkew+1.96*logSkewSE)) +
  geom_col(fill="red") + geom_errorbar(color = "black", width = 0.2) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        aspect.ratio = 0.7) +
  xlab("") + ylab("rs3808619 allelic skew") +
  scale_y_continuous(expand = expansion(mult = c(0, .05)))




### Figure 7g PU1 KO
pu1_ko_atac_zc2hc1a <- read.csv('../data/Fig7/pu1_ko_atac_zc2hc1a.txt', header = T, sep = '\t')

p_pu1ko <- pu1_ko_atac_zc2hc1a %>% mutate(condition = factor(condition, levels=c("SPI1 +/+", "SPI1 -/-"))) %>%
  ggplot(aes(x=as.factor(condition), y=cpm)) +
  geom_boxplot(fill='white', color="black") + theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.1), size=3) +
  labs(y="Accessibility (CPM)") + ylim(0, NA) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        aspect.ratio = 0.8)


p_7_fg <- p_zc2hc1a_mpra + p_rs3808619_mpra + p_pu1ko + plot_layout(nrow = 3, heights = c(2,1,1))


ggplot2::ggsave('../figures/Fig7fg.pdf',
                plot = p_7_fg,
                device='pdf',
                width=180, height=210, units="mm")
