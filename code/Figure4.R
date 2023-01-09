library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)

source("misc/locuscompare_updated.R")

chr = "5"
start <- 123952000
end <- 124520000

txdb <- AnnotationDbi::loadDb("../data/misc/txdb_v19_hg19.sqlite")

gr = GenomicRanges::GRanges(seqnames = "chr5", ranges = IRanges(123800000, 125000000))


##
# Size of red curves determined by B cell PCHiC data (Javierre et al. 2016)
p_znf608 <- ggplot() + theme_classic() +
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
  ylim(c(0.75,1.75)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2) +
  geom_curve(aes(x = 124081118, y = 1, xend = 124341259, yend = 1), size=2.297, curvature = -0.16, color = "#FF6379") +
  geom_curve(aes(x = 124081118, y = 1, xend = 124285447, yend = 1), size=0.963, curvature = -0.1, color = "#FF6379") +
  geom_point(aes(x = 124341259, y=1), shape=23, size=3, fill="purple") +
  geom_point(aes(x = 124285447, y=1), shape=23, size=3, fill="yellow")



##### Start Here
snp = NULL
chr='5'
population = 'EUR'


## lym count GWAS
gwas_file = "../data/Fig4/ZNF608.lym_count.metal.txt"
gwas_stat = read_metal(gwas_file, marker_col = 'rsid', pval_col = 'pval')
gwas_stat = get_position(gwas_stat)
snp = 'rs12517864'
ld = retrieve_LD(chr, snp, population)
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
p_gwas <- p_gwas + geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2)+
  geom_point(aes(x=124341259, y=gwas_stat[gwas_stat$pos ==124341259,]$logp), shape=23, size=3, fill="purple") +
  geom_point(aes(x=124285447, y=gwas_stat[gwas_stat$pos ==124285447,]$logp), shape=23, size=3, fill="yellow")


## PU1_27777
pu1_file = "../data/Fig4/PU1_27777.metal.txt"
pu1_stat = read_metal(pu1_file, marker_col = 'rsid', pval_col = 'pval')
pu1_stat = get_position(pu1_stat)
snp = 'rs12517864'
ld = retrieve_LD(chr, snp, population)
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
#title = expression(paste(italic('PU.1'), " bQTL"))
p_pu1 <- p_pu1 + geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2)+
  geom_point(aes(x=124341259, y=pu1_stat[pu1_stat$pos ==124341259,]$logp), shape=23, size=3, fill="purple") +
  geom_point(aes(x=124285447, y=pu1_stat[pu1_stat$pos ==124285447,]$logp), shape=23, size=3, fill="yellow") + labs(title = title)


## ZNF608
eqtl_file = "../data/Fig4/ZNF608.LCL_eqtl.metal.txt"
eqtl_stat = read_metal(eqtl_file, marker_col = 'rsid', pval_col = 'pval')
eqtl_stat = get_position(eqtl_stat)
#snp = 'rs12517864'
snp = 'rs2028854'
ld = retrieve_LD(chr, snp, population)
color = assign_color2(eqtl_stat$rsid, snp, ld)

shape = ifelse(eqtl_stat$rsid == snp, 23, 21)
names(shape) = eqtl_stat$rsid

size = ifelse(eqtl_stat$rsid == snp, 3, 2)
names(size) = eqtl_stat$rsid

eqtl_stat$label = ifelse((eqtl_stat$rsid == 'rs2028854'), eqtl_stat$rsid, '')
#eqtl_stat$label = ifelse((eqtl_stat$rsid == snp | eqtl_stat$rsid == 'rs2028854'), eqtl_stat$rsid, '')

metal = eqtl_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
#colnames(metal)[which(colnames(metal) == 'logp1')] = 'logp'
title = "ZNF608 eQTL"
p_eqtl <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)
title <- expression(paste(italic('ZNF608'), " eQTL (LCLs)"))
#title <- expression(paste(italic('ZNF608'), " eQTL"))
p_eqtl <- p_eqtl + geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2)+
  geom_point(aes(x=124285447, y=eqtl_stat[eqtl_stat$pos ==124285447,]$logp), shape=23, size=3, fill="yellow")+
  geom_point(aes(x=124341259, y=eqtl_stat[eqtl_stat$pos ==124341259,]$logp), shape=23, size=3, fill="purple") + labs(title = title)


## ZNF608 conditioned on rs2028854
eqtl_cond_file = "../data/Fig4/ZNF608.LCL_cond_eqtl.metal.txt"
eqtl_cond_stat = read_metal(eqtl_cond_file, marker_col = 'rsid', pval_col = 'pval')
eqtl_cond_stat = get_position(eqtl_cond_stat)
snp = 'rs12517864'
ld = retrieve_LD(chr, snp, population)
color = assign_color2(eqtl_cond_stat$rsid, snp, ld)

shape = ifelse(eqtl_cond_stat$rsid == snp, 23, 21)
names(shape) = eqtl_cond_stat$rsid

size = ifelse(eqtl_cond_stat$rsid == snp, 3, 2)
names(size) = eqtl_cond_stat$rsid

eqtl_cond_stat$label = ifelse(eqtl_cond_stat$rsid == snp, eqtl_cond_stat$rsid, '')


metal = eqtl_cond_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
#colnames(metal)[which(colnames(metal) == 'logp1')] = 'logp'
title <- "ZNF608 eQTL (LCLs) conditioned on rs2028854"
p_eqtl_cond <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)
title <- expression(paste(italic('ZNF608'), " eQTL (LCLs) conditioned on rs2028854"))
p_eqtl_cond <- p_eqtl_cond + geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2)+
  geom_point(aes(x=124341259, y=eqtl_cond_stat[eqtl_cond_stat$pos ==124341259,]$logp), shape=23, size=3, fill="purple") +
  geom_point(aes(x=124285447, y=eqtl_cond_stat[eqtl_cond_stat$pos ==124285447,]$logp), shape=23, size=3, fill="yellow") + labs(title = title)

## SUSIE fine-mapping
znf608_susie <- read.csv('../data/Fig4/ZNF608.LCL_eqtl.susie.txt', header = T, sep ='\t')
znf608_susie$credible_set <- as.factor(znf608_susie$credible_set)  # credible set id
znf608_susie$color <- ifelse(znf608_susie$credible_set == 1, "#00AFBB", "#FC4E07")

title = "ZNF608 eQTL credible sets (SuSiE)"

p_fm <- ggplot(znf608_susie, aes(x=pos, y=pip, fill=credible_set)) +
  geom_point(alpha=1, shape=21, size=2) + theme_classic() +
  labs(y="PIP", x = "chr5", fill="Credible set ID", title = title) +
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  scale_x_continuous(lim = c(start, end), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(lim = c(0,1), breaks = seq(0,1, 0.5)) +
  theme(legend.position = c(0.9, 0.55),
        legend.title=element_text(size=11), legend.text = element_text(size=11), # legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(
          hjust = 1e-2,
          margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1)))
        ) +
  geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2) +
  geom_point(aes(x=124341259, y=znf608_susie[znf608_susie$pos ==124341259,]$pip), shape=23, size=3, fill="purple") +
  geom_point(aes(x=124285447, y=znf608_susie[znf608_susie$pos ==124285447,]$pip), shape=23, size=3, fill="yellow")

title <- expression(paste(italic('ZNF608'), " eQTL (LCLs) credible sets (SuSiE)"))
p_fm <- p_fm + labs(title = title)


## B cells
bcell_eqtl_file = "../data/Fig4/ZNF608.Bcell_naive_eqtl.metal.txt"
bcell_eqtl_stat = read_metal(bcell_eqtl_file, marker_col = 'rsid', pval_col = 'pval')
bcell_eqtl_stat = get_position(bcell_eqtl_stat)
snp = 'rs12517864'
ld = retrieve_LD(chr, snp, population)
color = assign_color2(bcell_eqtl_stat$rsid, snp, ld)

shape = ifelse(bcell_eqtl_stat$rsid == snp, 23, 21)
names(shape) = bcell_eqtl_stat$rsid

size = ifelse(bcell_eqtl_stat$rsid == snp, 3, 2)
names(size) = bcell_eqtl_stat$rsid

bcell_eqtl_stat$label = ifelse((bcell_eqtl_stat$rsid == 'rs12517864'), bcell_eqtl_stat$rsid, '')

metal = bcell_eqtl_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "ZNF608 eQTL"
p_bcell_eqtl <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)
title <- expression(paste(italic('ZNF608'), " eQTL (Naive B cells)"))
p_bcell_eqtl <- p_bcell_eqtl + geom_vline(xintercept = 124341259, linetype=2) + geom_vline(xintercept = 124285447, linetype=2)+
  geom_point(aes(x=124285447, y=bcell_eqtl_stat[bcell_eqtl_stat$pos ==124285447,]$logp), shape=23, size=3, fill="yellow") +
  geom_point(aes(x=124341259, y=bcell_eqtl_stat[bcell_eqtl_stat$pos ==124341259,]$logp), shape=23, size=3, fill="purple") + labs(title = title)






#### Plotting together
p_4_acdeg <- p_pu1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + scale_y_continuous(breaks = c(0,5,10)) +
  p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + scale_y_continuous(expand = expansion(mult=c(0.02,0.1)), breaks = c(0,5,10)) +
  p_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + scale_y_continuous(expand = expansion(mult=c(0.02,0.15)),breaks = c(0,5,10)) +
  p_eqtl_cond +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + scale_y_continuous(expand = expansion(mult=c(0.02,0.4))) +
  p_fm + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14), legend.key.height = unit(0.3,'cm')) +
  p_bcell_eqtl +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + scale_y_continuous(expand = expansion(mult=c(0.02,0.1))) +
  p_znf608 +
  plot_layout(nrow = 7, heights = c(3, 3, 3, 3, 3, 3, 2))


ggplot2::ggsave('../figures/Fig4acdeg.pdf',
       plot = p_4_acdeg,
       device='pdf',
       width=200, height=200, units="mm")


## Figure h
### Regulatory QTLs

## PU.1 peak

PU1_27777 <- read.csv('../data/Fig4/qtl/PU1_27777.cpm.boxplot.txt', header = T, sep ='\t')

# Adding allele dosage 2 because there were no samples with homozygous alternate alleles
p_pu1_qtl_4 <- PU1_27777 %>% mutate(genotype = factor(genotype, levels=c("0", "1", "2"))) %>%
ggplot(aes(x=genotype, y=cpm)) +
 geom_boxplot(fill='#46ACC8', color="black") +
 theme_classic() +
 geom_jitter(shape=16, position=position_jitter(0.1)) +
 labs(y="Read per million", title = "PU.1 ChIP") +
 ylim(0, NA) +
 scale_x_discrete(breaks=factor(0:2), drop = F) +
 theme(plot.title = element_text(size=14),
       axis.title.x = element_blank(),
       axis.title.y = element_text(size=12),
       axis.text.x = element_blank(),
       axis.text.y = element_text(size=12))

## ATAC peak over PU.1 peak

ATAC_59004 <- read.csv('../data/Fig4/qtl/ATAC_59004.cpm.boxplot.txt', header = T, sep ='\t')

p_atac_qtl_4 <- ggplot(ATAC_59004, aes(x=as.factor(genotype), y=cpm)) +
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

H3K4me1_5 <- read.csv('../data/Fig4/qtl/H3K4me1_5_124340065_124345334.cpm.boxplot.txt', header = T, sep ='\t')

p_h3k4me1_qtl_4 <- ggplot(H3K4me1_5, aes(x=as.factor(genotype), y=cpm)) +
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

## H3K27ac

H3K27ac_5 <- read.csv('../data/Fig4/qtl/H3K27ac_5_124339830_124348198.cpm.boxplot.txt', header = T, sep ='\t')

p_h3k27ac_qtl_4 <- ggplot(H3K27ac_5, aes(x=as.factor(genotype), y=cpm)) +
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


p_4_h <- p_pu1_qtl_4 +
 p_atac_qtl_4 +
 p_h3k4me1_qtl_4 +
 p_h3k27ac_qtl_4 +
 plot_layout(nrow = 2, ncol = 2, widths = c(1,1), heights = c(1,1))

## Figure 4i  - ZNF608 expression across blood cell types (related to lymphocytes)
heme_RNA <- read.csv('~/Projects/PU1_gwas/plots/heme_RNA.txt', header = T, sep ='\t', row.names = 1)
ZNF608_RNA <- transpose((heme_RNA / colSums(heme_RNA) * 10^6)["ZNF608",1:49])
celltype <- c("HSC", "HSC", "HSC", "HSC", "MPP", "MPP", "MPP", "MPP", "LMPP", "LMPP", "LMPP", "CMP", "CMP", "CMP", "CMP", "GMP", "GMP", "GMP", "GMP", "MEP", "MEP", "MEP", "MEP", "Mono", "Mono", "Mono", "Mono", "CD4T", "CD4T", "CD4T", "CD4T", "CD8T", "CD8T", "CD8T", "CD8T", "NK", "NK", "NK", "NK", "B", "B", "B", "B", "CLP", "CLP", "CLP", "Ery", "Ery", "Ery")
ZNF608_RNA["celltype"] <- celltype

p_4_i <- ZNF608_RNA %>% filter(celltype %in% c("HSC", "MPP", "LMPP", "CLP",  "B", "CD4T", "CD8T", "NK")) %>%
 mutate(celltype = factor(celltype, levels=c("NK", "CD8T",  "CD4T", "B", "CLP", "LMPP", "MPP", "HSC"))) %>%
 ggplot(aes(x=celltype, y=(V1+0.1), fill=celltype)) +
 geom_boxplot(color="black") + geom_jitter(shape=16, position=position_jitter(0.1)) +
 theme_minimal() +
 labs(y=expression(paste(italic('ZNF608'), " CPM")), x="") +
 scale_fill_manual(values = c("#FF7373", "#FF7373", "#FF7373", "#FF7373", "#E56767", "#B25050", "#7F3939", "#333333")) +
 theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=14),
       axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
       aspect.ratio = 1.4,
       legend.position = "none") +
 coord_flip()





## Plotting Fig4h,i together
p_4_hi <- cowplot::plot_grid(p_4_h, p_4_i, nrow = 2, rel_heights = c(5,6))

ggplot2::ggsave('../figures/Fig4hi.pdf',
                plot = p_4_hi,
                device='pdf',
                width=105, height=180, units="mm")
