library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(cowplot, quietly = T)
library(data.table, quietly = T)
library(AnnotationDbi, quietly = T)
library(ggbio, quietly = T)
library(GenomicRanges, quietly = T)
library(locuscomparer, quietly = T)
library(patchwork, quietly = T)

source("misc/locuscompare_updated.R")


## Figure 5a - Colocalization plot

gwas_file="../data/Fig5/coloc/mono_p.LRRC25.withbeta.filtered.txt"
pu1_file="../data/Fig5/coloc/PU1_71930.bqtl.withbeta.filtered.txt"
ld_file="../data/Fig5/coloc/rs5827412.ld.txt"
snp = 'rs5827412'

gwas_stat <- read.table(gwas_file, header=T)
#colnames(gwas_stat) <- c("rsid", "chr", "pos", "pval", "beta", "se")
gwas_stat$z <- gwas_stat$beta / gwas_stat$se

pu1_stat <- read.table(pu1_file, header=T)
#colnames(pu1_stat) <- c("rsid", "chr", "pos", "pval", "beta")
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
#rs5827412_MPRA <- data.frame(source=as.factor(c("Tewhey2016", "Abell2022")), logSkew=c(-0.390574802, -0.342679482),logSkewSE=c(0.09816138,0.307302851))

rs5827412_mpra <- read.csv('../data/Fig5/rs5827412_mpra.txt', header=T, sep='\t')

p_rs5827412_mpra <- ggplot(rs5827412_mpra, aes(x= source, y = logSkew, ymin=logSkew-1.96*logSkewSE, ymax=logSkew+1.96*logSkewSE)) +
  geom_col(fill="gray", color="black") +
  theme_classic() + scale_x_discrete(position='top') +
  geom_errorbar(color = "black", width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(-1,0.3)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=14),
        aspect.ratio = 1.1) +
  geom_hline(yintercept=0, size=1) + ylab("rs5827412 \n allelic skew")


### Figure 5d - PU1KO
pu1ko_atac_lrrc25 <- read.csv('../data/Fig5/pu1_ko_atac_lrrc25.txt', header = T, sep ='\t')

p_pu1ko <- pu1ko_atac_lrrc25 %>% mutate(condition = factor(condition, levels=c("SPI1+/+", "SPI-/-"))) %>%
  ggplot(aes(x=as.factor(condition), y=cpm)) +
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

p_lrrc25_rna_mono_diff <- LRRC25_RNA %>% filter(celltype %in% c("HSC", "MPP", "CMP", "GMP", "Mono")) %>%
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


p_5_cde <- p_rs5827412_mpra + p_pu1ko + p_lrrc25_rna_mono_diff + plot_layout(ncol=3, widths = c(1,1,1.4))


ggplot2::ggsave('../figures/Fig5_cde.pdf',
                plot = p_5_cde,
                device='pdf',
                width=250, height=75, units="mm")




#### Figure 5f,g - LRRC25 eQTL (monocyte) and gene plot

# Common coordinate for 5f,g
chr <- "19"
start <- 18495845
end <- 18529789


## LRRC25 eQTL plot
bp_lrrc25_file="../data/Fig5/blueprint.LRRC25.eQTL.txt"
bp_lrrc25_stat <- read.csv(bp_lrrc25_file, header=T, sep='\t')
#colnames(bp_lrrc25_stat) <- c("chr", "pos", "rsid", "pval", "beta")

bp_lrrc25_stat <- bp_lrrc25_stat %>% filter(pos >= start & pos < end)

## Defining shape depending on estimated effect direction
shape = ifelse(bp_lrrc25_stat$pval > 1e-3, 21, ifelse(bp_lrrc25_stat$beta > 0, 24, 25))
names(shape) = bp_lrrc25_stat$rsid

snp = 'rs5827412'
size = ifelse(bp_lrrc25_stat$rsid == snp, 3, 2)
names(size) = bp_lrrc25_stat$rsid


p_bp_lrrc25 <- ggplot(bp_lrrc25_stat) +
  geom_point(aes(x=pos, y=-log10(pval)), shape=shape, fill="#FF7373",size=size) +
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
        axis.title.y = element_text(size=14))



### Figure 5g - LRRC25 Gene plot
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



p_5_fg <- p_bp_lrrc25 + p_lrrc25 + plot_layout(nrow=2, heights = c(2,1))

ggplot2::ggsave('../figures/Fig5_fg.pdf',
                plot = p_5_fg,
                device='pdf',
                width=200, height=72, units="mm")
