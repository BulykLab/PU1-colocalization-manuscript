library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(cowplot, quietly = T)
library(data.table, quietly = T)
library(patchwork, quietly = T)


#### Figure S8 - LRRC25 expression across all blood cell types

heme_RNA <- read.csv('../data/misc/heme_RNA.txt', header = T, sep ='\t', row.names = 1)
LRRC25_RNA <- transpose((heme_RNA / colSums(heme_RNA) * 10^6)["LRRC25",1:49])
celltype <- c("HSC", "HSC", "HSC", "HSC", "MPP", "MPP", "MPP", "MPP", "LMPP", "LMPP", "LMPP", "CMP", "CMP", "CMP", "CMP", "GMP", "GMP", "GMP", "GMP", "MEP", "MEP", "MEP", "MEP", "Mono", "Mono", "Mono", "Mono", "CD4T", "CD4T", "CD4T", "CD4T", "CD8T", "CD8T", "CD8T", "CD8T", "NK", "NK", "NK", "NK", "B", "B", "B", "B", "CLP", "CLP", "CLP", "Ery", "Ery", "Ery")
LRRC25_RNA["celltype"] <- celltype

p_lrrc25_rna_all <- ggplot(LRRC25_RNA, aes(x=reorder(as.factor(celltype), -V1), y=(V1+0.1), fill=(celltype=='Mono'))) +
  geom_boxplot(color="black") + geom_jitter(shape=16, position=position_jitter(0.1)) +
  theme_classic() + scale_y_continuous(trans='log10') +
  scale_fill_manual(values = c("#B8B4B4", "#FF7373")) +
  labs(y=expression(paste(italic('LRRC25'), " count per million")), x = "Cell type") +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position="none")

ggplot2::ggsave('../figures/FigS8.pdf',
                plot = p_lrrc25_rna_all,
                device='pdf',
                width=250, height=150, units="mm")

