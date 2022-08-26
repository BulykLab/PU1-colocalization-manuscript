library(dplyr)
library(ggplot2)

## Supp Fig 1 PU.1 expression across blood cell types
heme_RNA <- read.csv('../data/misc/heme_RNA.txt', header = T, sep ='\t', row.names = 1)

SPI1_RNA <- as.data.frame(t((heme_RNA / colSums(heme_RNA) * 10^6)["SPI1",1:49]))

celltype <- c("HSC", "HSC", "HSC", "HSC", "MPP", "MPP", "MPP", "MPP", "LMPP", "LMPP",
              "LMPP", "CMP", "CMP", "CMP", "CMP", "GMP", "GMP", "GMP", "GMP", "MEP",
              "MEP", "MEP", "MEP", "Mono", "Mono", "Mono", "Mono", "CD4T", "CD4T",
              "CD4T", "CD4T", "CD8T", "CD8T", "CD8T", "CD8T", "NK", "NK", "NK", "NK",
              "B", "B", "B", "B", "CLP", "CLP", "CLP", "Ery", "Ery", "Ery")

progenitor <- c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP")

SPI1_RNA["celltype"] <- celltype
SPI1_RNA$SPI1 = ifelse(SPI1_RNA$SPI1 == 0, 0.1, SPI1_RNA$SPI1)

p_pu1_RNA <- SPI1_RNA %>% mutate(celltype = factor(celltype, levels=c("HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "Mono", "B", "NK", "CD4T", "CD8T","Ery"))) %>%
  ggplot(aes(x=celltype, y=SPI1, fill=!(celltype %in% progenitor))) +
  geom_boxplot(color="black") + theme_classic() + geom_jitter(shape=16, position=position_jitter(0.1)) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values =c( '#B8B4B4', '#F8766D'), labels=c("Progenitor", "Differentiated")) +
  labs(y=expression(paste(italic("PU.1")," count per million")), x = "Cell type") +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=16), legend.position = 'bottom')


ggsave('../figures/SuppFig1.pdf',
       plot = p_pu1_RNA,
       device='pdf',
       width=200, height=150, units="mm")
