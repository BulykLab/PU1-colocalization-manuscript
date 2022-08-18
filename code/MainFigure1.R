library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(stringr)
library(patchwork)

## Figure 1b  -- Enrichment of # blood cell trait association tagging variants
tag_enrich <- read.csv('../data/Fig1/PU1_100kb_tagging_enrichment.txt', header=T, sep='\t')
#tag_enrich$P <- ifelse(tag_enrich$P == 0, 0.001, tag_enrich$P)
tag_enrich$label <- rep("", length(tag_enrich[1]))

tag_enrich$P <- (tag_enrich$P * 250 +1) / 251


for (i in 1:length(tag_enrich[,1])) {
  temp <- as.vector(tag_enrich$Trait)
  if (tag_enrich[i,2] < 0.05) {
    tag_enrich[i,6] <- temp[i]
  }
}


p_1b <- tag_enrich %>% mutate(Category = factor(Category, levels=c("Granulocyte", "Monocyte", "Lymphocyte", "MatureRed", "ImmatureRed","Platelet"))) %>%
  ggplot(aes(x=obs/exp, y=-log10(P), color=Category)) +
  geom_hline(yintercept = -log10(0.05), linetype=2, alpha = 0.4) +
  geom_point(size=3) +
  theme_classic() +
  xlim(c(0,3.5)) + ylim(c(0,2.5)) +
  scale_color_manual(values=c('#F98400', '#351C75', '#40A2EB', '#FFB6C1', '#A61C00', '#008000'), name="Trait group") +
  xlab("Fold enrichment") + ylab(expression(paste(-log[10],"(",italic(p),")"))) +
  ggrepel::geom_text_repel(aes(label=label), min.segment.length = 0, box.padding = 0.7, ylim=c(1.4,2.5), seed=1) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=16),
        legend.title=element_text(size=12), legend.text=element_text(size=11),
        legend.position = c(0.2, 0.78),
        legend.background = element_rect(fill = "white", color = "black"),
        aspect.ratio=1)


## Figure 1c - gkmSVM vs per-allele effect size

pu1_gkmsvm <- read.table('../data/Fig1/PU1.pwm_gkmsvm.nominal.snps.txt')
pu1_gkmsvm <- filter(pu1_gkmsvm, (V13 > 0 | V14 > 0) & V13 != V14)


# gkmsvm
p_1c <- ggplot(pu1_gkmsvm, aes(x=V14-V13, y=V17, color=-log10(V18))) +
  geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
  geom_point(alpha=0.8) +
  labs(y="PU.1 bQTL effect size", x = expression(Delta~"gkmSVM") , color=expression(paste(-log[10],"(",italic(p),")"))) +
  scale_colour_gradient(high = "black", low = "#E9E9E9") +
  theme_classic() + ylim(c(-1.9,1.9)) + xlim(c(-9,9)) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text=element_text(size=16),
        legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = "white", color = "black"),
        aspect.ratio=1) +
  geom_point(data=pu1_gkmsvm %>% filter(V4=='PU1_21900'), aes(x=V14-V13, y=V17), shape=23, size = 3, fill='purple') +
  geom_point(data=pu1_gkmsvm %>% filter(V4=='PU1_53081'), aes(x=V14-V13, y=V17), shape=23, size = 3, fill='purple') +
  geom_point(data=pu1_gkmsvm %>% filter(V4=='PU1_70939'), aes(x=V14-V13, y=V17), shape=23, size = 3, fill='purple')







## Figure 1d  --  Number of motif-altering variants at the center of peak

## number of motif-altering variants by location
#location_number <- data.frame(location = rep(c('Center 50bp', 'Elsewhere'),2),
#                              pu1motif = c(rep('motif',2), rep('none', 2)),
#                              number = c(316,38,453,241))

motif_location <- read.csv('../data/Fig1/motif_location.txt', header=T, sep='\t')

p_1d <- motif_location %>% mutate(pu1motif = factor(pu1motif, levels=c("none","motif"))) %>%
ggplot(aes(x=location, y=number)) +
  geom_col(aes(fill = location, alpha=pu1motif, linetype=pu1motif), width = 0.8, color = 'black') +
  theme_classic() +
  ylab(expression(paste("Number of PU.1 bQTLs"))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        aspect.ratio = 2.2) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values =c('#F8766D', '#B8B4B4')) +
  scale_alpha_manual(values = c(0,1),labels = c('no','yes'), name = "PU.1 motif altered") +
  scale_linetype_manual(values = c('dashed', 'solid'), labels = c('no', 'yes'), name = expression(paste(italic("PU.1")," motif altered"))) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  guides(fill = "none", alpha = "none", linetype="none")


p_1_merged <- p_1b + p_1c + p_1d + plot_layout(ncol = 3, widths = c(1, 1, 0.4))


ggsave('../figures/Fig1.pdf',
       plot = p_1_merged,
       device='pdf',
       width=360, height=130, units="mm")
