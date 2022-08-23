library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(forcats, quietly = TRUE)
library(patchwork, quietly = TRUE)



## Figure 1b  -- Enrichment of # blood cell trait association tagging variants
tag_enrich <- read.csv('../data/Fig1/PU1_100kb_tagging_enrichment.txt', header=T, sep='\t')


p_2a <- tag_enrich %>% mutate(Category = factor(Category, levels=c("Granulocyte", "Monocyte", "Lymphocyte", "MatureRed", "ImmatureRed","Platelet"))) %>%
  ggplot(aes(x=obs/exp, y=-log10(qval), color=Category)) +
  geom_hline(yintercept = -log10(0.05), linetype=2, alpha = 0.4) +
  geom_point(size=3) +
  theme_classic() +
  xlim(c(0,3.5)) + ylim(c(0,2)) +
  scale_color_manual(values=c('#F98400', '#351C75', '#40A2EB', '#FFB6C1', '#A61C00', '#008000'), name="Trait group") +
  xlab("Fold enrichment") + ylab(expression(paste(-log[10],"(",italic(p)[adjusted],")"))) +
  ggrepel::geom_text_repel(aes(label=label), min.segment.length = 0, box.padding = 0.3, seed=1) + # , ylim=c(1.4,2.5)
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=16),
        legend.title=element_text(size=12), legend.text=element_text(size=11),
        legend.position = c(0.2, 0.75),
        legend.background = element_rect(fill = "white", color = "black"),
        aspect.ratio=1)



#### Figure 2b - JLIM vs Coloc scatterplot ####

# Reading JLIM and Coloc results
jlim_v_coloc <- read.csv('../data/Fig2/JLIM_COLOC.all.stat.txt', header = T, sep ='\t')

## FDR threshold for JLIM statistics
# (0.01172 * 1621) / sum(jlim_v_coloc$jlimp < 0.01172) ~ 0.05
# Therefore, JLIM p value threshold 0.01172 will be used for FDR 5%
jlim_threshold <- 0.01172

### scatter plot of significance ###
p_2b <- ggplot(jlim_v_coloc) + theme_classic() +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_vline(xintercept = -log10(jlim_threshold), linetype = 2) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.15, ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0.01) +
  geom_point(aes(x=-log10(jlim_p), y = coloc_H4, color = status), alpha=1) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.15, ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0, color='black', size=1.1) +
  scale_color_manual("Significance", values = c('red', 'orange', 'skyblue', 'gray'), labels=c("Both", "Coloc only", "JLIM only", "Neither")) +
  xlab(expression(paste(-log[10],"(JLIM ",italic(p),")"))) + ylab("Coloc PP(Colocalization)") +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text = element_text(size=14)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .04))) + scale_x_continuous(expand = expansion(mult = c(.04, .04))) +
  theme(aspect.ratio=1)


#### Figure 2b - Colocalized loci PU.1 motif altering variants bar plot ####

# Reading PU.1 motif alteration summary
motif_altered <- read.csv('../data/Fig2/colocalized_motif_altered_summary.txt', header = T, sep ='\t')

# Bar plot summarizing PU.1 motif alteration in colocalized loci
p_2c <- motif_altered %>% mutate(type = factor(type, levels=c("SNP", "Indel", "CNV", "Multi", "Unk"))) %>%
  mutate(motif = factor(motif, levels=c("+", "-",  "?"))) %>%
  ggplot(aes(x = type, y= value, fill=motif)) +
  theme_classic() +
  geom_col() + geom_text(aes(label = value), position = position_stack(vjust = 0.5), color="black", size=5) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.position = c(0.7, 0.8), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  xlab("Variant type") + ylab(expression(paste("Number of PU.1 bQTLs"))) +
  scale_fill_manual(values=c('#F8766D', '#00A5FF','gray'), name="Binding change",
                    labels=c("Gained", "Lost", "Unknown")) +
  scale_y_continuous(expand = expansion(mult = c(0, .05)))



#### Figure 2c - PU.1 motif-altering variant location ####

# Reading variant position along the PU.1 motif
motif_pos <- read.csv('../data/Fig2/jlim_coloc.SNP.PU1_motif_positions.table.txt', header = T, sep ='\t')

# Bar plot showing the number of variants at each motif position in colocalized loci
p_2d <- ggplot(motif_pos, aes(x=Motif_position, fill=Direction)) +
  geom_histogram(binwidth = 1, col='black') +
  theme_minimal_hgrid() +
  theme(aspect.ratio = 0.3, axis.title.x=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y= element_text(size=18), axis.text.y = element_text(size=14),
        legend.title=element_blank(), legend.text = element_text(size=14), legend.position = 'top', legend.justification = "center") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) + scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  scale_fill_manual(values=c('#F8766D', '#00A5FF'))


#### Figure 2d - Credible set size (GWAS) ####

# Reading credible set size table
cred_size <- read.csv('../data/Fig2/jlim_coloc.cred_size.table.txt', header = T, sep ='\t')

p_2e <- cred_size %>% mutate(size = factor(size, levels=c("1", "2-5", "6-10", "11-20", "21-50", ">50"))) %>%
ggplot(aes(x=size, y=count)) +
  geom_col(col='black', fill="black") +
  theme_classic() + background_grid(major = 'y') + #theme_minimal_hgrid() +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab("Credible set size") + ylab("# Colocalized PU.1 binding sites")



#### Putting all figure panels together and saving image

p_2_top <- p_2a + p_2b + plot_layout(ncol=2, widths = c(1,1))

p_2_bottom <- plot_grid(p_2c, p_2d, p_2e, ncol = 3, rel_widths = c(1,1.5,1.16))

p_2_bottom <- p_2c + p_2d + p_2e + plot_layout(ncol=3, widths = c(1,1.5,1.16))

p_2_all <- plot_grid(p_2_top, p_2_bottom, nrow = 2, rel_heights = c(1,1))


ggsave('../figures/Fig2.pdf',
       plot = p_2_all,
       device='pdf',
       width=380, height=285, units="mm")
