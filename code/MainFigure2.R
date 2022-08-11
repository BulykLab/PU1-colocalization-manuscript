library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(forcats, quietly = TRUE)
library(patchwork, quietly = TRUE)


#### Figure 2a - JLIM vs Coloc scatterplot ####

# Reading JLIM and Coloc results
jlim_v_coloc <- read.csv('../data/Fig2/JLIM_COLOC.all.stat.txt', header = F, sep ='\t')

## FDR threshold for JLIM statistics
# (0.01172 * 1621) / sum(jlim_v_coloc$jlimp < 0.01172) ~ 0.05
# Therefore, JLIM p value threshold 0.01172 will be used for FDR 5%
jlim_threshold <- 0.01172

# Determining statistical significance
jlim_v_coloc["jlim"] <- ifelse(jlim_v_coloc$V11 <=  jlim_threshold, 1, 0)
jlim_v_coloc["coloc"] <- ifelse(jlim_v_coloc$V17 >= 0.5, 1, 0)
jlim_v_coloc["status"] <- ifelse(jlim_v_coloc$jlim & jlim_v_coloc$coloc, "both", ifelse(jlim_v_coloc$jlim & !jlim_v_coloc$coloc,"jlim_only", ifelse(jlim_v_coloc$coloc, "coloc_only", "neither")))
jlim_v_coloc["jlimp"] <- ifelse(jlim_v_coloc$V11 == 0, 10^(-5.1), jlim_v_coloc$V11)


### scatter plot of significance ###
p_2a <- ggplot(jlim_v_coloc) + theme_classic() +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_vline(xintercept = -log10(jlim_threshold), linetype = 2) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.28, ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0.01) +
  geom_point(aes(x=-log10(jlimp), y = V17, color = status), alpha=1) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.28, ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0, color='black', size=1.1) +
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
p_2b <- motif_altered %>% mutate(type = factor(type, levels=c("SNP", "Indel", "CNV", "Multi", "Unk"))) %>%
  ggplot(aes(x = type, y= value, fill=motif)) + theme_classic() +
  geom_col() + geom_text(aes(label = value), position = position_stack(vjust = 0.5), color="black", size=5) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.position = c(0.7, 0.8), legend.title = element_text(size=16), legend.text = element_text(size=14)) +
  xlab("Variant type") + ylab(expression(paste("Number of PU.1 bQTLs"))) +
  scale_fill_manual(values=c('#00A5FF','gray','#F8766D'), name="Binding change",
                    labels=c("Lost", "Unknown", "Gained")) +
  scale_y_continuous(expand = expansion(mult = c(0, .05)))



#### Figure 2c - PU.1 motif-altering variant location ####

# Reading variant position along the PU.1 motif
motif_pos <- read.csv('../data/Fig2/jlim_coloc.SNP.PU1_motif_positions.table.txt', header = T, sep ='\t')

# Bar plot showing the number of variants at each motif position in colocalized loci
p_2c <- ggplot(motif_pos, aes(x=position, fill=change)) +
  geom_histogram(binwidth = 1, col='black') +
  theme_minimal_hgrid() +
  theme(aspect.ratio = 0.3, axis.title.x=element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y= element_text(size=18), axis.text.y = element_text(size=14),
        legend.title=element_blank(), legend.text = element_text(size=14), legend.position = 'top', legend.justification = "center") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) + scale_x_continuous(expand = expansion(mult = c(0.05, 0))) +
  scale_fill_manual(values=c('#F8766D', '#00A5FF'))


#### Figure 2d - Credible set size (GWAS) ####

# Reading credible set size table
cred_size <- read.csv('../data/Fig2/jlim_coloc.cred_size.table.txt', header = T, sep ='\t')

p_2d <- cred_size %>% mutate(size = factor(size, levels=c("1", "2-5", "6-10", "11-20", "21-50", ">50"))) %>%
ggplot(aes(x=size, y=count)) +
  geom_col(col='black', fill="black") +
  theme_minimal_hgrid() +
  theme(aspect.ratio = 1,
        axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab("Credible set size") + ylab("# Colocalized PU.1 binding sites")



#### Putting all figure panels together and saving image

p_2_top <- p_2a + p_2b + plot_layout(ncol=2, widths = c(4.7,3.3))
p_2_bottom <- plot_grid(p_2c, NULL, p_2d, ncol = 3, rel_widths = c(4,1.1,2.9))
p_2_all <- plot_grid(p_2_top, p_2_bottom, nrow = 2, rel_heights = c(7,6))


ggsave('../figures/Fig2.pdf',
       plot = p_2_all,
       device='pdf',
       width=320, height=270, units="mm")
