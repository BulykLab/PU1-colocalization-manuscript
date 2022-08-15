library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(forcats, quietly = TRUE)
library(patchwork, quietly = TRUE)


#### Figure 2a - JLIM vs Coloc scatterplot ####

# Reading JLIM and Coloc results
jlim_v_coloc <- read.csv('../data/Fig2/JLIM_COLOC.all.stat.txt', header = T, sep ='\t')

## FDR threshold for JLIM statistics
# (0.01172 * 1621) / sum(jlim_v_coloc$jlimp < 0.01172) ~ 0.05
# Therefore, JLIM p value threshold 0.01172 will be used for FDR 5%
jlim_threshold <- 0.01172

# Determining statistical significance
jlim_v_coloc["jlim"] <- ifelse(jlim_v_coloc$jlim_p <=  jlim_threshold, 1, 0)
jlim_v_coloc["coloc"] <- ifelse(jlim_v_coloc$coloc_H4 >= 0.5, 1, 0)
jlim_v_coloc["status"] <- ifelse(jlim_v_coloc$jlim & jlim_v_coloc$coloc, "both", ifelse(jlim_v_coloc$jlim & !jlim_v_coloc$coloc,"jlim_only", ifelse(jlim_v_coloc$coloc, "coloc_only", "neither")))
jlim_v_coloc["jlimp"] <- ifelse(jlim_v_coloc$jlim_p == 0, 1/(10^5+1), jlim_v_coloc$jlim_p)

## Supp Figure 2
## JLIM X Coloc O
p_supp_2a <- ggplot(jlim_v_coloc) + theme_classic() +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_vline(xintercept = -log10(jlim_threshold), linetype = 2) +
  geom_rect(xmin=-0.1, xmax=-log10(jlim_threshold), ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0.01) +
  geom_point(aes(x=-log10(jlimp), y = coloc_H4, color = status), alpha=1) +
  geom_rect(xmin=-0.1, xmax=-log10(jlim_threshold), ymin=0.5, ymax=1.02, fill="#FFF4C0", alpha=0, color='black', size=1.1) +
  scale_color_manual("Significance", values = c('red', 'orange', 'skyblue', 'gray'), labels=c("Both", "Coloc only", "JLIM only", "Neither")) +
  xlab(expression(paste(-log[10],"(JLIM ",italic(p),")"))) + ylab("Coloc PP(Colocalization)") +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=17), legend.text = element_text(size=15)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .04))) + scale_x_continuous(expand = expansion(mult = c(.04, .04))) +
  theme(aspect.ratio=1)

## JLIM O Coloc X
p_supp_2b <- ggplot(jlim_v_coloc) + theme_classic() +
  geom_hline(yintercept = 0.5, linetype = 2) + geom_vline(xintercept = -log10(jlim_threshold), linetype = 2) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.15, ymin=-0.035, ymax=0.5, fill="#FFF4C0", alpha=0.01) +
  geom_point(aes(x=-log10(jlimp), y = coloc_H4, color = status), alpha=1) +
  geom_rect(xmin=-log10(jlim_threshold), xmax=5.15, ymin=-0.035, ymax=0.5, fill="#FFF4C0", alpha=0, color='black', size=1.1) +
  scale_color_manual("Significance", values = c('red', 'orange', 'skyblue', 'gray'), labels=c("Both", "Coloc only", "JLIM only", "Neither")) +
  xlab(expression(paste(-log[10],"(JLIM ",italic(p),")"))) + ylab("Coloc PP(Colocalization)") +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=17), legend.text = element_text(size=15)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .04))) + scale_x_continuous(expand = expansion(mult = c(.04, .04))) +
  theme(aspect.ratio=1)


ggsave('../figures/Supp2_a.pdf',
       plot = p_supp_2a,
       device='pdf',
       width=200, height=160, units="mm")

ggsave('../figures/Supp2_b.pdf',
       plot = p_supp_2b,
       device='pdf',
       width=200, height=160, units="mm")
