library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(ComplexUpset)

### Figure 4a Number of overlap
pu1_qtl_numbers <- read.table('../data/Fig4/QTL_numbers_table.txt', header = T)

p_4a <- ggplot(pu1_qtl_numbers, aes(x=QTL, y=Num)) +
  geom_col(aes(fill = QTL, alpha=Type), width = 0.8, color = 'black') +
  theme_classic() +
  ylab("Number of \n colocalized PU.1 bQTLs") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
        legend.title=element_text(size=16), legend.text=element_text(size=14),
        aspect.ratio = 1.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values =c('#4051A3', '#007E5C', '#B89300')) +
  scale_alpha_manual(values = c(0, 0.3, 1), labels = c('No overlap', 'Just overlap', 'QTL')) +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title=element_blank(), legend.text=element_text(size=12))


### Figure 4b Upset plot
pu1_upset_input <- read.table('../data/Fig4/QTL_upset_input_wPU1.txt', header = T, sep='\t')

colnames(pu1_upset_input) = c("PU1_peak","PU.1 bQTL", "caQTL","H3K27ac hQTL", "H3K4me1 hQTL")
QTLs <- colnames(pu1_upset_input)[2:5]

p_4b <-upset(pu1_upset_input, QTLs, sort_intersections_by='degree',
              name = "",
              base_annotations=list('Intersection size'=intersection_size(counts=FALSE)),
              sort_sets=FALSE,
              queries=list(
                upset_query(set='PU.1 bQTL', fill='#46ACC8'),
                upset_query(set='caQTL', fill='#4051A3'),
                upset_query(set="H3K27ac hQTL", fill='#007E5C'),
                upset_query(set="H3K4me1 hQTL", fill='#B89300')
              ),
              themes=upset_modify_themes(
                list(
                  'Intersection size'=theme(axis.title.y=element_blank(), axis.text.y=element_text(size=12)),
                  'overall_sizes'=theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12)),
                  'intersections_matrix'=theme(axis.text.y=element_text(size=14))
                )
              )
              )


### Figure 4c QTL effect comparison
## PU1 and ATAC
pu1_atac_data <- read.table('../data/Fig4/PU1_ATAC_comparison.txt', header = T)

p_pu1_atac <- ggplot(pu1_atac_data, aes(x=beta_atac, y=beta_pu1,
                                        xmin=beta_atac - se_atac, xmax=beta_atac + se_atac,
                                        ymin=beta_pu1 - se_pu1, ymax=beta_pu1 + se_pu1,
                                        alpha=(significant == 'O'))) +
  geom_point(size = 2, color = '#4051A3') + geom_errorbar(color = '#4051A3') + geom_errorbarh(color = '#4051A3') +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_alpha_manual(values=c(0.3,1)) +
  lims(x = c(-2,2), y=c(-2,2)) +
  labs(x="caQTL effect", y="PU.1 bQTL effect") +
  theme_minimal() + 
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.position="none", aspect.ratio = 1)

#cor.test((pu1_atac_data %>% filter(significant == 'O'))$beta_pu1, (pu1_atac_data %>% filter(significant == 'O'))$beta_atac)$estimate
# 0.9392347
#cor.test((pu1_atac_data %>% filter(significant == 'O'))$beta_pu1, (pu1_atac_data %>% filter(significant == 'O'))$beta_atac)$p.value
# 4.023967e-14


## PU1 and H3K27ac
pu1_h3k27ac_data <- read.table('../data/Fig4/PU1_H3K27ac_comparison.txt', header = T)

p_pu1_h3k27ac <- ggplot(pu1_h3k27ac_data, aes(x=beta_histone, y=beta_pu1, 
                                              xmin=beta_histone - se_histone, xmax=beta_histone + se_histone,
                                              ymin=beta_pu1 - se_pu1, ymax=beta_pu1 + se_pu1,
                                              alpha=(significant == 'O'),
                                              color=(significant != 'Z'))) +
  geom_point(size = 2) + geom_errorbar() + geom_errorbarh() +
  guides(color = "none") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_alpha_manual(values =c(0.3,1), labels = c('Just overlap', 'QTL')) +
  scale_color_manual(values = c('#007E5C', 'white')) +
  lims(x = c(-2,2), y=c(-2,2)) +
  labs(x="H3K27ac hQTL effect", y="") +
  theme_minimal() + 
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.position="top", legend.title = element_blank(),
        legend.text=element_text(size=14), aspect.ratio = 1)

#cor.test((pu1_h3k27ac_data %>% filter(significant == 'O'))$beta_pu1, (pu1_h3k27ac_data %>% filter(significant == 'O'))$beta_histone)$estimate
# 0.9392347
#cor.test((pu1_h3k27ac_data %>% filter(significant == 'O'))$beta_pu1, (pu1_h3k27ac_data %>% filter(significant == 'O'))$beta_histone)$p.value
# 4.023967e-14


## PU1 and H3K4me1
pu1_h3k4me1_data <- read.table('../data/Fig4/PU1_H3K4me1_comparison.txt', header = T)

p_pu1_h3k4me1 <- ggplot(pu1_h3k4me1_data, aes(x=beta_histone, y=beta_pu1, 
                                              xmin=beta_histone - se_histone, xmax=beta_histone + se_histone,
                                              ymin=beta_pu1 - se_pu1, ymax=beta_pu1 + se_pu1,
                                              alpha=(significant == 'O'))) +
  geom_point(size = 2, color ='#B89300') + geom_errorbar(color ='#B89300') + geom_errorbarh(color ='#B89300') +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_alpha_manual(values = c(0.3,1)) +    ## Modify colors
  lims(x = c(-2,2), y=c(-2,2)) +
  labs(x="H3K4me1 hQTL effect") +
  theme_minimal() +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.position="none", aspect.ratio = 1)

#cor.test((pu1_h3k4me1_data %>% filter(significant == 'O'))$beta_pu1, (pu1_h3k4me1_data %>% filter(significant == 'O'))$beta_histone)$estimate
# 0.9145257 
#cor.test((pu1_h3k4me1_data %>% filter(significant == 'O'))$beta_pu1, (pu1_h3k4me1_data %>% filter(significant == 'O'))$beta_histone)$p.value
# 4.023967e-14


### Figure 4d Allelic imbalance plot
pu1_qtl_as <- read.table('../data/Fig4/PU1_QTL_AS_comparison.txt', header = T)

p_qtl_as <- ggplot(pu1_qtl_as, aes(x=beta_as, y=beta_pu1,
                                   xmin=beta_as - se_as, xmax=beta_as + se_as,
                                   ymin=beta_pu1 - se_pu1, ymax=beta_pu1 + se_pu1
                                   )) +
  geom_point(size = 2, color = '#46ACC8') + geom_errorbar(color = '#46ACC8') + geom_errorbarh(color = '#46ACC8') +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_alpha_manual(values =c(0.3,1)) +  
  lims(x=c(-4.5,4.5), y = c(-2,2)) +
  labs(x="PU.1 ChIP-seq \n allelic imbalance effect", y="PU.1 bQTL effect") +
  theme_minimal() + 
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.position="none", aspect.ratio = 1)

#cor.test(pu1_qtl_as_merged$beta_pu1, pu1_qtl_as_merged$beta_as)$estimate
# 0.9437686 
#cor.test(pu1_qtl_as_merged$beta_pu1, pu1_qtl_as_merged$beta_as)$p.value
# 4.156833e-43


### Figure 4e PU1 and eQTL
pu1_eqtl_data <- read.table('../data/Fig4/PU1_eqtl_comparison.txt', header = T)

p_pu1_eqtl <- ggplot(pu1_eqtl_data, aes(x=beta_eqtl, y=beta_pu1,
                                        xmin=beta_eqtl - se_eqtl, xmax=beta_eqtl + se_eqtl,
                                        ymin=beta_pu1 - se_pu1, ymax=beta_pu1 + se_pu1
                                        )) +
  geom_point(size = 2, color = 'black') + geom_errorbar(color = 'black') + geom_errorbarh(color = 'black') +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  lims(x = c(-2,2), y=c(-2,2)) +
  labs(x="eQTL effect", y="PU.1 bQTL effect") +
  theme_minimal() + 
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.position="none", aspect.ratio = 1)

#cor.test(pu1_eqtl_data$beta_pu1, pu1_eqtl_data$slope)$estimate
# 0.7971235 
#cor.test(pu1_eqtl_data$beta_pu1, pu1_eqtl_data$slope)$p.value
# 4.156833e-43


### Putting the panels together
p_4_top <- p_4a + p_4b + plot_layout(ncol = 2, widths = c(1,2))

p_4_middle <- p_pu1_atac + p_pu1_h3k27ac + p_pu1_h3k4me1 + plot_layout(ncol = 3)

p_4_bottom <- p_qtl_as + theme(plot.margin = unit(c(0,50,0,0), "pt")) +
  p_pu1_eqtl + plot_layout(ncol = 2)

p_4_all <- plot_grid(p_4_top, p_4_middle, p_4_bottom, nrow = 3, rel_heights = c(1,1,0.95))

ggplot2::ggsave('../figures/Fig4.pdf',
                plot = p_4_all,
                device='pdf',
                width=300, height=300, units="mm")