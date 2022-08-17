library(ggplot2)
library(patchwork)

### Need to run Rscript PU1-colocalization-manuscript/code/analysis/PU1KO_ATAC_analysis.R to generate DESeq2 results file

res_pu1ko <- read.csv('../data/SuppFig6/res_pu1ko.txt', header = T, sep = '\t')


rs411_ATAC_PU1 <- read.table("../data/SuppFig6/rs411_ATAC.PU1.bed", header=FALSE)
rs411_ATAC_noPU1 <- read.table("../data/SuppFig6/rs411_ATAC.noPU1.bed", header=FALSE)


## Plotting

## accessibility at PU.1 binding sites
p_atac_pu1 <- ggplot(data=res_pu1ko[rownames(res_pu1ko) %in% rs411_ATAC_PU1$V4,]) +
  geom_hline(yintercept = 0, alpha=0.5, ) +
  geom_point(aes(x=baseMean, y=log2FoldChange, color=(padj < 0.05)), size = 0.5) +
  scale_x_continuous(trans='log10', limits = c(3,2000) ) +
  scale_y_continuous(limits = c(-5,5) ) +
  scale_color_manual(values = c("gray", "red")) +
  theme_classic() +
  theme(legend.position = 'none', axis.title = element_text(size=14),
        axis.text = element_text(size=14), plot.title = element_text(size=14)) +
  labs(x="Mean accessibility", y="log2 fold change", title="With PU.1 binding")

## accessibility at sites without PU.1 binding
p_atac_nopu1 <- ggplot(data=res_pu1ko[rownames(res_pu1ko) %in% rs411_ATAC_noPU1$V4,]) +
  geom_hline(yintercept = 0, alpha=0.5, ) +
  geom_point(aes(x=baseMean, y=log2FoldChange, color=(padj < 0.05)), size = 0.5) +
  scale_x_continuous(trans='log10', limits = c(3,2000) ) +
  scale_y_continuous(limits = c(-5,5) ) +
  scale_color_manual(values = c("gray", "red")) +
  theme_classic() +
  theme(legend.position = 'none', axis.title = element_text(size=14),
        axis.text = element_text(size=14), plot.title = element_text(size=14)) +
  labs(x="Mean accessibility", y="log2 fold change", title="No PU.1 binding")



p_supp_6 <- p_atac_pu1 +
  p_atac_nopu1 + plot_layout(ncol = 2, widths = c(1, 1))

ggplot2::ggsave('../figures/Supp6.pdf',
                plot = p_supp_6,
                device='pdf',
                width=250, height=125, units="mm")
