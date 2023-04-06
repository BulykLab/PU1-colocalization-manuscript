library(ggplot2)


### Linear regression for PU1_67321
pu1_data <- read.table('PU1_67321.rs74267027.data.txt', header = T)

pu1_p <- cor.test(pu1_data$genotype, pu1_data$phenotype)$p.value
# 1.196314e-08
pu1_beta <- cor(pu1_data$genotype, pu1_data$phenotype) * sd(pu1_data$phenotype) / sd(pu1_data$genotype)
# -0.8381416


### Linear regression for ATAC_137925
atac_data <- read.table('ATAC_137925.rs74267027.data.txt', header = T)

atac_p <- cor.test(atac_data$genotype, atac_data$phenotype)$p.value
# 1.004786e-27
atac_beta <- cor(atac_data$genotype, atac_data$phenotype) * sd(atac_data$phenotype) / sd(atac_data$genotype)
# -0.9974773



### Linear regression for H3K27ac_17_16170296_16172772
h3k27ac_data <- read.table('H3K27ac_17_16170296_16172772.rs74267027.data.txt', header = T)

h3k27ac_p <- cor.test(h3k27ac_data$genotype, h3k27ac_data$phenotype)$p.value
# 7.644327e-47
h3k27ac_beta <- cor(h3k27ac_data$genotype, h3k27ac_data$phenotype) * sd(h3k27ac_data$phenotype) / sd(h3k27ac_data$genotype)
# -0.7646476



### Linear regression for H3K4me1_17_16170245_16172722
h3k4me1_data <- read.table('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/PU1_bloodcelltrait/CG_submission/reresubmission/analysis/QTLs/PU1_67321/H3K4me1_17_16170245_16172722.rs74267027.data.txt', header = T)

h3k4me1_p <- cor.test(h3k4me1_data$genotype, h3k4me1_data$phenotype)$p.value
# 1.241711e-09
h3k4me1_beta <- cor(h3k4me1_data$genotype, h3k4me1_data$phenotype) * sd(h3k4me1_data$phenotype) / sd(h3k4me1_data$genotype)
# -0.2986415


rs74267027_df <- data.frame(
  pu1_p = pu1_p,
  pu1_beta = pu1_beta,
  atac_p = atac_p,
  atac_beta = atac_beta,
  h3k27ac_p = h3k27ac_p,
  h3k27ac_beta = h3k27ac_beta,
  h3k4me1_p = h3k4me1_p,
  h3k4me1_beta = h3k4me1_beta
)

write.table(rs74267027_df, 'rs74267027.QTL_results.txt', sep='\t', quote=F, row.names=F)
