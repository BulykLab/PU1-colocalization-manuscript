library(DESeq2)
library(patchwork)

atac_counts <-  read.table("../../data/FigS7/rs411_ATAC.readcount.031722.txt", header=TRUE, row.names=1)

atac_counts <- atac_counts[,6:14]

coldata = data.frame(
  row.names= colnames(atac_counts),
  condition = c(rep("rs411_Edited", 3), rep("rs411_Het", 3), rep("rs411_WT", 3)))


dds <- DESeqDataSetFromMatrix(countData = atac_counts,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)


res_pu1ko <- results(dds,contrast=c("condition","rs411_Edited", "rs411_WT"))

write.table(as.data.frame(res_pu1ko), sep = '\t', file = '../../data/FigS7/res_pu1ko.txt', quote=F)
