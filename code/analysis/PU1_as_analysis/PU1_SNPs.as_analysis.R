library(dplyr)

### Performing allelic imbalance analysis

pu1_as_data <- read.table('PU1_SNPs.data_matrix.txt', header = T)
pu1_as_snps <- read.table('PU1_SNPs.rsid.txt', header = T)


# Add pseudocount of 0.5
pu1_as_data$lafc <- pu1_as_data$het_value * log2((pu1_as_data$alt_count + 0.5) / (pu1_as_data$ref_count + 0.5))

pu1_as_data$wt <- ((pu1_as_data$alt_count + 0.5) * (pu1_as_data$ref_count + 0.5)) / ((pu1_as_data$alt_count + 0.5) + (pu1_as_data$ref_count + 0.5))


pu1_as_results <- data.frame(rsid = character(),
                             beta_as = double(),
                             se_as = double(),
                             p_as = double(),
                             stringsAsFactors = FALSE)

for (snp in pu1_as_snps$rsid){
  result <- summary(lm(formula=lafc ~ 0+ het_value, data=pu1_as_data %>% filter(rsid == snp), weights = wt))
  pu1_as_results[nrow(pu1_as_results) + 1, ] <- c(snp, as.double(result$coefficients[1]), as.double(result$coefficients[2]), as.double(result$coefficients[4]))
}


## Saving to file
filename <- 'PU1_SNPs.as_results.txt'
write.table(pu1_as_results, file = filename, quote = FALSE, sep = '\t', row.names=F)