# PU.1 motif gkm-SVM model
This repository contains data for gkm-SVM motif score prediction.

Updated August 22nd 2022.

### Note
- <a href=https://github.com/Dongwon-Lee/lsgkm>LS-GKM</a> must be installed to replicate gkm-SVM score predictions.

### Code
For instance, this code makes predictions for 30bp sequences centered at SNPs within PU.1 binding motifs.
```
lsgkm/bin/gkmpredict PU1.snps.30bp.fa PU1.gkmSVM.model.txt PU1.snps.30bp.predictions.txt
```
