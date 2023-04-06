# PU1-colocalization-manuscript
This repository contains code and data for generating figures to Jeong and Bulyk (in prep).

Updated April 6th 2023. Raehoon Jeong (rjeong@g.harvard.edu), Martha Bulyk Lab (mlbulyk@genetics.med.harvard.edu).


### Note
- Genome track images of PU.1 ChIP-seq, ATAC-seq, etc. were generated using <a href="https://software.broadinstitute.org/software/igv/">IGV browser</a>.
- Fuji plot in Figure 3 was generated using code from <a href="https://github.com/mkanai/fujiplot">Masahiro Kanai</a>.
- Canonical transcripts of genes were plotted using database from <a href="https://github.com/mkanai/locusviz/tree/master/inst/extdata">Locusviz</a>.  
- The <a href="https://github.com/BulykLab/PU1-colocalization-manuscript/blob/main/figures/r2_panel.pdf">legend for colors reflecting LD with selected variants</a> were added individually.
- For the other R scripts, the following packages need to be installed.
```
ggplot2
dplyr
cowplot
forcats
patchwork
data.table
AnnotationDbi
Ggbio
GenomicRanges
locuscomparer
ComplexUpset
```

### Workflow
1) Clone the repository
```
git clone git@github.com:BulykLab/PU1-colocalization-manuscript.git
```
2) Move to codes directory
```
cd PU1-colocalization-manuscript/codes
```
3) Run each code to generate figures using R
```
Rscript Figure2.R   # or any other code in the directory
```
