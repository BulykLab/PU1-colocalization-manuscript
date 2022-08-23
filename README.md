# PU1-colocalization-manuscript
This repository contains codes and data for generating figures to Jeong and Bulyk (in prep).

Updated August 23rd 2022. Raehoon Jeong (rjeong@g.harvard.edu), Martha Bulyk Lab (mlbulyk@genetics.med.harvard.edu).


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
Rscript MainFigure2.R   # or any other code in the directory
```

### Note
- Genome track images of PU.1 ChIP-seq, ATAC-seq, etc. were generated using <a href="https://software.broadinstitute.org/software/igv/">IGV browser</a>.
- Fuji plot in Figure 3 was generated using code from <a href="https://github.com/mkanai/fujiplot">Masahiro Kanai</a>.
- The <a href="https://github.com/BulykLab/PU1-colocalization-manuscript/blob/main/figures/r2_panel.pdf">legend for colors reflecting LD with selected variants</a> were added individually.
