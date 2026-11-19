# TCGA Glioblastoma survival analysis

### overview

Analysis of gene expression and clinical survival data from glioblastoma patients using R.

The analysis includes Kaplan–Meier survival curves (median and quartile splits), gene expression comparisons across molecular subtypes, and survival analysis according to temozolomide (TMZ) treatment and radiotherapy status.

### data sources

Publicly available data from The Cancer Genome Atlas (TCGA), accessed through UCSC Xena.

**gene expression dataset:**
- dataset: Gene expression array – AffyU133a
- description: TCGA glioblastoma multiforme (GBM) gene expression data by AffyU133a array
- hub: https://tcga.xenahubs.net
- dataset page: https://xenabrowser.net/datapages/?dataset=TCGA.GBM.sampleMap%2FHT_HG-U133A&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443


### analysis

- Kaplan–Meier survival analysis (median and quartile splits)
- log-rank tests
- Cox proportional hazards regression
- gene expression comparison across molecular subtypes
- Kruskal–Wallis and pairwise Wilcoxon tests
- survival analysis by TMZ and radiotherapy status

### tools

R, survival, survminer, ggplot2, ggpubr, tidyverse.
