# BenchMarkDeconv

## Requirements

```
install.packages("BiocManager")
install.packages("DWLS")
BiocManager::install("Biobase")
BiocManager::install("TOAST")
BiocManager::install("SingleCellExperiment")
devtools::install_github("xuranw/MuSiC") 
devtools::install_github("rosedu1/deconvSeq")
devtools::install_github("cozygene/bisque")
```
## Dataset 
[Single cell RNA-seq data of pancreatic islets from healthy individuals (Segerstolpe et al.)](https://xuranw.github.io/MuSiC/data/EMTABsce_healthy.rds)

#### Bulk Simulated Data
- 10 samples 
- 11 cell types ("delta", "alpha", "gamma", "ductal", "acinar", "beta", "MHC class II", "PSC", "endothelial", "epsilon", "mast")
- 25453 genes

1. Simulated Bulk Data: bulk_counts.csv
2. Simulated Bulk Data TPM Normalized: bulk_tpm.csv
3. True Proportions of Bulk Simulated Data: true_prop.csv
