# BenchMarkDeconv

## Dataset 
[Single cell RNA-seq data of pancreatic islets from healthy individuals (Segerstolpe et al.)](https://xuranw.github.io/MuSiC/data/EMTABsce_healthy.rds)

#### Bulk Simulated Data
- 10 samples 
- 11 cell types ("delta", "alpha", "gamma", "ductal", "acinar", "beta", "MHC class II", "PSC", "endothelial", "epsilon", "mast")
- 25453 genes

1. Simulated Bulk Data: bulk_counts.csv
2. Simulated Bulk Data TPM Normalized: bulk_tpm.csv
3. True Proportions of Bulk Simulated Data: true_prop.csv

#### final result 
- final result is in final_result.rds
- the rmse value of each cell type is in final_results$stats 
- the rmse values are given for 1) each deconvolution method 2) signature - this is the reference matrix used (sig1 - all 11 cell types // sig2 - 4 cell types) 
- final_results$stats$deconvolutionMethod_signature
- ex) final_results$stats$ols_sig1 should give you a matrix of the rmse values of each cell type deconvolution results when we use sig1(4cell types) and use ols (ordinary least squares) for deconvolution. The rmse column shows the rmse values, but there are also pcc (pearson coefficient) values availalbe. I would prefer to use rmse if possible
