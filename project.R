if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("granulator")

library(granulator)
load_ABIS()

# get cell type reference matrices and perform deconvolution
sigList = list(
  ABIS_S0 = sigMatrix_ABIS_S0,
  ABIS_S1 = sigMatrix_ABIS_S1, 
  ABIS_S2 = sigMatrix_ABIS_S2, 
  ABIS_S3 = sigMatrix_ABIS_S3)
plot_similarity(sigMatrix=sigList)

decon <- deconvolute(m = bulkRNAseq_ABIS, sigMatrix = sigList)
# print decon results
decon$proportions$svr_ABIS_S0[1:5, 1:5]
plot_proportions(deconvoluted = decon, method = 'svr', signature = 'ABIS_S0')

# plot results for each model (reference matrix )
plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)


# deconvolution benchmarking 
bench <- benchmark(deconvoluted = decon, ground_truth = groundTruth_ABIS)
head(bench$summary)
head(bench$rank)

plot_benchmark(benchmarked = bench, metric = 'pcc')


# correlation of deconvolution analysis 
# use s2 signature matrix
# deconvolute input data using selected methods and reference profile matrix
methods <- c('ols','nnls','qprog','rls','svr', 'dtangle')
decon <- deconvolute(bulkRNAseq_ABIS, list(ABIS_S2 = sigMatrix_ABIS_S2), methods)

# correlation analysis
correl <- correlate(deconvoluted = decon)

# correlation heatmap
plot_correlate(correlated = correl, method="heatmap", legend=TRUE)
