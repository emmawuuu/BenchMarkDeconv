# These are deconvolution methods that use scRNA-seq as reference 

# Read in single cell and bulk data
sc_data <- readRDS("EMTABsce_healthy.rds")
bulk_counts <- read.csv("bulk_counts.csv")
bulk_counts <- t(as.matrix(bulk_counts))

# Run MuSiC
library("MuSiC")
library("Biobase")
RESULTS.sim.MuSic <- MuSiC::music_prop(
  bulk.mtx = bulk_counts,
  sc.sce = sc_data,
  clusters = 'cellType',
  samples = 'sampleID'
)

# deConvSeq

# DWLS

# Bisque
require(BisqueRNA)
library(Biobase)
library(SingleCellExperiment)
bulk.eset <- ExpressionSet(assayData = bulk_counts)
count_matrix <- assays(sc_data)$counts
rownames(count_matrix) <- rownames(sc_data)
colnames(count_matrix) <- colnames(sc_data)
metadata <- as.data.frame(colData(sc_data))
sc.eset <- ExpressionSet(assayData = count_matrix, phenoData = AnnotatedDataFrame(metadata))
RESULTS.BisqueRNA <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)$bulk.props 
