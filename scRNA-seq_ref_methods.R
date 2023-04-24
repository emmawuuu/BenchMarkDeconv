# These are deconvolution methods that use scRNA-seq as reference 

# Read in single cell and bulk data
library(SingleCellExperiment)
sc_data <- readRDS("EMTABsce_healthy.rds")
`%notin%` <- function(x, y) !(x %in% y)
cell_types_to_remove <- c("unclassified", "unclassified endocrine", "co-expression")
sc_data <- sc_data[, colData(sc_data)$cellType %notin% cell_types_to_remove]
bulk_counts <- readRDS("bulk_counts.rds")
true_prop <- readRDS("true_prop.rds")
true_prop_matrix <- as.matrix(true_prop)

# Run MuSiC
require(MuSiC)
library("Biobase")
RESULTS.sim.MuSic <- MuSiC::music_prop(
  bulk.mtx = as.matrix(bulk_counts),
  sc.sce = sc_data,
  clusters = 'cellType',
  samples = 'sampleID',
  normalize = FALSE,
  markers = NULL)$Est.prop.weighted

#select.ct = c("delta", "alpha","gamma","ductal","acinar",
#              "beta","MHC class II","PSC", "endothelial","epsilon","mast"),

# Bisque
require(BisqueRNA)
library(Biobase)
library(SingleCellExperiment)

bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk_counts))
sc_pDdata <- as.data.frame(colData(sc_data))
sc.eset <- Biobase::ExpressionSet(assayData = as.matrix(assay(sc_data)), phenoData = Biobase::AnnotatedDataFrame(sc_pDdata))

RESULTS.BisqueRNA <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)$bulk.props 
BisqueRNA.correlations <- apply(RESULTS.BisqueRNA, 2, function(x) cor(x, as.numeric(true_prop_matrix[,1])))

# deConvSeq
sc_celltypes = as.character(sc_pDdata$cellType) 
singlecell = model.matrix(~ -1 + as.factor(sc_celltypes))
colnames(singlecell) = levels(as.factor(sc_celltypes))
rownames(singlecell) = colnames(sc.eset)

dge = deconvSeq::getdge(sc.eset, singlecell, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
b0 = deconvSeq::getb0.rnaseq(dge, singlecell, ncpm.min =1, nsamp.min = 4)
dge.bulk = deconvSeq::getdge(bulk.eset, NULL, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")

RESULTS.deconvSeq = t(deconvSeq::getx1.rnaseq(NB0 = "top_fdr",b0, dge.bulk)$x1)

# DWLS

# SCDC
require(SCDC)
RESULTS.SCDC <- SCDC::SCDC_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, ct.varname = "cellType", sample = "SubjectName", ct.sub = unique(as.character(sc_pDdata$cellType)), iter.max = 200)$prop.est.mvw




