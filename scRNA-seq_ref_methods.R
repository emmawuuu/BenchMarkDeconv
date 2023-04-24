# These are deconvolution methods that use scRNA-seq as reference 

# Read in single cell and bulk data
library(SingleCellExperiment)
sc_data <- readRDS("EMTABsce_healthy.rds")
`%notin%` <- function(x, y) !(x %in% y)
cell_types_to_remove <- c("unclassified", "unclassified endocrine", "co-expression") # Clean single cell data
# 11 cell types
sc_data <- sc_data[, colData(sc_data)$cellType %notin% cell_types_to_remove]
# 4 cell types
sub_data <- sc_data[, colData(sc_data)$cellType %in% c("alpha", "beta", "gamma", "delta")]
bulk_counts <- readRDS("bulk_counts.rds")

# Read in true bulk prop
true_prop <- readRDS("true_prop.rds")
true_prop_matrix <- as.matrix(true_prop)
true_prop_matrix.sub <- true_prop_matrix[, c("sample", "alpha", "beta", "gamma", "delta")] 


# extract the cell types from the true proportion matrix
cell_types <- colnames(true_prop_matrix)[-1]
cell_types.sub <- colnames(true_prop_matrix.sub)[-1]

#########
# MuSiC # 
#########
# ISSUE: Not all cell types show up?
require(MuSiC)
library("Biobase")

common_genes <- intersect(rownames(as.matrix(bulk_counts)), rownames(counts(sc_data)))

RESULTS.sim.MuSic <- MuSiC::music_prop(
  bulk.mtx = as.matrix(bulk_counts),
  sc.sce = sc_data,
  clusters = 'cellType',
  samples = 'sampleID',
  select.ct = c("delta", "alpha","gamma","ductal","acinar",
               "beta","MHC class II","PSC", "endothelial","epsilon","mast"),
  normalize = FALSE,
  markers = common_genes)$Est.prop.weighted


##########
# Bisque #
##########
require(BisqueRNA)
library(Biobase)
library(SingleCellExperiment)

# Create expression set of bulk rna seq data
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk_counts))

# Create expression set of single cell reference with 11 cell types
sc_pDdata <- as.data.frame(colData(sc_data))
sc.eset <- Biobase::ExpressionSet(assayData = as.matrix(assay(sc_data)), 
                                  phenoData = Biobase::AnnotatedDataFrame(sc_pDdata))
# Run 11 cell types
RESULTS.BisqueRNA <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, 
                                                            sc.eset, markers=NULL,
                                                            use.overlap=FALSE)$bulk.props 

# Create expression set of single cell reference with 4 cell types
sc_pDdata.sub <- as.data.frame(colData(sub_data))
sc.eset.sub <- Biobase::ExpressionSet(assayData = as.matrix(assay(sub_data)), 
                                  phenoData = Biobase::AnnotatedDataFrame(sc_pDdata.sub))

# Run 4 cell types
RESULTS.BisqueRNA.sub <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, 
                                                                sc.eset.sub, markers=NULL,
                                                            use.overlap=FALSE)$bulk.props

saveRDS(RESULTS.BisqueRNA, "BisqueRNA_results.rds")
saveRDS(RESULTS.BisqueRNA.sub, "BisqueRNA_sub_results.rds")

#############
# deConvSeq #
#############

# ISSUE - comes up as all nans?
sc_celltypes = as.character(sc_pDdata$cellType) 
design = model.matrix(~ -1 + as.factor(sc_celltypes))
colnames(design) = levels(as.factor(sc_celltypes))
colnames(design) <- gsub(" ", "_", colnames(design))
rownames(design) = colnames(sc.eset)

dge = deconvSeq::getdge(counts(sc_data), design, ncpm.min = 1, nsamp.min = 4)
set.seed(1234)
b0 = deconvSeq::getb0.rnaseq(dge, design, ncpm.min =1, nsamp.min = 4, sigg = NULL)

dge.bulk = deconvSeq::getdge(bulk.eset, NULL, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")

RESULTS.deconvSeq = deconvSeq::getx1.rnaseq(NB0 = "top_fdr",b0, dge.bulk)$x1

########
# DWLS #
########

# Issue with dimension for beta cell types? Idk
require(DWLS)

workdir <- getwd()
dir.create(file.path(workdir, "data"), showWarnings = FALSE) #folder to store data
dir.create(file.path(workdir, "results"), showWarnings = FALSE) #folder to save results

Signature <- DWLS::buildSignatureMatrixMAST(scdata = counts(sc_data), id = as.character(sc_pDdata$cellType), 
                                      path = "results", diff.cutoff = 0.5, pval.cutoff = 0.01)

########
# SCDC #
########
require(SCDC)
# Run SCDC for 11 cell types
RESULTS.SCDC <- SCDC::SCDC_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, 
                                ct.varname = "cellType", sample = "SubjectName", 
                                ct.sub = unique(as.character(sc_pDdata$cellType)), 
                                iter.max = 200)$prop.est.mvw
# Run SCDC for 4 cell types
RESULTS.SCDC.sub <- SCDC::SCDC_prop(bulk.eset = bulk.eset, sc.eset = sc.eset.sub, 
                                ct.varname = "cellType", sample = "SubjectName", 
                                ct.sub = unique(as.character(sc_pDdata.sub$cellType)), 
                                iter.max = 200)$prop.est.mvw

saveRDS(RESULTS.SCDC, "SCDC_results.rds")
saveRDS(RESULTS.SCDC.sub, "SCDC_results_sub.rds")





