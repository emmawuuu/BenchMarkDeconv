library(SimBu)
library(dplyr)
library(tidyr)


# Load in single cell data
sc_data <- readRDS("EMTABsce_healthy.rds")
`%notin%` <- function(x, y) !(x %in% y)
cell_types_to_remove <- c("unclassified", "unclassified endocrine", "co-expression")
sc_data_filtered <- sc_data[, colData(sc_data)$cellType %notin% cell_types_to_remove]

# Normalize the single-cell data using TPM normalization
sc_data_tpm <- counts(sc_data_filtered) %>% 
  t() %>% 
  as.data.frame() %>% 
  apply(2, function(x) exp(x - log(sum(exp(x))))) %>% 
  t() %>% 
  as.matrix()

# Estimate cell type proportions
cell_type <- colData(sc_data_filtered)$cellType
cell_type_prop <- prop.table(table(colData(sc_data_filtered)$cellType))

# Combine cells of the same cell type using dplyr
sc_data_tbl <- as.data.frame(t(log1p(counts(sc_data_filtered))))
sc_data_tbl$cell_type <- colData(sc_data_filtered)$cellType
bulk_data <- sc_data_tbl %>%
  group_by(cell_type) %>%
  summarize(across(everything(), ~ mean(exp(.)), .names = "mean_{.col}")) %>%
  select(-cell_type) %>%
  t() %>%
  as.matrix()

annotation <- data.frame(ID = colnames(counts(sc_data_filtered)),
                         cell_type = colData(sc_data_filtered)$cellType)

# Create the dataset
ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts(sc_data_filtered),
  tpm_matrix = sc_data_tpm,
  filter_genes = TRUE
)

# Simulate pseudo bulk datasets
simulation <- SimBu::simulate_bulk(
  data = ds,
  scenario = "random",
  scaling_factor = "NONE",
  ncells = 100,
  nsamples = 10,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4),
  run_parallel = TRUE,
)

plot_object <- SimBu::plot_simulation(simulation = simulation)

# Reformat and write true prop to csv
spread_data <- spread(data = plot_object$data, key = cell_type, value = fraction, fill = 0)
write.csv(spread_data, "true_prop.csv", row.names = FALSE)

# Write simulated bulk counts data to csv 
bulk_counts <- SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]]
bulk_counts_matrix <- as.matrix(bulk_counts)
write.csv(bulk_counts_matrix, "bulk_counts.csv", row.names = FALSE)

# Write simultated bulk tpm normalized data to csv
bulk_tpm <- SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]]
bulk_tpm_matrix <- as.matrix(bulk_tpm)
write.csv(bulk_tpm_matrix, "bulk_tpm.csv", row.names = FALSE)
