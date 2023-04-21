library(SimBu)
library(dplyr)

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
  filter_genes = TRUE,
  variance_cutoff = 0.1,
  type_abundance_cutoff = 10
)

custom_scenario_dataframe <- data.frame(
  "delta" = c(0.1124, 0.0562, 0.2810, 0.1686, 0.0562, 0.4496, 0.6744, 0.2248, 0.8992, 0.4496),
  "alpha" = c(0.1739, 0.2319, 0.0579, 0.1158, 0.1739, 0.1158, 0.0579, 0.1739, 0.0579, 0.1158),
  "gamma" = c(0.0870, 0.1739, 0.2609, 0.1739, 0.1739, 0.0870, 0.1739, 0.0870, 0.0870, 0.2609),
  "ductal" = c(0.0702, 0.0702, 0.1404, 0.1404, 0.0702, 0.1404, 0.0702, 0.2106, 0.1404, 0.0702),
  "acinar" = c(0.0413, 0.0826, 0.0413, 0.0826, 0.0413, 0.0413, 0.0413, 0.0826, 0.0413, 0.0826),
  "beta" = c(0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385),
  "PSC" = c(0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385),
  "endothelial" = c(0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385, 0.0385),
  row.names = paste0("sample", 1:10)
)

# Simulate pseudo bulk datasets
simulation <- SimBu::simulate_bulk(
  data = ds,
  scenario = "custom",
  scaling_factor = "NONE",
  ncells = 100,
  nsamples = 10,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4),
  run_parallel = TRUE,
  custom_scenario_data = custom_scenario_dataframe
)
