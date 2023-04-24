library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

# Read in true bulk prop
true_prop <- readRDS("true_prop.rds")
true_prop_matrix <- as.matrix(true_prop)
true_prop_matrix.sub <- true_prop_matrix[, c("sample", "alpha", "beta", "gamma", "delta")] 

# extract the cell types from the true proportion matrix
cell_types <- colnames(true_prop_matrix)[-1]
cell_types.sub <- colnames(true_prop_matrix.sub)[-1]

# Read in bulk deconvolution method results
bulk_deconv_results = readRDS("final_result.rds")

# Read in scRNA-reference deconvolution method results
RESULTS.BisqueRNA = readRDS("BisqueRNA_results.rds")
RESULTS.BisqueRNA.sub = readRDS("BisqueRNA_sub_results.rds")
RESULTS.SCDC = readRDS("SCDC_results.rds")
RESULTS.SCDC.sub = readRDS("SCDC_results_sub.rds")

# Calculate the Pearson correlation for each 11 cell types
BisqueRNA.correlations <- sapply(cell_types, 
                                 function(ct) cor(as.numeric(t(as.matrix(RESULTS.BisqueRNA))[, ct]), 
                                                  as.numeric(true_prop_matrix[, ct])))

# Calculate the RMSE for each 11 cell types
BisqueRNA.RMSEs <- sapply(cell_types, function(ct) {
  estimated <- as.numeric(t(as.matrix(RESULTS.BisqueRNA))[, ct])
  true <- as.numeric(true_prop_matrix[, ct])
  sqrt(mean((true - estimated)^2))
})

# Calculate the Pearson correlation for each 4 cell types
BisqueRNA.correlations.sub <- sapply(cell_types.sub, 
                                     function(ct) cor(as.numeric(t(as.matrix(RESULTS.BisqueRNA.sub))[, ct]), 
                                                      as.numeric(true_prop_matrix.sub[, ct])))

# Calculate the RMSE for each 4 cell types
BisqueRNA.RMSEs.sub <- sapply(cell_types.sub, function(ct) {
  estimated <- as.numeric(t(as.matrix(RESULTS.BisqueRNA.sub))[, ct])
  true <- as.numeric(true_prop_matrix.sub[, ct])
  sqrt(mean((true - estimated)^2))
})

# Calculate the Pearson correlation for each 11 cell types
SCDC.correlations <- sapply(cell_types, function(ct) cor(as.numeric(as.matrix(RESULTS.SCDC)[, ct]), 
                                                         as.numeric(true_prop_matrix[, ct])))
# Calculate the RMSE for each 11 cell types
SCDC.RMSEs <- sapply(cell_types, function(ct) {
  estimated <- as.numeric((as.matrix(RESULTS.SCDC))[, ct])
  true <- as.numeric(true_prop_matrix[, ct])
  sqrt(mean((true - estimated)^2))
})

# Calculate the Pearson correlation for each 4 cell types
SCDC.correlations.sub <- sapply(cell_types.sub, function(ct) cor(as.numeric(as.matrix(RESULTS.SCDC.sub)[, ct]), 
                                                                 as.numeric(true_prop_matrix.sub[, ct])))
# Calculate the RMSE for each 4 cell types
SCDC.RMSEs.sub <- sapply(cell_types.sub, function(ct) {
  estimated <- as.numeric((as.matrix(RESULTS.SCDC.sub))[, ct])
  true <- as.numeric(true_prop_matrix.sub[, ct])
  sqrt(mean((true - estimated)^2))
})

# Plot pearson correlation for scRNA-seq refernece methods

##########################
# All 11 cell types RMSE #
##########################
df_scRef_rmse <- data.frame(
  celltype = c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma",
                "mast", "MHC class II", "PSC"),
  SCDC = as.numeric(SCDC.RMSEs),
  Bisque = as.numeric(BisqueRNA.RMSEs)
)

# Reshape the data into a long format
df_scRef_long_rmse <- reshape2::melt(df_scRef_rmse, id.vars = "celltype", variable.name = "method", value.name = "rmse")

# Create a box plot of RMSE values with different colors for each method
ggplot(df_scRef_long_rmse, aes(x = celltype, y = rmse, color = method)) +
  geom_boxplot(aes(group = celltype), outlier.shape = NA) +
  geom_jitter(aes(shape = method), position = position_jitter(width = 0.3), size = 3, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "RMSE by Cell Type and Method",
       x = "Cell Type",
       y = "RMSE",
       color = "Method",
       shape = "Method")

#############################
# All 11 cell types Pearson #
#############################

df_scRef_pcc <- data.frame(
  celltype = c("acinar","alpha","beta","delta","ductal","endothelial","epsilon","gamma",
               "mast", "MHC class II", "PSC"),
  SCDC = as.numeric(SCDC.correlations),
  Bisque = as.numeric(BisqueRNA.correlations)
)

# Reshape the data into a long format
df_scRef_long_pcc <- reshape2::melt(df_scRef_rmse, id.vars = "celltype", variable.name = "method", value.name = "pcc")

# Calculate mean for each cell type and method
summary_df <- df_scRef_long_pcc %>%
  group_by(celltype, method) %>%
  summarise(mean_pcc = mean(pcc))

# Create a new dataset with pairwise comparisons between cell types
pairwise_comparison_df <- expand_grid(CellType1 = unique(summary_df$celltype),
                                      CellType2 = unique(summary_df$celltype),
                                      method = unique(summary_df$method))

# Calculate the absolute difference between PCC values
pairwise_comparison_df <- pairwise_comparison_df %>%
  left_join(summary_df, by = c("CellType1" = "celltype", "method")) %>%
  rename(PCC1 = mean_pcc) %>%
  left_join(summary_df, by = c("CellType2" = "celltype", "method")) %>%
  rename(PCC2 = mean_pcc) %>%
  mutate(Difference = abs(PCC1 - PCC2))

# Create the heatmap
plot <- ggplot(pairwise_comparison_df, aes(x = CellType1, y = CellType2, fill = Difference)) +
  geom_tile() +
  facet_wrap(~ Method, ncol = 1) +
  scale_fill_gradient2(low = "white", mid = "blue", high = "darkred", midpoint = median(pairwise_comparison_df$Difference)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PCC Difference between Cell Types",
       x = "Cell Type 1",
       y = "Cell Type 2",
       fill = "PCC Difference")

################
# 4 cell types #
################
df_scRef_rmse.sub <- data.frame(
  celltype = c("alpha","beta","gamma","delta"),
  SCDC = as.numeric(SCDC.RMSEs.sub),
  Bisque = as.numeric(BisqueRNA.RMSEs.sub)
)

# Reshape the data into a long format
df_scRef_long_rmse.sub <- reshape2::melt(df_scRef_rmse.sub, id.vars = "celltype", variable.name = "method", value.name = "rmse")

# Create a box plot of RMSE values with different colors for each method
ggplot(df_scRef_long_rmse.sub, aes(x = celltype, y = rmse, color = method)) +
  geom_boxplot(aes(group = celltype), outlier.shape = NA) +
  geom_jitter(aes(shape = method), position = position_jitter(width = 0.3), size = 3, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "RMSE by Cell Type and Method 4 Cell Types",
       x = "Cell Type",
       y = "RMSE",
       color = "Method",
       shape = "Method")


########################################
# Combine Methods RMSE - 11 Cell Types #
########################################

# Add OLS
df_rmse <- df_scRef_long_rmse %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$ols_sig1 %>% select(celltype, method, rmse))

# Add dtangle
df_rmse <- df_rmse %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$dtangle_sig1 %>% select(celltype, method, rmse))

# Add nnls
df_rmse <- df_rmse %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$nnls_sig1 %>% select(celltype, method, rmse))

# Add qprog
df_rmse <- df_rmse %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$qprog_sig1 %>% select(celltype, method, rmse))

# Add rls
df_rmse <- df_rmse %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$rls_sig1 %>% select(celltype, method, rmse))

# Create a box plot of RMSE values with different colors for each method
ggplot(df_rmse, aes(x = celltype, y = rmse, color = method)) +
  geom_boxplot(aes(group = celltype), outlier.shape = NA) +
  geom_jitter(aes(shape = method), position = position_jitter(width = 0.3), size = 3, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "RMSE by Cell Type and Method",
       x = "Cell Type",
       y = "RMSE",
       color = "Method",
       shape = "Method")


#######################################
# Combine Methods RMSE - 4 Cell Types #
#######################################

# Add OLS
df_rmse.sub <- df_scRef_long_rmse.sub %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$ols_sig2 %>% select(celltype, method, rmse))

# Add dtangle
df_rmse.sub <- df_rmse.sub %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$dtangle_sig2 %>% select(celltype, method, rmse))

# Add nnls
df_rmse.sub <- df_rmse.sub %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$nnls_sig2 %>% select(celltype, method, rmse))

# Add qprog
df_rmse.sub <- df_rmse.sub %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$qprog_sig2 %>% select(celltype, method, rmse))

# Add rls
df_rmse.sub <- df_rmse.sub %>%
  select(celltype, method, rmse) %>%
  rbind(final_results$stats$rls_sig2 %>% select(celltype, method, rmse))

# Create a box plot of RMSE values with different colors for each method
ggplot(df_rmse.sub, aes(x = celltype, y = rmse, color = method)) +
  geom_boxplot(aes(group = celltype), outlier.shape = NA) +
  geom_jitter(aes(shape = method), position = position_jitter(width = 0.3), size = 3, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "RMSE by Cell Type and Method - 4 Cell Types",
       x = "Cell Type",
       y = "RMSE",
       color = "Method",
       shape = "Method")

