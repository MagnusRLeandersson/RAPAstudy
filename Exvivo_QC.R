#______________________________ LIBRARIES ______________________________________

library(tibble)         # version 3.3.1
library(dplyr)          # version 1.2.1
library(tidyr)          # version 1.3.2
library(ggplot2)        # version 4.0.3
library(ComplexHeatmap) # version 2.26.1
library(circlize)       # version 0.4.18
library(RColorBrewer)   # version 1.1.3
library(Rtsne)          # version 0.17

#______________________________ DATA REQUIREMENTS ______________________________
# Ensure the following files are available in the working directory:
# - Table_S9_exvivo_sample_key.xlsx         -> Imported as 'Sample_key_eFT508'
# - Table_S10_exvivo_phos.xlsx              -> Imported as 'eFT508_phos'

#______________________________ FILTER DATA ____________________________________

# Format data frame and set rownames to GenePhos if present
eFT508_phos <- as.data.frame(eFT508_phos)
if ("GenePhos" %in% colnames(eFT508_phos)) {
  rownames(eFT508_phos) <- eFT508_phos$GenePhos
  eFT508_phos$GenePhos <- NULL
}

# 60% valid-value filtering across all samples
threshold_filt <- 0.6 * ncol(eFT508_phos)
Phospho_data_filtered <- eFT508_phos[rowSums(!is.na(eFT508_phos)) >= threshold_filt, ]

# Remove duplicate phosphosites
Phospho_data_filtered <- distinct(Phospho_data_filtered)


#______________________________ SAMPLE CORRELATION MATRIX ______________________

# Pairwise Pearson correlation and complete-linkage hierarchical clustering
cor_matrix <- cor(Phospho_data_filtered, use = "pairwise.complete.obs", method = "pearson")
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "complete")

# Color ramp: 0.85 to 1.0 (white to red)
col_fun_cor <- colorRamp2(c(0.85, 1), c("white", "red"))

# Generate Correlation Heatmap
Heatmap(
  cor_matrix,
  name = "Pearson\nCorr",
  col = col_fun_cor,
  cluster_rows = hc,
  cluster_columns = hc,
  column_dend_height = unit(20, "mm"),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  column_title = "Sample Correlation",
  border = TRUE
)


#______________________________ PCA PLOTS ______________________________________

# Clean complete cases and transpose
Phospho_data_clean <- Phospho_data_filtered[complete.cases(Phospho_data_filtered), ]
data_t <- t(Phospho_data_clean)

# Compute PCA with scaling
pca_result <- prcomp(data_t, scale. = TRUE)

# Prepare dataframe and merge metadata
pca_df <- as.data.frame(pca_result$x)
pca_df$Subject_Insulin_Drug <- rownames(pca_df)
pca_df <- merge(pca_df, Sample_key_eFT508, by = "Subject_Insulin_Drug")

# Calculate variance explained per principal component
percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# Base theme configuration for PCA plots
pca_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )

# Plot PC1 vs PC2 (with legend)
pca_p1_p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subject)) +
  geom_point(size = 3) +
  geom_text(aes(label = Insulin_Drug), vjust = -1.2, size = 2.5, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Subject"
  ) +
  pca_theme +
  theme(legend.position = "right")

print(pca_p1_p2)

# Plot PC3 vs PC4 (legend removed)
pca_p3_p4 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = Subject)) +
  geom_point(size = 3) +
  geom_text(aes(label = Insulin_Drug), vjust = -1.2, size = 2.5, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = paste0("PC3 (", percent_var[3], "%)"),
    y = paste0("PC4 (", percent_var[4], "%)")
  ) +
  pca_theme +
  theme(legend.position = "none")

print(pca_p3_p4)


#______________________________ t-SNE PLOT _____________________________________

# Set seed for reproducibility and prepare complete-case input
set.seed(42)
tsne_input <- t(na.omit(Phospho_data_filtered))

# Compute t-SNE (Perplexity 7)
tsne_out_combined <- Rtsne(tsne_input, perplexity = 7)

# Extract coordinates and merge with metadata
tsne_plot_combined <- data.frame(
  x = tsne_out_combined$Y[, 1],
  y = tsne_out_combined$Y[, 2],
  Subject_Insulin_Drug = rownames(tsne_input)
)
tsne_plot_combined <- merge(tsne_plot_combined, Sample_key_eFT508, by = "Subject_Insulin_Drug")

# Plot t-SNE
ggplot(tsne_plot_combined, aes(x = x, y = y, color = Subject)) +
  geom_point(size = 3) +
  geom_text(aes(label = Insulin_Drug), vjust = -1.2, size = 2.5, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "t-SNE dimension 1",
    y = "t-SNE dimension 2",
    color = "Subject"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )


#______________________________ DELTA eFT508 RESHAPING _________________________

# Transpose and merge with Sample_key_eFT508
df_Phos_subset <- as.data.frame(t(Phospho_data_filtered)) %>%
  rownames_to_column(var = "Subject_Insulin_Drug")

Phos <- merge(Sample_key_eFT508, df_Phos_subset, by = "Subject_Insulin_Drug")

# Reshape to long format
melted_Phos <- Phos %>%
  pivot_longer(
    cols = -(1:which(colnames(Phos) == "Sample_ID")),
    names_to = "Site",
    values_to = "Phos"
  )

# Calculate paired delta (eFT508 - DMSO)
delta_eFT508 <- melted_Phos %>%
  group_by(Site, Subject_Insulin) %>%
  summarize(
    DeltaPhos = Phos[Drug == "eFT508"] - Phos[Drug == "DMSO"],
    .groups = "drop"
  )

# Join delta back to long-format dataset
melted_Phos <- left_join(melted_Phos, delta_eFT508, by = c("Site", "Subject_Insulin"))


#______________________________ ALL CONDITIONS Z-SCORE HEATMAP _________________

# Pivot to matrix
melted_eFT508 <- melted_Phos %>%
  dplyr::select(Subject_Insulin_Drug, Site, Phos) %>%
  distinct() %>%
  pivot_wider(names_from = Subject_Insulin_Drug, values_from = Phos) %>%
  column_to_rownames("Site")

# Retain sites with >= 8 valid values
filtered_eFT508 <- melted_eFT508[rowSums(!is.na(melted_eFT508)) >= 8, ]

# Matrix conversion and row-wise Z-score scaling
eFT508_matrix <- as.matrix(filtered_eFT508)
eFT508_scaled <- t(scale(t(eFT508_matrix)))
eFT508_scaled <- eFT508_scaled[complete.cases(eFT508_scaled), ]

# Clustering
col_dist_all <- dist(t(eFT508_scaled), method = "euclidean")
col_clustering_all <- hclust(col_dist_all, method = "complete")

row_cor_all <- cor(t(eFT508_scaled), use = "pairwise.complete.obs")
row_dist_all <- as.dist(1 - row_cor_all)
row_dist_all[is.na(row_dist_all)] <- max(row_dist_all, na.rm = TRUE)
row_clustering_all <- hclust(row_dist_all, method = "ward.D2")

# Color ramp (-3 to 3)
col_fun_z <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

Heatmap(
  eFT508_scaled,
  name = "Z-score",
  col = col_fun_z,
  heatmap_legend_param = list(at = c(-3, 0, 3), labels = c("-3", "0", "3")),
  na_col = "grey80",
  cluster_rows = row_clustering_all,
  cluster_columns = col_clustering_all,
  row_split = NULL,
  column_split = 2,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_title = NULL,
  row_title = paste(nrow(eFT508_scaled), "sites"),
  column_names_rot = 90,
  border = TRUE,
  use_raster = FALSE
)


#______________________________ DELTA eFT508 Z-SCORE HEATMAP ___________________

# Pivot delta values to matrix
melted_eFT508_delta <- melted_Phos %>%
  dplyr::select(Subject_Insulin, Site, DeltaPhos) %>%
  distinct() %>%
  pivot_wider(names_from = Subject_Insulin, values_from = DeltaPhos) %>%
  column_to_rownames("Site")

# Retain sites with >= 4 valid delta values
filtered_delta <- melted_eFT508_delta[rowSums(!is.na(melted_eFT508_delta)) >= 4, ]

# Matrix conversion and row-wise Z-score scaling
delta_matrix <- as.matrix(filtered_delta)
delta_scaled <- t(scale(t(delta_matrix)))
delta_scaled <- delta_scaled[complete.cases(delta_scaled), ]

# Clustering
col_dist_delta <- dist(t(delta_scaled), method = "euclidean")
col_clustering_delta <- hclust(col_dist_delta, method = "complete")

row_cor_delta <- cor(t(delta_scaled), use = "pairwise.complete.obs")
row_dist_delta <- as.dist(1 - row_cor_delta)
row_dist_delta[is.na(row_dist_delta)] <- max(row_dist_delta, na.rm = TRUE)
row_clustering_delta <- hclust(row_dist_delta, method = "ward.D2")

Heatmap(
  delta_scaled,
  name = "Z-score",
  col = col_fun_z,
  heatmap_legend_param = list(at = c(-3, 0, 3), labels = c("-3", "0", "3")),
  na_col = "grey80",
  cluster_rows = row_clustering_delta,
  cluster_columns = col_clustering_delta,
  row_split = NULL,
  column_split = 2,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_title = expression(paste(Delta, " eFT508")),
  row_title = paste(nrow(delta_scaled), "sites"),
  column_names_rot = 90,
  border = TRUE,
  use_raster = FALSE
)




