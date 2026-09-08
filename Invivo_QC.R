#______________________________ LIBRARIES ______________________________________

library(tibble)          # version 3.3.1
library(qvalue)          # version 2.42.0
library(dplyr)           # version 1.2.1
library(ggplot2)         # version 4.0.3
library(tidyr)           # version 1.3.2
library(ComplexHeatmap)  # version 2.26.1
library(circlize)        # version 0.4.18
library(clusterProfiler) # version 4.18.4
library(magick)          # version 2.9.1
library(eulerr)          # version 8.0.0
library(grid)            # version 4.5.2 (base)
library(RColorBrewer)    # version 1.1.3
library(Rtsne)           # version 0.17


#______________________________ DATA REQUIREMENTS ______________________________
# Ensure the following files are available in the working directory:
# - Table_S3_invivo_glucoseuptake.xlsx      -> Imported as 'Sample_key_RAPA'
# - Table_S4_invivo_phos.xlsx               -> Imported as 'RAPA_phos'

#______________________________ FILTER DATA ____________________________________
# Filter for 70% detected values
RAPA_phos <- as.data.frame(RAPA_phos)
rownames(RAPA_phos) <- RAPA_phos$GenePhos
RAPA_phos$GenePhos <- NULL

threshold_filt <- 0.7 * ncol(RAPA_phos)
Phospho_data_filtered <- RAPA_phos[rowSums(!is.na(RAPA_phos)) >= threshold_filt, ]

Phospho_data_filtered <- distinct(Phospho_data_filtered)


#______________________________ CORRELATION MATRIX _____________________________

# Calculate Correlation
cor_matrix <- cor(Phospho_data_filtered, use = "pairwise.complete.obs", method = "pearson")
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "complete")

# Create Custom Labels from Sample_key_RAPA and ensure order matches cor_matrix
mapping_df <- Sample_key_RAPA[match(colnames(cor_matrix), Sample_key_RAPA$SampleID), ]

new_labels <- paste0(
  mapping_df$Subject, "_", 
  substr(mapping_df$Drug, 1, 4), "_", 
  substr(mapping_df$Leg, 1, 2), "_", 
  substr(mapping_df$Clamp, 1, 3)
)

# Plot matrix
col_fun_cor = colorRamp2(c(0.75, 1), c("white", "red"))
Heatmap(cor_matrix, 
        name = "Pearson\nCorr", 
        col = col_fun_cor,
        cluster_rows = hc, 
        cluster_columns = hc,
        show_row_dend = FALSE,          
        show_column_dend = TRUE,        
        column_dend_height = unit(20, "mm"),
        row_labels = new_labels,
        row_names_side = "left",        
        row_names_gp = gpar(fontsize = 5),
        show_column_names = FALSE,      
        column_title = "Sample Correlation",
        border = TRUE)



#______________________________ t-SNE PLOT _____________________________

# Set seed for reproducibility and prepare transposed, complete data
set.seed(42)
tsne_input <- t(na.omit(Phospho_data_filtered))

# Calculate t-SNE (Perplexity 30)
tsne_out_combined <- Rtsne(tsne_input, perplexity = 30) 

# Extract coordinates and merge with metadata
tsne_plot_combined <- data.frame(
  x = tsne_out_combined$Y[,1], 
  y = tsne_out_combined$Y[,2], 
  SampleID = rownames(tsne_input)
)
tsne_plot_combined <- merge(tsne_plot_combined, Sample_key_RAPA, by = "SampleID")

# Factorize variables to ensure consistent alphanumeric sorting and shape mapping
tsne_plot_combined$Subject <- factor(tsne_plot_combined$Subject, levels = paste0("S", 1:13))
tsne_plot_combined$Drug <- factor(tsne_plot_combined$Drug)

# Define color palette for subjects
subject_colors <- c(RColorBrewer::brewer.pal(12, "Paired"), "#999999")

# Generate t-SNE plot mapping color to Subject and shape to Drug
ggplot(tsne_plot_combined, aes(x = x, y = y, color = Subject, shape = Drug)) + 
  geom_point(size = 4, alpha = 0.8) + 
  scale_color_manual(values = subject_colors) +
  scale_shape_manual(values = c("Placebo" = 16, "Rapamycin" = 17)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank()
  ) +
  guides(shape = guide_legend(override.aes = list(color = "black"))) +
  labs(
    x = "t-SNE dimension 1", 
    y = "t-SNE dimension 2",
    color = "Subject ID",
    shape = "Treatment"
  )





