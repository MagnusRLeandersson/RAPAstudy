#______________________________ LIBRARIES ______________________________________

library(limma)     # version 3.66.0
library(tibble)    # version 3.3.1
library(qvalue)    # version 2.42.0
library(dplyr)     # version 1.2.1
library(tidyr)     # version 1.3.2
library(ggplot2)   # version 4.0.3
library(purrr)     # version 1.2.2
library(patchwork) # version 1.3.2
library(openxlsx)  # version 4.2.8.1

#______________________________ DATA REQUIREMENTS ______________________________
# Ensure the following files are available in the working directory:
# - Table_S9_exvivo_sample_key.xlsx         -> Imported as 'Sample_key_eFT508'
# - Table_S10_exvivo_phos.xlsx              -> Imported as 'eFT508_phos'
# - Table_S3_invivo_glucoseuptake.xlsx      -> Imported as 'Sample_key_RAPA'
# - Table_S4_invivo_phos.xlsx               -> Imported as 'RAPA_phos'
# - Table_S6_invivo_rapamycin_sites.xlsx    -> Imported as 'RAPA_sites' (Sheet: "Rapa_sign")


#______________________________ FILTER DATA _________________________

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


#______________________________ DATA PREPARATION & OUTLIER REMOVAL _____________

# Reshape phospho data to long format and merge with Sample_key
df_Phos_subset <- as.data.frame(t(Phospho_data_filtered)) %>%
  rownames_to_column(var = "Subject_Insulin_Drug")

Phos <- merge(Sample_key_eFT508, df_Phos_subset, by = "Subject_Insulin_Drug")

melted_Phos <- Phos %>%
  pivot_longer(
    cols = -(1:which(colnames(Phos) == "Sample_ID")),
    names_to = "Site",
    values_to = "Phos"
  )

# Exclude outlier sample S4_Basal (identified during QC)
melted_Phos <- melted_Phos %>%
  filter(!Subject_Insulin %in% c("S4_Basal"))

Sample_key_eFT508 <- Sample_key_eFT508 %>%
  filter(!Subject_Insulin %in% c("S4_Basal"))


#______________________________ MODULAR LIMMA PIPELINE FUNCTION ________________

run_limma_eft508 <- function(data_long, meta_key, min_delta_obs = 6) {
  
  # Calculate paired eFT508 delta (eFT508 - DMSO) per Subject_Insulin pair
  delta_df <- data_long %>%
    group_by(Site, Subject_Insulin) %>%
    summarize(
      DeltaPhos = Phos[Drug == "eFT508"] - Phos[Drug == "DMSO"],
      .groups = "drop"
    )
  
  # Join delta back to long-format data
  data_long_joined <- left_join(data_long, delta_df, by = c("Site", "Subject_Insulin"))
  
  # Filter sites by minimum non-NA delta observations
  data_long_filtered <- data_long_joined %>%
    group_by(Site) %>%
    filter(sum(!is.na(DeltaPhos)) >= min_delta_obs) %>%
    ungroup()
  
  # Compute median logFC across paired deltas
  median_deltas <- delta_df %>%
    group_by(Site) %>%
    summarize(median_logFC = median(DeltaPhos, na.rm = TRUE), .groups = "drop")
  
  # Reconstruct wide expression matrix
  Phos_matrix <- data_long_filtered %>%
    dplyr::select(Site, Subject_Insulin_Drug, Phos) %>%
    distinct() %>%
    pivot_wider(names_from = Subject_Insulin_Drug, values_from = Phos) %>%
    column_to_rownames("Site") %>%
    as.matrix()
  
  Phos_matrix <- Phos_matrix[!duplicated(Phos_matrix), ]
  
  # Align metadata with matrix columns
  Design_sub <- as.data.frame(meta_key) %>%
    filter(Subject_Insulin_Drug %in% colnames(Phos_matrix)) %>%
    arrange(match(Subject_Insulin_Drug, colnames(Phos_matrix)))
  
  # Linear model design and contrast matrix
  Design_matrix <- model.matrix(~ 0 + Drug, data = Design_sub)
  Contrasts <- makeContrasts(DrugeFT508 - DrugDMSO, levels = Design_matrix)
  
  # Run limma with Subject_Insulin blocking for paired design
  corfit <- duplicateCorrelation(Phos_matrix, Design_matrix, block = Design_sub$Subject_Insulin)
  fit <- lmFit(Phos_matrix, Design_matrix, block = Design_sub$Subject_Insulin, correlation = corfit$consensus)
  fit <- contrasts.fit(fit, Contrasts)
  fit <- eBayes(fit)
  
  # Extract statistics and append median logFC and Storey q-values
  Results <- topTable(fit, number = Inf, genelist = rownames(Phos_matrix))
  colnames(Results)[1] <- "Site"
  Results <- left_join(Results, median_deltas, by = "Site")
  Results$qvalue <- qvalue(p = Results$P.Value)$qvalues
  
  return(list(
    results = Results,
    melted_filtered = data_long_filtered
  ))
}


#______________________________ RUNNING LIMMA FOR CONDITIONS ___________________

# Combined (Basal + Insulin)
out_comb <- run_limma_eft508(
  data_long = melted_Phos,
  meta_key = Sample_key_eFT508,
  min_delta_obs = 10
)
Results_eFT508_comb <- out_comb$results
melted_Phos_eFT508  <- out_comb$melted_filtered

# Basal Only
out_basal <- run_limma_eft508(
  data_long = melted_Phos %>% filter(Insulin != "Insulin"),
  meta_key = Sample_key_eFT508 %>% filter(Insulin != "Insulin"),
  min_delta_obs = 6
)
Results_eFT508_basal <- out_basal$results

# Insulin Only
out_insulin <- run_limma_eft508(
  data_long = melted_Phos %>% filter(Insulin != "Basal"),
  meta_key = Sample_key_eFT508 %>% filter(Insulin != "Basal"),
  min_delta_obs = 6
)
Results_eFT508_insulin <- out_insulin$results


#______________________________ MERGE RESULTS & DEFINE SIGNIFICANCE ____________

# Combine all 3 limma condition results into a single table
df1_eFT508 <- Results_eFT508_comb    %>% rename_with(~paste0(., "_comb"), -Site)
df2_eFT508 <- Results_eFT508_basal   %>% rename_with(~paste0(., "_basal"), -Site)
df3_eFT508 <- Results_eFT508_insulin %>% rename_with(~paste0(., "_insulin"), -Site)

Results_eFT508_all <- list(df1_eFT508, df2_eFT508, df3_eFT508) %>%
  purrr::reduce(full_join, by = "Site")

# Merge phosphosite annotation metadata if available
if (exists("Phos_info")) {
  Results_eFT508_all <- left_join(Results_eFT508_all, Phos_info, by = c("Site" = "GenePhos"))
}

# Filter significant eFT508-regulated sites (P < 0.05 and |Median logFC| > log2(1.5))
Results_eFT508_sign <- Results_eFT508_all %>%
  filter(
    (P.Value_comb < 0.05    & abs(median_logFC_comb)    > log2(1.5)) |
      (P.Value_basal < 0.05   & abs(median_logFC_basal)   > log2(1.5)) |
      (P.Value_insulin < 0.05 & abs(median_logFC_insulin) > log2(1.5))
  ) %>%
  separate(Site, into = c("Gene", "Residue", "M"), sep = "_", remove = FALSE)

cat("Significant Genes:", paste(unique(Results_eFT508_sign$Gene), collapse = ", "), "\n")


#______________________________ TRANSLATIONAL OVERLAP WITH IN VIVO RAPA ________

# Standardize RAPA_sites column identifiers using Median_logFC
rapa_ref <- RAPA_sites
if (!"ID" %in% colnames(rapa_ref) & "Site" %in% colnames(rapa_ref)) {
  rapa_ref <- rename(rapa_ref, ID = Site)
}
if ("Median_logFC" %in% colnames(rapa_ref)) {
  rapa_ref <- rename(rapa_ref, Median_logFC_Rapa = Median_logFC)
}

# Overlap using combined eFT508 results (opposite directionality: eFT508 * Rapa < 0)
eFT508_Rapa_overlap <- inner_join(
  Results_eFT508_sign %>% dplyr::select(Site, logFC_eFT508 = median_logFC_comb),
  rapa_ref %>% dplyr::select(ID, Median_logFC_Rapa),
  by = c("Site" = "ID")
) %>%
  filter((logFC_eFT508 * Median_logFC_Rapa) < 0)

# Comprehensive overlap across any condition with opposite directionality
eFT508_overlap <- Results_eFT508_all %>%
  inner_join(rapa_ref, by = c("Site" = "ID"))

opposite_direction <- eFT508_overlap %>%
  filter(
    (P.Value_comb < 0.05    & abs(median_logFC_comb)    > log2(1.5) & (median_logFC_comb * Median_logFC_Rapa) < 0) |
      (P.Value_basal < 0.05   & abs(median_logFC_basal)   > log2(1.5) & (median_logFC_basal * Median_logFC_Rapa) < 0) |
      (P.Value_insulin < 0.05 & abs(median_logFC_insulin) > log2(1.5) & (median_logFC_insulin * Median_logFC_Rapa) < 0)
  )


#______________________________ PREPARE HUMAN IN VIVO RAPA TARGET DATA _________

# Extract target sites from overlap
overlap_sites <- unique(opposite_direction$Site)

# Filter raw RAPA dataset directly for overlapping sites without 70% cutoff
RAPA_phos_df <- as.data.frame(RAPA_phos)
if ("GenePhos" %in% colnames(RAPA_phos_df)) {
  rownames(RAPA_phos_df) <- RAPA_phos_df$GenePhos
  RAPA_phos_df$GenePhos <- NULL
}

RAPA_filtered_targets <- RAPA_phos_df[rownames(RAPA_phos_df) %in% overlap_sites, , drop = FALSE]

Phos_long_base_RAPA <- RAPA_filtered_targets %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Site",
    values_to = "Phos"
  ) %>%
  mutate(Phos = as.numeric(Phos)) %>%
  inner_join(Sample_key_RAPA, by = "SampleID")


#______________________________ VISUALIZATION OF OVERLAPPING TARGET SITES _______

# Generate unified comparison for each overlapping site
for (curr_site in overlap_sites) {
  
  # --- Ex Vivo eFT508 Plot ---
  plot_data_eft508 <- melted_Phos_eFT508 %>%
    filter(Site == curr_site) %>%
    mutate(
      Condition = factor(Insulin, levels = c("Basal", "Insulin")),
      Drug = factor(Drug, levels = c("DMSO", "eFT508"))
    )
  
  p_eft508 <- ggplot(plot_data_eft508, aes(x = Drug, y = Phos)) +
    stat_summary(
      fun = mean, geom = "errorbar",
      aes(ymax = after_stat(y), ymin = after_stat(y)),
      color = "black", width = 0.6, linewidth = 1.2
    ) +
    geom_line(aes(group = Subject_Insulin), color = "grey60", linewidth = 0.4) +
    geom_point(aes(color = Drug), size = 3.5, alpha = 0.9) +
    geom_text(aes(label = Subject), size = 2, color = "black", fontface = "bold") +
    facet_wrap(~Condition, strip.position = "bottom") +
    scale_color_manual(values = c("DMSO" = "skyblue3", "eFT508" = "firebrick2")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black"),
      strip.placement = "outside",
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = paste0("Ex Vivo eFT508: ", curr_site),
      x = NULL,
      y = "Phos Intensity"
    )
  
  # --- In Vivo Rapamycin Raw Data Plot ---
  site_data_rapa <- Phos_long_base_RAPA %>% filter(Site == curr_site)
  
  plot_data_orig_rapa <- site_data_rapa %>%
    mutate(
      Condition_Combined = factor(
        paste(Clamp, Leg),
        levels = c("Basal Rest", "Basal Exercise", "Insulin Rest", "Insulin Exercise")
      ),
      Drug = factor(Drug, levels = c("Placebo", "Rapamycin"))
    )
  
  p_rapa_raw <- ggplot(plot_data_orig_rapa, aes(x = Drug, y = Phos)) +
    stat_summary(
      fun = mean, geom = "errorbar",
      aes(ymax = after_stat(y), ymin = after_stat(y)),
      color = "black", width = 0.6, linewidth = 1.2
    ) +
    geom_line(aes(group = Subject), color = "grey60", linewidth = 0.4, alpha = 0.6) +
    geom_point(aes(color = Drug), size = 3.5, alpha = 0.9) +
    geom_text(aes(label = Subject), size = 2, color = "black", fontface = "bold") +
    facet_wrap(~Condition_Combined, nrow = 1, strip.position = "bottom") +
    scale_color_manual(values = c("Placebo" = "lightblue3", "Rapamycin" = "darkorange2")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black"),
      strip.placement = "outside",
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = paste0("In Vivo RAPA (Raw): ", curr_site),
      x = NULL,
      y = "Phos Intensity"
    )
  
  # --- In Vivo Rapamycin Delta Plot (Exercise - Rest) ---
  plot_data_delta_rapa <- site_data_rapa %>%
    dplyr::select(Subject, Drug, Clamp, Leg, Phos) %>%
    pivot_wider(names_from = Leg, values_from = Phos) %>%
    mutate(
      Delta_Phos = Exercise - Rest,
      Condition_Combined = factor(
        if_else(Clamp == "Basal", "∆PEX Basal", "∆PEX Insulin"),
        levels = c("∆PEX Basal", "∆PEX Insulin")
      ),
      Drug = factor(Drug, levels = c("Placebo", "Rapamycin"))
    )
  
  p_rapa_delta <- ggplot(plot_data_delta_rapa, aes(x = Drug, y = Delta_Phos)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
    stat_summary(
      fun = mean, geom = "errorbar",
      aes(ymax = after_stat(y), ymin = after_stat(y)),
      color = "black", width = 0.6, linewidth = 1.2
    ) +
    geom_line(aes(group = Subject), color = "grey60", linewidth = 0.4, alpha = 0.6) +
    geom_point(aes(color = Drug), size = 3.5, alpha = 0.9) +
    geom_text(aes(label = Subject), size = 2, color = "black", fontface = "bold") +
    facet_wrap(~Condition_Combined, nrow = 1, strip.position = "bottom") +
    scale_color_manual(values = c("Placebo" = "lightblue3", "Rapamycin" = "darkorange2")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black"),
      strip.placement = "outside",
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = paste0("In Vivo RAPA (Delta): ", curr_site),
      x = NULL,
      y = "∆ Phos Intensity"
    )
  
  # Combined Side-by-Side Panel: Ex Vivo (left) | RAPA Raw (middle) | RAPA Delta (right)
  combined_panel <- p_eft508 + p_rapa_raw + p_rapa_delta + plot_layout(widths = c(1, 1.8, 1))
  
  print(combined_panel)
}

### Visually evaluate the regulation of these sites and exclude based on the following criteria: 
      ### 1. Opposite regulation by drug in different conditions conditions making it unreliable (e.g. Basal and Insulin)
      ### 2. Large variation between samples and in effect sizes, i.e. probable false positive
      ### 3. Too large a fraction of samples regulated in opposite direction from main effect, i.e., 
      ### if 6/13 samples goes slightly down, but main effect goes up because of larger effect sizes.

#______________________________ EXPORT DATA TO EXCEL ___________________________

p_thresh   <- 0.05
lfc_thresh <- log2(1.5)

# Split significant and non-significant sites for volcano dataset
Results_eFT508_sign_export <- subset(
  Results_eFT508_comb,
  P.Value < p_thresh & abs(median_logFC) > lfc_thresh
)

Results_eFT508_nonsign_export <- subset(
  Results_eFT508_comb,
  !(P.Value < p_thresh & abs(median_logFC) > lfc_thresh)
)

# Export all result tables and opposite direction overlap to Excel
write.xlsx(
  list(
    "Significant_eFT508"      = Results_eFT508_sign_export,
    "Non_Significant_eFT508"  = Results_eFT508_nonsign_export,
    "Opposite_Direction_Rapa" = opposite_direction
  ),
  file = "~/Desktop/Table_S11_exvivo_eFT508_sites.xlsx",
  rowNames = FALSE
)




