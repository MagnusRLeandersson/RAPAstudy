#______________________________ LIBRARIES ______________________________________

library(limma)           # version 3.66.0
library(tibble)          # version 3.3.1
library(qvalue)          # version 2.42.0
library(ReactomePA)      # version 1.54.0
library(AnnotationDbi)   # version 1.72.0
library(org.Hs.eg.db)    # version 3.22.0
library(dplyr)           # version 1.2.1
library(ComplexHeatmap)  # version 2.26.1
library(circlize)        # version 0.4.18
library(ggplot2)         # version 4.0.3
library(clusterProfiler) # version 4.18.4
library(rmcorr)          # version 0.7.0
library(doParallel)      # version 1.0.17
library(foreach)         # version 1.5.2
library(openxlsx)        # version 4.2.8.1

#______________________________ DATA REQUIREMENTS ______________________________
# Ensure the following files are available in the working directory:
# - Table_S3_invivo_glucoseuptake.xlsx      -> Imported as 'Sample_key_RAPA'
# - Table_S4_invivo_phos.xlsx               -> Imported as 'RAPA_phos'

#______________________________ FILTER DATA _____________________________

# Filter for 70% detected values
RAPA_phos <- as.data.frame(RAPA_phos)
rownames(RAPA_phos) <- RAPA_phos$GenePhos
RAPA_phos$GenePhos <- NULL

threshold_filt <- 0.7 * ncol(RAPA_phos)
Phospho_data_filtered <- RAPA_phos[rowSums(!is.na(RAPA_phos)) >= threshold_filt, ]

Phospho_data_filtered <- distinct(Phospho_data_filtered)


#______________________________ PERSONALIZED PHOSPHO STEP 1 - RAPAMYCIN EFFECT ____________________________________________

#___________ DATA PREPARATION FOR LIMMA ____________________

# Format long-format phosphoproteomic data merged with metadata
Phos_long_base <- Phospho_data_filtered %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  inner_join(Sample_key_RAPA, by = "SampleID") %>%
  reshape2::melt(
    id.vars = c("SampleID", "Subject", "Drug", "Leg", "Clamp", "GU_AUC", "GU_MEAN"),
    variable.name = "Site",
    value.name = "Phos"
  ) %>%
  mutate(Phos = as.numeric(Phos))


#_______________ LIMMA FUNCTION __________________

# Function to run paired delta calculations, limma modeling, and summary statistics
run_limma_pipeline <- function(data_long, meta_key, min_pairs = 8) {
  
  # Calculate paired Rapamycin delta (Rapamycin - Placebo) per Subject and Site
  Phos_delta <- data_long %>%
    group_by(Subject, Site) %>%
    mutate(Delta_Phos = Phos[Drug == "Rapamycin"] - Phos[Drug == "Placebo"]) %>%
    ungroup() %>%
    filter(Drug != "Placebo") %>%
    reshape2::dcast(Site ~ Subject, value.var = "Delta_Phos") %>%
    column_to_rownames("Site") %>%
    t()
  
  # Filter for phosphosites with at least N valid paired delta measurements
  valid_sites <- colnames(Phos_delta)[colSums(!is.na(Phos_delta)) >= min_pairs]
  Phos_delta_filtered <- Phos_delta[, valid_sites]
  
  # Construct matrix for limma modeling using only valid sites
  Phos_matrix <- data_long %>%
    filter(Site %in% valid_sites) %>%
    reshape2::dcast(Site ~ SampleID, value.var = "Phos") %>%
    column_to_rownames("Site")
  
  # Align design metadata with matrix columns
  Design <- meta_key %>%
    filter(SampleID %in% colnames(Phos_matrix)) %>%
    arrange(SampleID)
  
  # Linear model design and contrast matrix
  Design_matrix <- model.matrix(~ 0 + Drug, data = Design)
  Contrasts <- makeContrasts(DrugRapamycin - DrugPlacebo, levels = Design_matrix)
  
  # Run limma with subject blocking for repeated measures
  corfit <- duplicateCorrelation(Phos_matrix, Design_matrix, block = Design$Subject)
  fit <- lmFit(Phos_matrix, Design_matrix, block = Design$Subject, correlation = corfit$consensus)
  fit <- contrasts.fit(fit, Contrasts)
  fit <- eBayes(fit)
  
  # Extract statistics and calculate Storey q-values
  Results <- topTable(fit, number = Inf, genelist = rownames(Phos_matrix))
  Results$qvalue <- qvalue(p = Results$P.Value)$qvalues
  
  # Compute summary statistics across subjects for delta values
  Stats <- data.frame(
    ID = colnames(Phos_delta_filtered),
    Median_logFC = apply(Phos_delta_filtered, 2, median, na.rm = TRUE),
    Mean_logFC = apply(Phos_delta_filtered, 2, mean, na.rm = TRUE)
  )
  
  # Merge summary stats with limma output
  Results_all <- merge(Results, Stats, by = "ID", all.x = TRUE)
  
  # Filter significant sites (p < 0.05 and |Median_logFC| > 0.585)
  Results_sig <- Results_all %>%
    filter(P.Value < 0.05 & (Median_logFC < -0.5849625 | Median_logFC > 0.5849625))
  
  return(list(all = Results_all, sig = Results_sig))
}


#______________RUNNING LIMMA FOR 4 DEFINED BIOPSY CONDITIONS__________________

#___________1. PEX BASAL

# Filter long data for Basal clamp and Exercise leg
Phos_PEx_basal_long <- Phos_long_base %>%
  filter(Clamp != "Insulin", Leg != "Rest")

# Run pipeline
PEx_basal_out <- run_limma_pipeline(Phos_PEx_basal_long, Sample_key_RAPA)
Results_PEx_basal <- PEx_basal_out$all
Results_PEx_basal_sig <- PEx_basal_out$sig


#__________ 2. DELTA PEX BASAL

# Filter long data for Basal clamp and calculate Exercise response (Exercise - Rest)
Phos_Delta_PEx_basal_long <- Phos_long_base %>%
  filter(Clamp != "Insulin") %>%
  group_by(Subject, Site, Drug) %>%
  mutate(Phos = Phos[Leg == "Exercise"] - Phos[Leg == "Rest"]) %>%
  ungroup() %>%
  filter(Leg != "Rest")

# Run pipeline
Delta_PEx_basal_out <- run_limma_pipeline(Phos_Delta_PEx_basal_long, Sample_key_RAPA)
Results_Delta_PEx_basal <- Delta_PEx_basal_out$all
Results_Delta_PEx_basal_sig <- Delta_PEx_basal_out$sig


#___________ 3. PEX INSULIN

# Filter long data for Insulin clamp and Exercise leg
Phos_PEx_insulin_long <- Phos_long_base %>%
  filter(Clamp != "Basal", Leg != "Rest")

# Run pipeline
PEx_insulin_out <- run_limma_pipeline(Phos_PEx_insulin_long, Sample_key_RAPA)
Results_PEx_insulin <- PEx_insulin_out$all
Results_PEx_insulin_sig <- PEx_insulin_out$sig


#___________ 4. DELTA PEX INSULIN

# Filter long data for Insulin clamp and calculate Exercise response (Exercise - Rest)
Phos_Delta_PEx_insulin_long <- Phos_long_base %>%
  filter(Clamp != "Basal") %>%
  group_by(Subject, Site, Drug) %>%
  mutate(Phos = Phos[Leg == "Exercise"] - Phos[Leg == "Rest"]) %>%
  ungroup() %>%
  filter(Leg != "Rest")

# Run pipeline
Delta_PEx_insulin_out <- run_limma_pipeline(Phos_Delta_PEx_insulin_long, Sample_key_RAPA)
Results_Delta_PEx_insulin <- Delta_PEx_insulin_out$all
Results_Delta_PEx_insulin_sig <- Delta_PEx_insulin_out$sig







#______________________________ PERSONALIZED PHOSPHO STEP 2 - EXERCISE x INSULIN EFFECT _____________________________

# Prepare expression matrix and ensure samples align with metadata
Phos_matrix_all <- as.matrix(Phospho_data_filtered)
Design_all <- Sample_key_RAPA %>%
  filter(SampleID %in% colnames(Phos_matrix_all)) %>%
  arrange(match(SampleID, colnames(Phos_matrix_all)))

# Create combined factor and individual block ID (Subject x Drug)
Design_all$Leg.Clamp <- factor(paste(Design_all$Leg, Design_all$Clamp, sep = "."))
Design_all$Subject_Day <- factor(paste(Design_all$Subject, Design_all$Drug, sep = "_"))

# Construct linear model design matrix and contrasts
Design_matrix_all <- model.matrix(~ 0 + Leg.Clamp, data = Design_all)

Contrasts_all <- makeContrasts(
  LegEffect = (Leg.ClampExercise.Basal + Leg.ClampExercise.Insulin)/2 - 
    (Leg.ClampRest.Basal + Leg.ClampRest.Insulin)/2,
  ClampEffect = (Leg.ClampExercise.Insulin + Leg.ClampRest.Insulin)/2 - 
    (Leg.ClampExercise.Basal + Leg.ClampRest.Basal)/2,
  Interaction = (Leg.ClampExercise.Insulin - (Leg.ClampRest.Insulin - Leg.ClampRest.Basal) - (Leg.ClampExercise.Basal - Leg.ClampRest.Basal) - Leg.ClampRest.Basal),
  levels = Design_matrix_all
)

# Run limma with Subject_Day blocking for repeated measures across legs and clamp states
corfit_all <- duplicateCorrelation(Phos_matrix_all, Design_matrix_all, block = Design_all$Subject_Day)
fit_all <- lmFit(Phos_matrix_all, Design_matrix_all, block = Design_all$Subject_Day, correlation = corfit_all$consensus)
fit_all <- contrasts.fit(fit_all, Contrasts_all)
fit_all <- eBayes(fit_all)

# Extract statistics for main effects and interaction
results_Exercise <- topTable(fit_all, coef = "LegEffect", number = Inf, genelist = rownames(Phos_matrix_all))
results_Insulin <- topTable(fit_all, coef = "ClampEffect", number = Inf, genelist = rownames(Phos_matrix_all))
results_Interaction <- topTable(fit_all, coef = "Interaction", number = Inf, genelist = rownames(Phos_matrix_all))

# Merge all contrast statistics into a single data frame
results_Exercise_sub <- results_Exercise %>% 
  dplyr::select(ID, logFC_Exercise = logFC, AveExpr_Exercise = AveExpr, t_Exercise = t, P.Value_Exercise = P.Value, adj.P.Val_Exercise = adj.P.Val, B_Exercise = B)

results_Insulin_sub <- results_Insulin %>% 
  dplyr::select(ID, logFC_Insulin = logFC, AveExpr_Insulin = AveExpr, t_Insulin = t, P.Value_Insulin = P.Value, adj.P.Val_Insulin = adj.P.Val, B_Insulin = B)

results_Interaction_sub <- results_Interaction %>% 
  dplyr::select(ID, logFC_Interaction = logFC, AveExpr_Interaction = AveExpr, t_Interaction = t, P.Value_Interaction = P.Value, adj.P.Val_Interaction = adj.P.Val, B_Interaction = B)

Results_exins_all <- results_Exercise_sub %>%
  inner_join(results_Insulin_sub, by = "ID") %>%
  inner_join(results_Interaction_sub, by = "ID")

# Filter for phosphosites regulated by Exercise, Insulin, or Interaction (p < 0.05)
Results_exins_sig <- Results_exins_all %>%
  filter(P.Value_Exercise < 0.05 | P.Value_Insulin < 0.05 | P.Value_Interaction < 0.05)






#______________________________________ PERSONALIZED PHOSPHO STEP 3 - PHOSPHO VS. GLUCOSE CORRELATION ____________________________ 

#_______FILTER SIGNIFICANT SITES AT BOTH RAPAMYCIN AND EXERCISE/INSULIN LEVEL __________

# 1. PEx Basal
Results_PEx_basal_sig_exins <- Results_PEx_basal_sig %>%
  inner_join(Results_exins_sig, by = "ID")

# 2. Delta PEx Basal
Results_Delta_PEx_basal_sig_exins <- Results_Delta_PEx_basal_sig %>%
  inner_join(Results_exins_sig, by = "ID")

# 3. PEx Insulin
Results_PEx_insulin_sig_exins <- Results_PEx_insulin_sig %>%
  inner_join(Results_exins_sig, by = "ID")

# 4. Delta PEx Insulin
Results_Delta_PEx_insulin_sig_exins <- Results_Delta_PEx_insulin_sig %>%
  inner_join(Results_exins_sig, by = "ID")






#____________PREPARE GLUCOSE UPTAKE DATA _____________________________

# Calculate delta exercise glucose uptake (Exercise - Rest) per Subject, Drug, and Clamp
corr_GU <- Sample_key_RAPA %>%
  group_by(Subject, Drug, Clamp) %>%
  mutate(GU_AUC = GU_AUC[Leg == "Exercise"] - GU_AUC[Leg == "Rest"]) %>%
  ungroup() %>%
  filter(Leg == "Exercise")


#____________PREPARE DATASETS FOR 4 BIOPSY CONDITIONS __________________

# Transpose filtered phospho matrix into sample rows
df_Phos_wide <- Phospho_data_filtered %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID")

# 1. PEx Basal
b_melted_Phos <- corr_GU %>%
  filter(Clamp != "Insulin") %>%
  inner_join(df_Phos_wide, by = "SampleID") %>%
  reshape2::melt(
    id.vars = c("SampleID", "Subject", "Drug", "Leg", "Clamp", "GU_AUC", "GU_MEAN"),
    variable.name = "Site",
    value.name = "Phos"
  ) %>%
  mutate(Phos = as.numeric(Phos)) %>%
  filter(Site %in% Results_PEx_basal_sig_exins$ID) %>%
  mutate(Site = paste0("b_", Site))

# 2. Delta PEx Basal (Exercise - Rest for Phos)
db_melted_Phos <- Phospho_data_filtered %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  inner_join(Sample_key_RAPA, by = "SampleID") %>%
  reshape2::melt(
    id.vars = c("SampleID", "Subject", "Drug", "Leg", "Clamp", "GU_AUC", "GU_MEAN"),
    variable.name = "Site",
    value.name = "Phos"
  ) %>%
  mutate(Phos = as.numeric(Phos)) %>%
  filter(Clamp != "Insulin") %>%
  group_by(Subject, Site, Drug) %>%
  mutate(Phos = Phos[Leg == "Exercise"] - Phos[Leg == "Rest"]) %>%
  ungroup() %>%
  filter(Leg == "Exercise") %>%
  inner_join(corr_GU %>% dplyr::select(SampleID, GU_delta = GU_AUC), by = "SampleID") %>%
  mutate(GU_AUC = GU_delta) %>%
  dplyr::select(-GU_delta) %>%
  filter(Site %in% Results_Delta_PEx_basal_sig_exins$ID) %>%
  mutate(Site = paste0("db_", Site))

# 3. PEx Insulin
i_melted_Phos <- corr_GU %>%
  filter(Clamp != "Basal") %>%
  inner_join(df_Phos_wide, by = "SampleID") %>%
  reshape2::melt(
    id.vars = c("SampleID", "Subject", "Drug", "Leg", "Clamp", "GU_AUC", "GU_MEAN"),
    variable.name = "Site",
    value.name = "Phos"
  ) %>%
  mutate(Phos = as.numeric(Phos)) %>%
  filter(Site %in% Results_PEx_insulin_sig_exins$ID) %>%
  mutate(Site = paste0("i_", Site))

# 4. Delta PEx Insulin (Exercise - Rest for Phos)
di_melted_Phos <- Phospho_data_filtered %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  inner_join(Sample_key_RAPA, by = "SampleID") %>%
  reshape2::melt(
    id.vars = c("SampleID", "Subject", "Drug", "Leg", "Clamp", "GU_AUC", "GU_MEAN"),
    variable.name = "Site",
    value.name = "Phos"
  ) %>%
  mutate(Phos = as.numeric(Phos)) %>%
  filter(Clamp != "Basal") %>%
  group_by(Subject, Site, Drug) %>%
  mutate(Phos = Phos[Leg == "Exercise"] - Phos[Leg == "Rest"]) %>%
  ungroup() %>%
  filter(Leg == "Exercise") %>%
  inner_join(corr_GU %>% dplyr::select(SampleID, GU_delta = GU_AUC), by = "SampleID") %>%
  mutate(GU_AUC = GU_delta) %>%
  dplyr::select(-GU_delta) %>%
  filter(Site %in% Results_Delta_PEx_insulin_sig_exins$ID) %>%
  mutate(Site = paste0("di_", Site))


#______________ RUN REPEATED MEASURES CORRELATION (rmcorr) ________________

# Combine the 4 condition datasets and clean SampleID suffixes
melted_Phos_all <- bind_rows(b_melted_Phos, db_melted_Phos, i_melted_Phos, di_melted_Phos) %>%
  mutate(SampleID = sub("_(4hPEX\\+2hIns|4hPEX)$", "", SampleID))

# Extract unique phosphosites
unique_sites <- unique(melted_Phos_all$Site)

# Set up parallel cluster
cl <- makeCluster(4)
registerDoParallel(cl)

# Run repeated measures correlation with 1,000-iteration subject-level bootstrap
correlation_results <- foreach(
  curr_site = unique_sites, 
  .packages = c("rmcorr"), 
  .combine = rbind
) %dopar% {
  
  # Subset and clean complete cases for current site
  subset_data <- melted_Phos_all[melted_Phos_all$Site == curr_site & !is.na(melted_Phos_all$Phos), ]
  subset_data$Subject <- as.character(subset_data$Subject)
  
  # Evaluate threshold (minimum 13 phos observations)
  if (nrow(subset_data) >= 13) {
    
    # Standard rmcorr execution
    orig_res <- tryCatch({
      rmcorr(participant = Subject, measure1 = GU_AUC, measure2 = Phos, dataset = subset_data)
    }, error = function(e) NULL)
    
    if (!is.null(orig_res)) {
      
      unique_subs <- unique(subset_data$Subject)
      boot_r_values <- numeric(1000)
      
      # 1,000 subject-level resamples
      for (i in 1:1000) {
        samp_subs <- sample(unique_subs, replace = TRUE)
        
        boot_list <- lapply(seq_along(samp_subs), function(idx) {
          temp <- subset_data[subset_data$Subject == samp_subs[idx], ]
          temp$Subject <- paste0("S_", idx)
          return(temp)
        })
        boot_df <- do.call(rbind, boot_list)
        
        boot_r_values[i] <- tryCatch({
          rmcorr(Subject, GU_AUC, Phos, boot_df)$r
        }, error = function(e) NA)
      }
      
      # Clean bootstrap distribution and derive 95% CI
      valid_boot_r <- boot_r_values[!is.na(boot_r_values)]
      
      if (length(valid_boot_r) > 950) {
        ci_lower <- quantile(valid_boot_r, 0.025)
        ci_upper <- quantile(valid_boot_r, 0.975)
      } else {
        ci_lower <- NA
        ci_upper <- NA
      }
      
      data.frame(
        Site = curr_site, 
        r = orig_res$r, 
        p = orig_res$p, 
        CI_lower = ci_lower, 
        CI_upper = ci_upper
      )
      
    } else {
      data.frame(Site = curr_site, r = NA, p = NA, CI_lower = NA, CI_upper = NA)
    }
    
  } else {
    data.frame(Site = curr_site, r = NA, p = NA, CI_lower = NA, CI_upper = NA)
  }
}

# Terminate parallel cluster
stopCluster(cl)

# Classify statistical significance and confidence interval consistency
correlation_results <- correlation_results %>%
  mutate(Sign = ifelse(
    !is.na(p) & p < 0.05 & (CI_lower > 0 | CI_upper < 0), 
    "sign", 
    "non-sign"
  ))







#__________________ ADDITIONAL PERSONALIZED PHOSPHO (STEP 4) (GLUCOSE CORRELATION FROM REST TO EXERCISE IN PLACEBO) ___________________

# Filter for Placebo and Insulin, keeping GU_MEAN across both legs (Rest & Exercise)
corr_GU <- Sample_key_RAPA %>%
  filter(Drug == "Placebo", Clamp == "Insulin") %>%
  select(Subject, Leg, GU_MEAN) %>%
  distinct()

#___________ PREPARE PHOSPHO DATA __________________

# Define target sites
target_sites <- c("MKNK2_S74_M1", "EIF4G1_S1124_M1", "EIF4G1_S1209_M1")

# Filter Phos_long_base for Placebo, Insulin, and target sites, then join GU_MEAN
melted_Phos_target <- Phos_long_base %>%
  filter(
    Drug == "Placebo", 
    Clamp == "Insulin", 
    Site %in% target_sites
  ) %>%
  select(Subject, Leg, Site, Phos) %>%
  inner_join(corr_GU, by = c("Subject", "Leg"))

#_________ RUN REPEATED MEASURES CORRELATION (rmcorr) ___________

# Set up parallel cluster
cl <- makeCluster(4)
registerDoParallel(cl)

# Run repeated measures correlation with 1,000-iteration subject-level bootstrap
correlation_results_exercise <- foreach(
  curr_site = target_sites, 
  .packages = c("rmcorr"), 
  .combine = rbind
) %dopar% {
  
  # Subset and clean complete cases for current site
  subset_data <- melted_Phos_target[melted_Phos_target$Site == curr_site & !is.na(melted_Phos_target$Phos), ]
  subset_data$Subject <- as.character(subset_data$Subject)
  
  # Standard rmcorr execution
  orig_res <- tryCatch({
    rmcorr(participant = Subject, measure1 = GU_MEAN, measure2 = Phos, dataset = subset_data)
  }, error = function(e) NULL)
  
  if (!is.null(orig_res)) {
    
    unique_subs <- unique(subset_data$Subject)
    boot_r_values <- numeric(1000)
    
    # 1,000 subject-level resamples
    for (i in 1:1000) {
      samp_subs <- sample(unique_subs, replace = TRUE)
      
      boot_list <- lapply(seq_along(samp_subs), function(idx) {
        temp <- subset_data[subset_data$Subject == samp_subs[idx], ]
        temp$Subject <- paste0("S_", idx)
        return(temp)
      })
      boot_df <- do.call(rbind, boot_list)
      
      boot_r_values[i] <- tryCatch({
        rmcorr(Subject, GU_MEAN, Phos, boot_df)$r
      }, error = function(e) NA)
    }
    
    # Clean bootstrap distribution and derive 95% CI
    valid_boot_r <- boot_r_values[!is.na(boot_r_values)]
    
    if (length(valid_boot_r) > 950) {
      ci_lower <- quantile(valid_boot_r, 0.025)
      ci_upper <- quantile(valid_boot_r, 0.975)
    } else {
      ci_lower <- NA
      ci_upper <- NA
    }
    
    data.frame(
      Site = curr_site, 
      r = orig_res$r, 
      p = orig_res$p, 
      CI_lower = ci_lower, 
      CI_upper = ci_upper
    )
    
  } else {
    data.frame(Site = curr_site, r = NA, p = NA, CI_lower = NA, CI_upper = NA)
  }
}

# Terminate parallel cluster
stopCluster(cl)

# Classify statistical significance and confidence interval consistency
correlation_results_exercise <- correlation_results_exercise %>%
  mutate(Sign = ifelse(
    !is.na(p) & p < 0.05 & (CI_lower > 0 | CI_upper < 0), 
    "sign", 
    "non-sign"
  ))






#__________________ ADDITIONAL PERSONALIZED PHOSPHO (STEP 5) (RE-ANALYSIS OF EIF4G1 S1209) ___________________

# Calculate delta exercise glucose uptake (Exercise - Rest) per Subject, Drug, and Clamp
corr_GU <- Sample_key_RAPA %>%
  group_by(Subject, Drug, Clamp) %>%
  mutate(GU_AUC = GU_AUC[Leg == "Exercise"] - GU_AUC[Leg == "Rest"]) %>%
  ungroup() %>%
  filter(Leg == "Exercise")

#__________ PREPARE PEx INSULIN DATA FOR EIF4G1_S1209_M1 ____________

# Extract EIF4G1_S1209_M1 from Phos_long_base for PEx Insulin (Exercise leg, Insulin clamp)
target_site <- "EIF4G1_S1209_M1"

subset_data <- Phos_long_base %>%
  filter(
    Site == target_site,
    Clamp == "Insulin",
    Leg == "Exercise",
    !is.na(Phos)
  ) %>%
  select(Subject, Drug, Clamp, Site, Phos) %>%
  inner_join(
    corr_GU %>% select(Subject, Drug, Clamp, GU_AUC), 
    by = c("Subject", "Drug", "Clamp")
  ) %>%
  mutate(Subject = as.character(Subject))

#_____________________ RUN RMCORR ______________

set.seed(123) # For reproducible bootstrap sampling

# 1. Standard rmcorr execution
orig_res <- tryCatch({
  rmcorr(participant = Subject, measure1 = GU_AUC, measure2 = Phos, dataset = subset_data)
}, error = function(e) NULL)

if (!is.null(orig_res)) {
  
  unique_subs <- unique(subset_data$Subject)
  boot_r_values <- numeric(1000)
  
  # 2. 1,000 subject-level resamples
  for (i in 1:1000) {
    samp_subs <- sample(unique_subs, replace = TRUE)
    
    boot_list <- lapply(seq_along(samp_subs), function(idx) {
      temp <- subset_data[subset_data$Subject == samp_subs[idx], ]
      temp$Subject <- paste0("S_", idx)
      return(temp)
    })
    boot_df <- do.call(rbind, boot_list)
    
    boot_r_values[i] <- tryCatch({
      rmcorr(Subject, GU_AUC, Phos, boot_df)$r
    }, error = function(e) NA)
  }
  
  # 3. Clean bootstrap distribution and derive 95% CI
  valid_boot_r <- boot_r_values[!is.na(boot_r_values)]
  
  if (length(valid_boot_r) > 950) {
    ci_lower <- quantile(valid_boot_r, 0.025)
    ci_upper <- quantile(valid_boot_r, 0.975)
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }
  
  # 4. Construct final summary result
  correlation_result_EIF4G1 <- data.frame(
    Site = paste0("i_", target_site),
    r = orig_res$r,
    df = orig_res$df,
    p = orig_res$p,
    CI_lower = ci_lower,
    CI_upper = ci_upper
  ) %>%
    mutate(Sign = ifelse(
      !is.na(p) & p < 0.05 & (CI_lower > 0 | CI_upper < 0), 
      "sign", 
      "non-sign"
    ))
  
} else {
  correlation_result_EIF4G1 <- data.frame(
    Site = paste0("i_", target_site),
    r = NA,
    df = NA,
    p = NA,
    CI_lower = NA,
    CI_upper = NA,
    Sign = "non-sign"
  )
}







#_____________________________ VISUALISATION_____________________________________________________

#__________________ Z-SCORED RAPAMYCIN HEATMAP PIPELINE ____________


# Function to extract, directional-adjust, and z-score paired Rapamycin response per condition
get_zscored_rapa_deltas <- function(long_data, sig_results, prefix) {
  
  # 1. Calculate paired Rapamycin delta (Rapamycin - Placebo)
  delta_mat <- long_data %>%
    filter(Site %in% sig_results$ID) %>%
    group_by(Subject, Site) %>%
    mutate(Delta_Phos = Phos[Drug == "Rapamycin"] - Phos[Drug == "Placebo"]) %>%
    ungroup() %>%
    filter(Drug != "Placebo") %>%
    reshape2::dcast(Subject ~ Site, value.var = "Delta_Phos") %>%
    column_to_rownames("Subject") %>%
    as.matrix()
  
  # 2. Adjust columns: if median response < 0, invert (multiply by -1) for absolute magnitude
  col_medians <- apply(delta_mat, 2, median, na.rm = TRUE)
  neg_cols <- which(!is.na(col_medians) & col_medians < 0)
  if (length(neg_cols) > 0) {
    delta_mat[, neg_cols] <- delta_mat[, neg_cols] * -1
  }
  
  # 3. Z-score across subjects per site (columns) and transpose to Sites x Subjects
  z_mat <- scale(delta_mat) %>%
    t() %>%
    as.matrix()
  
  # 4. Attach condition-specific prefix to rownames
  rownames(z_mat) <- paste0(prefix, "_", rownames(z_mat))
  
  return(z_mat)
}


#______ PROCESS THE 4 BIOPSY CONDITIONS _____

# 1. PEx Basal
z_PEx_basal <- get_zscored_rapa_deltas(Phos_PEx_basal_long, Results_PEx_basal_sig, "b")

# 2. Delta PEx Basal
z_Delta_PEx_basal <- get_zscored_rapa_deltas(Phos_Delta_PEx_basal_long, Results_Delta_PEx_basal_sig, "db")

# 3. PEx Insulin
z_PEx_insulin <- get_zscored_rapa_deltas(Phos_PEx_insulin_long, Results_PEx_insulin_sig, "i")

# 4. Delta PEx Insulin
z_Delta_PEx_insulin <- get_zscored_rapa_deltas(Phos_Delta_PEx_insulin_long, Results_Delta_PEx_insulin_sig, "di")

# Combine all 4 conditions into one master Sites x Subjects matrix
Phos_rapa_z_master <- rbind(z_PEx_basal, z_Delta_PEx_basal, z_PEx_insulin, z_Delta_PEx_insulin)


#____________ GENERATE HEATMAP ______

# Define color ramp for Z-scores (-2 to +2)
col_fun_zscore <- colorRamp2(c(-2, 0, 2), c("green", "yellow", "red"))

# Generate clustered heatmap without legend
ht_zscore <- Heatmap(
  Phos_rapa_z_master, 
  name = "Z-score",
  col = col_fun_zscore,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_distance_columns = "pearson",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_split = 10,
  row_title = NULL,
  row_gap = unit(1.5, "mm"),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 0.8),
  show_heatmap_legend = FALSE, # Suppress legend
  na_col = "white"
)

# Draw heatmap
draw(ht_zscore)




#______________ CORRELATION MATRIX OF PERSONALIZED PHOSPHO SITES ___________________


# Extract significant phosphosites passing both limma and rmcorr filtering
sig_sites <- correlation_results$Site[correlation_results$Sign == "sign"]

# Reshape merged long data to wide format for significant sites
phos_wide <- melted_Phos_all %>%
  dplyr::filter(Site %in% sig_sites) %>%
  dplyr::select(Subject, SampleID, Site, Phos) %>%
  tidyr::pivot_wider(names_from = Site, values_from = Phos)

# Initialize symmetric pairwise correlation matrix
n_sites <- length(sig_sites)
rmcorr_matrix <- matrix(
  NA, 
  nrow = n_sites, 
  ncol = n_sites, 
  dimnames = list(sig_sites, sig_sites)
)

# Compute pairwise absolute repeated measures correlations
for (i in 1:n_sites) {
  for (j in 1:n_sites) {
    if (i == j) {
      rmcorr_matrix[i, j] <- 1
    } else if (i < j) {
      site1 <- sig_sites[i]
      site2 <- sig_sites[j]
      
      temp_data <- phos_wide[, c("Subject", site1, site2)]
      colnames(temp_data) <- c("Subject", "Measure1", "Measure2")
      temp_data <- temp_data[complete.cases(temp_data), ]
      
      if (length(unique(temp_data$Subject)) >= 3) {
        res <- tryCatch({
          rm_val <- rmcorr(participant = Subject, measure1 = Measure1, measure2 = Measure2, dataset = temp_data)$r
          abs(rm_val)
        }, error = function(e) NULL)
        
        if (!is.null(res)) {
          rmcorr_matrix[i, j] <- res
          rmcorr_matrix[j, i] <- res
        }
      }
    }
  }
}

# 6. Assign condition colors based on original prefix before cleaning
raw_labels <- rownames(rmcorr_matrix)
row_colors <- case_when(
  grepl("^db_", raw_labels) ~ "purple",
  grepl("^b_",  raw_labels) ~ "tan3",
  grepl("^di_", raw_labels) ~ "#D81B60",
  grepl("^i_",  raw_labels) ~ "forestgreen",
  TRUE                      ~ "black"
)

# 5. Clean labels by stripping condition prefixes, multiplicity tags, and underscores
clean_labels <- raw_labels %>%
  sub("^(db_|di_|b_|i_)", "", .) %>%
  sub("(_M1|_M2|_M3)", "", .) %>%
  gsub("_", " ", .)


# Define color palette for absolute correlation scale (0 to 1)
col_fun_rmcorr <- colorRamp2(c(0, 0.5, 1), c("white", "moccasin", "red"))

# Generate clustered correlation heatmap without any borders/lines
ht <- Heatmap(
  rmcorr_matrix, 
  name = "Abs rmcorr", 
  col = col_fun_rmcorr,
  rect_gp = gpar(col = NA),                 # No cell borders
  border = FALSE,                           # No slice borders
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,                    # Remove side dendrogram
  show_column_dend = TRUE,                  # Keep top dendrogram
  column_dend_height = unit(5, "mm"),       # 1/3 height
  row_split = 6,
  column_split = 6,
  row_gap = unit(1, "mm"),                  # 1.0 mm gap between rows
  column_gap = unit(1, "mm"),               # 1.0 mm gap between columns
  row_title = NULL,                         # Remove cluster numbers
  column_title = NULL,                      # Remove cluster numbers
  show_column_names = FALSE,
  row_labels = clean_labels,
  row_names_side = "left",                  # Labels on left
  row_names_gp = gpar(
    fontsize = 5.5, 
    fontface = "bold",                      # Bold text
    col = row_colors
  ),
  heatmap_legend_param = list(
    title = "Absolute\nrmcorr (r)", 
    title_gp = gpar(fontsize = 6, fontface = "bold"), 
    labels_gp = gpar(fontsize = 5.5),
    grid_width = unit(2.5, "mm"),           # 1/2 size
    legend_height = unit(15, "mm")
  ),
  na_col = "grey90"
)

# Draw heatmap with centered right legend
draw(ht, heatmap_legend_side = "right")






#______________________________ EXPORT DATA ___________________________________________________

#________ EXPORT ALL RAPAMYCIN RESULTS & Z-SCORED DATA TO A SINGLE EXCEL FILE __________

# Helper function to compute -log10(P.Value), assign significance flag, and select specific columns
format_results_df <- function(df, sig_df) {
  # 1. Calculate -log10 of P.Value
  df$neg_log10_P.Value <- -log10(df$P.Value)
  
  # 2. Add significance flag: "+" if present in sig_df$ID, else ""
  df$`significant threshold` <- ifelse(df$ID %in% sig_df$ID, "+", "")
  
  # 3. Subset and order the specific columns
  target_cols <- c("ID", "Median_logFC", "P.Value", "neg_log10_P.Value", "qvalue", "significant threshold")
  df <- df[, target_cols]
  
  return(df)
}

# Update the 4 data frames directly in place
Results_PEx_basal_export         <- format_results_df(Results_PEx_basal, Results_PEx_basal_sig)
Results_Delta_PEx_basal_export   <- format_results_df(Results_Delta_PEx_basal, Results_Delta_PEx_basal_sig)
Results_PEx_insulin_export       <- format_results_df(Results_PEx_insulin, Results_PEx_insulin_sig)
Results_Delta_PEx_insulin_export <- format_results_df(Results_Delta_PEx_insulin, Results_Delta_PEx_insulin_sig)

# Prepare z-score data frame by converting rownames into an explicit "ID" column
z_scored_export <- as.data.frame(Phos_rapa_z_master) %>%
  rownames_to_column(var = "ID")

# Compile all significant rapamycin sites across the 4 conditions
Rapa_sign <- bind_rows(
  Results_PEx_basal_sig %>% 
    dplyr::select(ID, Median_logFC) %>% 
    mutate(Condition = "PEx_Basal"),
  
  Results_Delta_PEx_basal_sig %>% 
    dplyr::select(ID, Median_logFC) %>% 
    mutate(Condition = "Delta_PEx_Basal"),
  
  Results_PEx_insulin_sig %>% 
    dplyr::select(ID, Median_logFC) %>% 
    mutate(Condition = "PEx_Insulin"),
  
  Results_Delta_PEx_insulin_sig %>% 
    dplyr::select(ID, Median_logFC) %>% 
    mutate(Condition = "Delta_PEx_Insulin")
) %>%
  dplyr::select(ID, Condition, Median_logFC)

# Export all 6 sheets directly to Desktop
write.xlsx(
  list(
    "PEx_Basal"                = Results_PEx_basal_export,
    "Delta_PEx_Basal"          = Results_Delta_PEx_basal_export,
    "PEx_Insulin"              = Results_PEx_insulin_export,
    "Delta_PEx_Insulin"        = Results_Delta_PEx_insulin_export,
    "z-scored rapamycin sites" = z_scored_export,
    "Rapa_sign"                = Rapa_sign
  ),
  file = "~/Desktop/Table_S6_invivo_rapamycin_sites.xlsx",
  rowNames = FALSE
)


#________ EXPORT EXINS SITES TO EXCEL __________

# 1. Define target columns in the specified order
target_cols_exins <- c(
  "ID", 
  "logFC_Exercise", "P.Value_Exercise", "adj.P.Val_Exercise", 
  "logFC_Insulin", "P.Value_Insulin", "adj.P.Val_Insulin", 
  "logFC_Interaction", "P.Value_Interaction", "adj.P.Val_Interaction"
)

# 2. Subset and order columns for the "all sites" sheet
Results_exins_formatted <- Results_exins_all[, target_cols_exins]

# 3. Combine rapamycin-significant sites with their respective condition tags
rapa_sig_by_condition <- bind_rows(
  data.frame(ID = Results_PEx_basal_sig$ID,         Condition = "PEx basal"),
  data.frame(ID = Results_Delta_PEx_basal_sig$ID,   Condition = "dPEx basal"),
  data.frame(ID = Results_PEx_insulin_sig$ID,       Condition = "PEx insulin"),
  data.frame(ID = Results_Delta_PEx_insulin_sig$ID, Condition = "dPEx insulin")
)

# 4. Join with exercise/insulin statistics (sites significant in multiple conditions are listed per condition)
Results_exins_rapa_sig <- rapa_sig_by_condition %>%
  inner_join(Results_exins_formatted, by = "ID") %>%
  dplyr::select(ID, Condition, all_of(target_cols_exins[-1]))

# 5. Export both sheets directly to Desktop
write.xlsx(
  list(
    "all sites"       = Results_exins_formatted,
    "rapamycin sites" = Results_exins_rapa_sig
  ),
  file = "~/Desktop/Table_S7_invivo_exins_sites.xlsx",
  rowNames = FALSE
)



#________ EXPORT CORRELATION RESULTS TO EXCEL __________
# Export the three correlation result data frames to Desktop
write.xlsx(
  list(
    "rapamycin effect"                = correlation_results,
    "exercise effect"                 = correlation_results_exercise,
    "EIF4G1 S1209 rapamycin effect"   = correlation_result_EIF4G1
  ),
  file = "~/Desktop/Table_S8_invivo_correlation.xlsx",
  rowNames = FALSE
)





