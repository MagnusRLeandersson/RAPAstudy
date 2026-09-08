# Human Skeletal Muscle Phosphoproteomic Analysis Pipeline

This repository contains R scripts for quality control, linear modeling (limma), repeated measures correlation (rmcorr), and cross-study translational integration of human skeletal muscle phosphoproteomics data.

## Pipeline Overview

### In Vivo Pipeline
- **Invivo_QC.R**: Data quality control.
- **Invivo_analyses.R**: Personalized phosphoproteomics workflow to identify downstream effectors of mTORC1 associated with skeletal muscle insulin sensitivity.

### Ex Vivo Pipeline
- **Exvivo_QC.R**: Data quality control.
- **Exvivo_analyses.R**: Linear modeling setup evaluating downstream signaling from MKNK2 (via the inhibitor eFT508) and determining the translational overlap of sites oppositely regulated by rapamycin in the in vivo model.

## Required Input Data
The scripts expect the following Excel source tables (available as Supplementary Data accompanying the published manuscript) located in the active working directory:

- Table_S3_invivo_glucoseuptake.xlsx: In vivo subject characteristics and muscle glucose uptake.
- Table_S4_invivo_phos.xlsx: In vivo phosphoproteomics dataset.
- Table_S6_invivo_rapamycin_sites.xlsx: In vivo rapamycin-regulated phosphosites across different biopsy conditions.
- Table_S9_exvivo_sample_key.xlsx: Ex vivo subject characteristics, experimental layout, and sample metadata.
- Table_S10_exvivo_phos.xlsx: Ex vivo phosphoproteomics dataset.

## Environment & Dependencies
The pipeline was developed and validated in R (version 4.5.2). Exact library versions and dependencies (e.g., limma, rmcorr, ComplexHeatmap, clusterProfiler) are documented in the header of each individual script.

## License
This project is licensed under the MIT License (LICENSE).
