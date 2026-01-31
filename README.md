# cp_pattern
Repository for exploring the evolutionary relationship between inheritance patterns and traits in organelle DNA of seed plants

## Background / Motivation
This repository contains analysis code used to examine the association between plastid genome inheritance (maternal vs. non-maternal) and plant functional traits across seed plants while accounting for phylogenetic relationships.

## Directory Structure
```text
.
├── input/
│   ├── family_replace.csv  # Family name harmonization table
│   ├── genus_family.csv         # Additional genus–family mapping
│   ├── inheritance_pattern.csv  # Plastid inheritance labels
│   ├── remove_Dataname.csv  # DataName entries to exclude
│   ├── unit_convert_factor.csv  # Unit conversion lookup table
│   └── zanne2014.new       # Phylogenetic tree (Zanne et al. 2014)
├── output/
│   ├── analysis_all_results.RData     # Full analysis workspace after running all analyses
│   ├── bhpmf/                         # Files related to BHPMF-based trait imputation
│   ├── glm_result.rds                 # GLM results (also included in analysis_all_results.RData)
│   ├── glm_result_rm_gymn.rds         # GLM results excluding gymnosperms
│   ├── plot/                          # Figures generated from the analyses
│   ├── tree/                          # Repeatedly generated phylogenetic trees for GLM
│   └── tree_rm_gymn/                  # Phylogenetic trees excluding gymnosperms
├── src/
│   ├── 1_process_data.R                 # Data cleaning and preprocessing
│   ├── 2_mean_target_trait.R            # Calculation of species-level mean trait values
│   ├── 3_make_dataframe_for_analysis.R  # Construction of analysis-ready datasets
│   ├── 4_regression.R                   # Phylogenetic GLM analysis
│   ├── 4_regression_testgymn.R          # Phylogenetic GLM excluding gymnosperms
│   ├── 5_make_figure.R                  # Visualization of main analysis results
│   ├── 5_make_figure_testgymn.R          # Visualization excluding gymnosperms
│   └── utils.R                           # Utility functions used across scripts
└── README.md
