# code-files
Analysis code for the manuscript: A Non-Invasive Urinary Diagnostic Signature for Diabetic Kidney Disease Revealed by Machine Learning and Single-Cell Analysis
# Non-Invasive Urinary Diagnostic Signature for Diabetic Kidney Disease

This repository contains the R code and analysis pipeline for the study: **"A Non-Invasive Urinary Diagnostic Signature for Diabetic Kidney Disease Revealed by Machine Learning and Single-Cell Analysis"**.

## 📋 Project Overview

This study integrates urinary and renal single-cell sequencing with machine learning to identify a three-gene biomarker panel (PDK4, RHCG, FBP1) for non-invasive diagnosis of diabetic kidney disease (DKD).

## 🛠️ Environment Setup

### Required R Packages

Please install the following R packages with specified versions:

```r
# Core analysis packages
install.packages("Seurat")        # v5.3.0
install.packages("caret")         # v7.0.1
install.packages("glmnet")        # v4.1.9
install.packages("pROC")          # v1.18.5
install.packages("randomForest")  # v4.7.1.2
install.packages("kernlab")       # v0.9.33
install.packages("xgboost")       # v1.7.11.1
install.packages("DALEX")         # For model explanation
install.packages("ggplot2")       # For visualization

# Additional dependencies
install.packages("harmony")       # v1.2.3 for batch correction

 Data Requirements
🔴 IMPORTANT: Download Required Data First!
Before running any scripts, you MUST download the following data files from Figshare:

Figshare Repository: https://doi.org/10.6084/m9.figshare.30218608

Download and place these files in the data/raw/ directory:

data.train.txt - Training dataset for machine learning models

data.test.txt - Testing dataset for model validation

interGenes.txt - LASSO-selected genes (8-gene panel)

GEO Datasets Used
This study utilizes the following public datasets (automatically downloaded by scripts):

GSE96804: Bulk RNA-seq (40 DKD vs 21 controls)

GSE104948/GSE104954: Validation cohort (30 DKD vs 42 controls)

GSE142025: Independent validation (27 DKD vs 9 controls)

GSE131882: DKD kidney tissue scRNA-seq (3 DKD vs 3 controls)

GSE266146: DKD urine sediment scRNA-seq (8 DKD patients)

GSE157640: Healthy urine sediment scRNA-seq (10 controls)

🚀 Reproduction Steps
Step 1: Environment Preparation
Install all required R packages (see Environment Setup above)

Download essential data files from Figshare to data/raw/

Set working directory to project root

Step 2: Data Preprocessing
r
source("scripts/01_data_preprocessing.R")
Loads and preprocesses training and testing data

Applies quality control and normalization

Extracts LASSO-selected genes

Step 3: Single-Cell Analysis
r
source("scripts/02_single_cell_analysis.R")
Processes single-cell RNA sequencing data

Performs cell type annotation and clustering

Identifies differentially expressed genes

Step 4: Machine Learning Modeling
r
source("scripts/03_machine_learning.R")
Trains four machine learning models:

Random Forest (RF)

Support Vector Machine (SVM)

XGBoost (XGB)

Logistic Regression (LR)

Performs cross-validation and model evaluation

Generates ROC curves and performance metrics

Step 5: Results Generation
r
source("scripts/04_results_visualization.R")
Creates publication-quality figures

Generates performance tables

Produces diagnostic model outputs

⚙️ Key Parameters
Machine Learning Parameters
Cross-validation: 5-fold repeated 3 times

Performance metric: ROC-AUC

Feature selection: LASSO regression

Final gene panel: PDK4, RHCG, FBP1

Single-Cell Analysis Parameters
Quality control: nFeature_RNA > 250, nCount_RNA > 500, mt% < 15%

Normalization: LogNormalize (scale factor = 10,000)

Integration: Harmony batch correction

Clustering: 50 principal components, UMAP visualization

📈 Expected Outputs
After successful execution, the following files will be generated:

results/tables/AUC_Results_with_CI2.csv - Model performance metrics

results/figures/ROC_curves.pdf - ROC curves for all models

results/figures/feature_importance.pdf - Feature importance plots

results/tables/model_comparison.csv - Comparative model performance

🐛 Troubleshooting
Common Issues
Data not found error: Ensure all required data files are downloaded from Figshare and placed in data/raw/

Package version conflicts: Use the exact package versions specified above

Memory issues: Single-cell analysis may require 16GB+ RAM

Download failures: Check GEO database accessibility and internet connection
