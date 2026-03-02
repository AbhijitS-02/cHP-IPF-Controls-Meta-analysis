# Meta-Analysis of Interstitial Lung Diseases: Chronic Hypersensitivity Pneumonitis (cHP) vs Idiopathic Pulmonary Fibrosis (IPF) vs Controls

## Biomarker Discovery via Integrative Transcriptomic Meta-Analysis and Machine Learning

---

## Table of Contents

1. [Overview](#1-overview)
2. [Study Design & Biological Context](#2-study-design--biological-context)
3. [Datasets](#3-datasets)
4. [Pairwise Comparisons](#4-pairwise-comparisons)
5. [Pipeline Overview](#5-pipeline-overview)
6. [Detailed Methodology](#6-detailed-methodology)
   - [Step 1 — Batch Correction](#step-1--batch-correction-01_batch_correctionr)
   - [Step 2 — TPM Normalisation](#step-2--tpm-normalisation-02_tpmr)
   - [Step 3 — xCell Immune Deconvolution](#step-3--xcell-immune-deconvolution-03_xcellr)
   - [Step 3b — Subsetting to Common Genes](#step-3b--subsetting-to-common-genes-03b_subset_common_genesr)
   - [Step 4 — Correlation Heatmaps](#step-4--correlation-heatmaps-04_correlation_heatmapsr)
   - [Step 5 — mRMR Feature Selection](#step-5--mrmr-feature-selection-05_mrmripynb)
   - [Step 6 — Ridge Regression & GLM Modelling](#step-6--ridge-regression--glm-modelling-06_lasso_refactoredr)
   - [Step 7 — Boxplots & Statistical Testing](#step-7--boxplots--statistical-testing-07_boxplotr)
   - [Step 8 — Independent Validation AUC-ROC](#step-8--independent-validation-auc-roc-08_independent_auc_rocr)
   - [Pairwise xCell Comparison](#pairwise-xcell-comparison-pairwise_comparisonr)
7. [Results Summary](#7-results-summary)
8. [Directory Structure](#8-directory-structure)
9. [How to Reproduce](#9-how-to-reproduce)
10. [Software & Dependencies](#10-software--dependencies)
11. [Key Output Files](#11-key-output-files)
12. [Author](#12-author)
13. [License](#13-license)

---

## 1. Overview

This project implements a **multi-step integrative transcriptomic meta-analysis pipeline** to identify gene expression biomarkers capable of distinguishing between three clinical conditions of interstitial lung disease (ILD):

- **Chronic Hypersensitivity Pneumonitis (cHP)**
- **Idiopathic Pulmonary Fibrosis (IPF)**
- **Healthy Controls**

The pipeline integrates RNA-seq data from two independent GEO datasets (GSE150910 and GSE184316), applies batch correction, performs immune cell deconvolution via xCell, and uses a machine learning–driven feature selection and classification workflow (mRMR → Ridge Regression → GLM) to discover minimal gene signatures with high discriminatory power. Final biomarkers are validated on independent dataset subsets (biopsy and explant tissue from GSE150910).

---

## 2. Study Design & Biological Context

### Clinical Motivation

cHP and IPF are two forms of progressive fibrosing interstitial lung disease that share overlapping clinical and radiological features, making differential diagnosis extremely challenging. However, their aetiology, prognosis, and treatment strategies differ significantly:

- **cHP** is caused by chronic inhalation of environmental antigens, leading to an immune-mediated inflammatory response.
- **IPF** is a progressive fibrotic disease of unknown cause with a median survival of 3–5 years.

Accurate molecular biomarkers could enable earlier and more precise diagnosis, improving patient outcomes.

### Approach

This study uses a **meta-analysis framework** that:

1. Identifies **meta-differentially expressed genes (metaDEGs)** that are consistently dysregulated across multiple independent datasets.
2. Removes platform-specific technical variation using **ComBat batch correction**.
3. Investigates **immune microenvironment differences** using xCell deconvolution.
4. Applies **machine learning feature selection** (mRMR + Ridge/GLM) to discover minimal diagnostic gene panels.
5. Validates candidate biomarkers on **independent tissue subsets** (biopsy vs explant).

---

## 3. Datasets

| Dataset      | GEO Accession | Platform           | Tissue Types          | Conditions Available      |
|:-------------|:--------------|:-------------------|:----------------------|:--------------------------|
| **Dataset 1** | [GSE184316](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184316) | Illumina | Lung tissue | cHP, IPF, Control |
| **Dataset 2** | [GSE150910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150910) | Ion Torrent | Lung biopsy & explant | cHP, IPF, Control |

### Input Data Files

Located in `./data/`:

| Subdirectory | File | Description |
|:-------------|:-----|:------------|
| `GSE184316/` | `GSE184316_raw_counts_HP_Control.csv` | Raw count matrix — cHP vs Control samples |
| `GSE184316/` | `GSE184316_raw_counts_HP_IPF.csv` | Raw count matrix — cHP vs IPF samples |
| `GSE184316/` | `GSE184316_raw_counts_IPF_Control.csv` | Raw count matrix — IPF vs Control samples |
| `GSE184316/` | `Sample_information_HP_Control_GSE184316.csv` | Sample metadata (sample ID, condition) |
| `GSE150910/` | `GSE150910_raw_counts_HP_Control.csv` | Raw count matrix — cHP vs Control samples |
| `GSE150910/` | `GSE150910_raw_counts_HP_IPF.csv` | Raw count matrix — cHP vs IPF samples |
| `GSE150910/` | `GSE150910_raw_counts_IPF_Control.csv` | Raw count matrix — IPF vs Control samples |
| `GSE150910/` | `Sample_information_HP_Control_GSE150910.csv` | Sample metadata |
| `GSE150910/` | `Sample_info_biopsy_HP_Ctrl.csv` | Biopsy-only subset metadata (with `Group_Encoded`) |
| `GSE150910/` | `Sample_info_explant_HP_Ctrl.csv` | Explant-only subset metadata (with `Group_Encoded`) |

### Meta-DEG Input Files

Located in the project root (`./`):

| File | Description | Gene Count |
|:-----|:------------|:-----------|
| `cHP_Ctrl_Common_genes_same_direction_expr.csv` | MetaDEGs for cHP vs Control | ~1,724 genes |
| `cHP_IPF_Common_genes_same_direction_expr.csv` | MetaDEGs for cHP vs IPF | ~389 genes |
| `IPF_Ctrl_Common_genes_same_direction_expr.csv` | MetaDEGs for IPF vs Control | ~1,558 genes |

**Selection criteria for metaDEGs:**
1. Mean absolute log₂FC > 1 across datasets
2. Meta p-adjusted < 0.05
3. Expressed in all datasets
4. Consistent direction of expression (up or down) across all datasets

---

## 4. Pairwise Comparisons

All analyses are conducted for three pairwise conditions:

| Comparison Tag | Condition 1 (Positive) | Condition 2 (Negative/Reference) |
|:---------------|:----------------------|:----------------------------------|
| `cHP_Ctrl` | Chronic Hypersensitivity Pneumonitis (HP) | Healthy Control |
| `cHP_IPF` | Chronic Hypersensitivity Pneumonitis (HP) | Idiopathic Pulmonary Fibrosis (IPF) |
| `IPF_Ctrl` | Idiopathic Pulmonary Fibrosis (IPF) | Healthy Control |

---

## 5. Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    META-ANALYSIS PIPELINE                       │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌─────────┐    ┌──────────┐    ┌──────────┐                   │
│  │GSE184316│    │GSE150910 │    │ MetaDEGs │                   │
│  │Raw Counts│   │Raw Counts│    │(FINM)    │                   │
│  └────┬────┘    └────┬─────┘    └────┬─────┘                   │
│       │              │               │                          │
│       └──────┬───────┘               │                          │
│              ▼                       ▼                          │
│  ┌──────────────────────────────────────────┐                   │
│  │ 01. Merge & Batch Correction (ComBat)   │                   │
│  │     VST → ComBat → PCA visualisation    │                   │
│  └──────────────────┬───────────────────────┘                   │
│                     │                                           │
│          ┌──────────┼──────────┐                                │
│          ▼          ▼          ▼                                │
│  ┌──────────┐ ┌──────────┐ ┌──────────────┐                   │
│  │02. TPM   │ │03b.Subset│ │04. Correlation│                   │
│  │Normalise │ │ Genes    │ │   Heatmaps    │                   │
│  └────┬─────┘ └────┬─────┘ └──────────────┘                   │
│       │             │                                           │
│       ▼             ▼                                           │
│  ┌──────────┐ ┌──────────────────────────────┐                 │
│  │03. xCell │ │05. mRMR Feature Selection    │                 │
│  │Deconvolve│ │    (Python / Jupyter)         │                 │
│  └────┬─────┘ └──────────────┬───────────────┘                 │
│       │                      │                                  │
│       ▼                      ▼                                  │
│  ┌──────────────┐  ┌──────────────────────────────┐            │
│  │ Pairwise     │  │06. Ridge + GLM Modelling     │            │
│  │ xCell Stats  │  │    Gene sweep → optimal k    │            │
│  │ + Bubble Plot│  │    Per-gene ROC + Combo ROC  │            │
│  └──────────────┘  └──────────────┬───────────────┘            │
│                                   │                             │
│                    ┌──────────────┼──────────────┐              │
│                    ▼              ▼              ▼              │
│           ┌──────────────┐ ┌──────────┐ ┌──────────────┐      │
│           │07. Boxplots  │ │Gene Stats│ │08. Indep.    │      │
│           │+ Wilcoxon    │ │  CSV     │ │  AUC-ROC     │      │
│           │  (Train/Test)│ │          │ │  Validation  │      │
│           └──────────────┘ └──────────┘ └──────────────┘      │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 6. Detailed Methodology

### Step 1 — Batch Correction (`01_batch_correction.R`)

**Purpose:** Merge raw count matrices from GSE184316 and GSE150910 and remove platform-specific batch effects while preserving biological variation.

**Method:**
1. **Load metaDEGs** from the `*_Common_genes_same_direction_expr.csv` files and compute their union across all three comparisons.
2. **Subset** each raw count matrix (per dataset, per comparison) to only the union of metaDEG genes.
3. **Merge** the two dataset matrices for each condition (column-bind by matching gene rows).
4. **Check for NA values** after merging (confirmed: none present).
5. **Create DESeq2 objects** with design formula `~ platform + condition`.
6. **VST (Variance Stabilizing Transformation)** with `blind = TRUE` for quality control PCA plots.
7. **PCA visualisation** (before correction):
   - Coloured by **platform** (GSE184316 vs GSE150910) — to visualise batch effects.
   - Coloured by **condition** (HP, IPF, Control) — to visualise biological separation.
   - 95% confidence ellipses added for each group.
8. **ComBat batch correction** (from the `sva` package):
   - Input: VST-transformed expression matrix.
   - Batch variable: `platform` (GSE184316 / GSE150910).
   - Model matrix: `~ condition` — protects biological variation.
   - Parametric empirical Bayes (`par.prior = TRUE`).
9. **Post-correction PCA** — confirms platform effects are removed while biological separation is maintained.
10. **Export** batch-corrected expression matrices as CSV.

**Key Outputs:**
- `results/PCA_plots/` — 24 PCA plots (PNG + PDF, before & after ComBat, by platform & condition)
- `results/Corrected_matrix/` — 3 batch-corrected expression matrices
- `results/Raw_merged_data/` — 3 raw merged expression matrices

---

### Step 2 — TPM Normalisation (`02_tpm.R`)

**Purpose:** Convert raw counts to Transcripts Per Million (TPM) for input to xCell, which requires TPM-normalised expression data.

**Method:**
1. **Load raw count matrices** for all 6 dataset × condition combinations.
2. **Retrieve gene lengths** from Ensembl (via `biomaRt`):
   - Query `hgnc_symbol` → `ensembl_gene_id` → `ensembl_transcript_id` → `transcript_length`.
   - For genes with multiple transcripts (isoforms), keep the **maximum transcript length**.
3. **Remove genes** with no matching Ensembl annotation (e.g., pseudogenes, miRNAs, snoRNAs).
4. **Compute TPM:**
   ```
   TPM = (counts / gene_length_kb) / sum(counts / gene_length_kb) × 10^6
   ```
5. **Export** TPM matrices as tab-separated files for xCell input.

**Key Outputs:**
- `results/xCell_input/` — 6 TPM-normalised matrices (one per dataset × condition)

---

### Step 3 — xCell Immune Deconvolution (`03_xCell.R`)

**Purpose:** Estimate immune and stromal cell type enrichment scores from gene expression data using xCell, then apply batch correction to the enrichment scores.

**Prerequisite:** TPM files from Step 2 must be uploaded to the [xCell web portal](https://xcell.ucsf.edu/) to obtain enrichment scores. Results are downloaded and placed in `results/xCell_results/`.

**Method:**
1. **Load xCell enrichment scores** for each dataset × condition (6 files total).
2. **Merge scores** across datasets for each condition (cell types × samples).
3. **Merge metadata** and verify sample alignment.
4. **Remove zero-variance cell types** (cell types with identical scores across all samples cannot contribute to PCA or statistical analysis).
5. **PCA before batch correction** — coloured by platform and condition.
6. **ComBat batch correction** applied to enrichment scores:
   - Input: filtered enrichment score matrix (cell types × samples).
   - Batch: platform (GSE184316 / GSE150910).
   - Protected: biological condition.
7. **PCA after batch correction** — confirms removal of platform effects.
8. **Export** corrected enrichment scores and PCA plots.

**Key Outputs:**
- `results/xCell_results/` — Merged scores, merged metadata, PCA plots (before & after ComBat)
- `results/xCell_results/Corrected_scores/` — 3 ComBat-corrected enrichment score matrices

---

### Step 3b — Subsetting to Common Genes (`03b_subset_common_genes.R`)

**Purpose:** Subset each batch-corrected expression matrix to retain only genes that appear in the corresponding metaDEG list, ensuring a consistent gene universe across all downstream analyses.

**Method:**
1. For each condition (`cHP_Ctrl`, `cHP_IPF`, `IPF_Ctrl`):
   - Read the `*_Common_genes_same_direction_expr.csv` to get the gene list.
   - Read the batch-corrected expression matrix from `results/Corrected_matrix/`.
   - Compute the intersection of genes.
   - Subset and export the reduced matrix.

**Key Outputs:**
- `results/Subsetted_corrected_matrix/` — 3 subsetted expression matrices

---

### Step 4 — Correlation Heatmaps (`04_correlation_heatmaps.R`)

**Purpose:** Visualise inter-sample and inter-gene correlations in the batch-corrected expression data to assess data quality and identify co-expression patterns.

**Method:**
1. **Sample-sample correlation** — Pearson correlation across all genes for each pair of samples.
   - Annotated by condition (Control, HP, IPF).
   - Hierarchical clustering using Ward's method (`ward.D2`).
   - Diverging blue–white–red colour palette.
2. **Gene-gene (feature) correlation** — For the top 100 most variable genes (selected by row-wise variance).
   - Pearson correlation across all samples.
   - Hierarchical clustering.

**Key Outputs:**
- `results/Correlation_Heatmaps/` — 12 correlation heatmaps (sample + feature, 3 conditions, PNG + PDF)

---

### Step 5 — mRMR Feature Selection (`05_mRMR.ipynb`)

**Purpose:** Apply Minimum Redundancy Maximum Relevance (mRMR) feature selection to identify a ranked list of genes that are highly informative for classification while minimising inter-gene redundancy.

**Method (Python / Jupyter Notebook):**
1. **Load** the subsetted batch-corrected expression matrices and metadata for each condition.
2. **Prepare** feature matrices (samples × genes) with binary condition labels.
3. **Split** data into train (70%) and test (30%) sets — stratified by condition.
4. **Run mRMR** feature selection on training data:
   - mRMR selects features that have high mutual information with the target variable (relevance) while having low mutual information with already-selected features (minimum redundancy).
5. **Rank genes** by cumulative mRMR score.
6. **Visualise** cumulative feature importance scores vs. number of genes.
7. **Export** ranked gene lists and train/test data splits.

**Key Outputs:**
- `train_data/` — 3 training CSVs + 3 mRMR gene lists

| Condition | Total mRMR Genes | Top 5 Genes |
|:----------|:-----------------|:------------|
| cHP vs Control | ~221 | ADRA1A, BTLA, BTNL9, PTPRB, AFF3 |
| cHP vs IPF | ~157 | RNF208, CALHM6, NKX6-2, PRKY, COL5A3 |
| IPF vs Control | ~221 | KLRG2, ADRA1A, ITLN2, BAAT, EPB41L5 |

- `test_data/` — 3 test CSVs
- `results/mRMR_plots/` — Cumulative feature score plots

---

### Step 6 — Ridge Regression & GLM Modelling (`06_lasso_refactored.R`)

**Purpose:** Build and evaluate classification models using the mRMR-selected features to identify the optimal minimal gene signature for each condition.

**Workflow (per condition):**

#### A. Incremental Gene Sweep (Ridge Regression)
1. Starting from the top 2 mRMR genes up to all N, fit a **Ridge regression** model (`alpha = 0`, 10-fold CV) at each step.
2. Ridge regression was chosen over LASSO because:
   - mRMR has already handled feature selection; Ridge keeps all features in the model.
   - Ridge handles multicollinearity more gracefully via coefficient shrinkage without variable elimination.
3. For each `k` (number of genes), compute **AUC** and **Accuracy** on the held-out test set.
4. Plot AUC & Accuracy vs. Number of Genes.

#### B. Interactive Optimal Gene Selection
- The pipeline **pauses** for the user to review the sweep plot and choose `optimal_k` — the number of top mRMR genes to carry forward.
- This human-in-the-loop step ensures domain expertise guides final gene panel selection.

**Selected gene panels:**

| Condition | Optimal k | Selected Genes |
|:----------|:----------|:---------------|
| cHP vs Control | 5 | ADRA1A, BTLA, BTNL9, PTPRB, AFF3 |
| cHP vs IPF | 3 | RNF208, CALHM6, NKX6-2 |
| IPF vs Control | 6 | KLRG2, ADRA1A, ITLN2, BAAT, EPB41L5, FHL2 |

#### C. Per-Gene ROC Analysis
- Individual ROC curves (with AUC + 95% CI) for each selected gene on both training and test sets.
- Panel and individual PNG plots generated.

#### D. Combinatorial Gene Search (GLM)
- Evaluate all 2–5 gene combinations from the selected panel using **unpenalised logistic regression (GLM)**.
- Rationale: mRMR genes are already filtered for independence; no further regularisation needed for an accurate assessment of joint predictive power.
- AUC and Accuracy computed for each combination on the test set.
- Top 10 combinations plotted with individual ROC curves (train + test).

**Feature Scaling:**
- Z-score standardisation of training features.
- Test features scaled using training mean and standard deviation (no data leakage).

**Key Outputs:**
- `results/lasso/<condition>/` — Per-condition outputs:
  - `gene_sweep_performance.csv` — AUC/Accuracy for each k
  - `LASSO_AUC_and_Accuracy_plot.png` — Sweep visualisation
  - `lasso_selected_genes.txt` — Final gene panel
  - `per_gene_AUC_train.csv` / `per_gene_AUC_test.csv` — Individual gene AUCs
  - `Predictor_genes_roc_train_panel.png` / `...test_panel.png` — ROC panel plots
  - `GLM_combinations_AUC_results.csv` — All combination results
  - `Combination_Top10_Test_ROC/` / `Combination_Top10_Train_ROC/` — Top-10 combination ROC plots

---

### Step 7 — Boxplots & Statistical Testing (`07_boxplot.R`)

**Purpose:** Generate publication-quality comparison boxplots for each selected gene and compute Mann-Whitney U test statistics.

**Method:**
1. For each selected gene in each condition:
   - Create a **boxplot with jittered data points** for both training and test sets.
   - Perform **Mann-Whitney U test (Wilcoxon rank-sum)** via `ggpubr::stat_compare_means`.
   - Significance brackets annotated directly on the plots.
2. **Per-gene statistics CSV** generated, containing:
   - Group means (train/test × positive/negative)
   - Raw p-values (Wilcoxon)
   - BH-adjusted p-values
   - log₂ Fold Change
   - Merged columns from the original metaDEG files (number of datasets, mean absolute logFC, per-dataset logFC, meta FDR)

**Plot Aesthetics:**
- Clean `theme_classic` style with semi-transparent boxplots, jittered individual points.
- Blue/coral colour scheme for negative/positive groups.
- 600 DPI PNG + PDF output.

**Key Outputs:**
- `results/boxplots/<condition>/train/` — Training set boxplots (PNG + PDF per gene)
- `results/boxplots/<condition>/test/` — Test set boxplots (PNG + PDF per gene)
- `results/boxplots/<condition>/<condition>_gene_statistics.csv` — Statistical summary

---

### Step 8 — Independent Validation AUC-ROC (`08_independent_AUC_ROC.R`)

**Purpose:** Validate the discovered biomarker genes and their combinations on three **independent** raw-count datasets that were not used during model training.

**Validation Datasets:**

| Dataset | Source | Description |
|:--------|:-------|:------------|
| GSE184316 (full) | Full condition raw counts | Independent sequencing platform |
| GSE150910 (Biopsy) | Biopsy-only subset | Tissue-specific subset |
| GSE150910 (Explant) | Explant-only subset | Tissue-specific subset |

**Method:**
1. For **single genes**: Compute AUC directly from raw expression values vs. binary labels using `pROC::roc`.
2. For **gene combinations**: Fit a logistic regression (GLM) on the validation dataset and compute AUC from predicted probabilities.
3. **No training data leakage** — models are fit and evaluated on each validation set independently.

**log₂FC Heatmaps:**
- For each condition, a heatmap of per-dataset log₂FC values is generated (genes × datasets).
- Blue → White → Red diverging colour scale.
- Column labels repositioned above the heatmap for clean visualisation.

**Key Outputs:**
- `results/independent_AUC_ROC/<condition>/`
  - `<condition>_single_gene_AUC.csv` — AUC for each gene across 3 validation datasets
  - `<condition>_combination_AUC.csv` — AUC for each gene combination
- `results/boxplots/<condition>/<condition>_logFC_heatmap.png` — Cross-dataset log₂FC heatmaps

---

### Pairwise xCell Comparison (`pairwise_comparison.R`)

**Purpose:** Perform pairwise statistical comparison of immune cell type enrichment scores between conditions and generate a combined bubble heatmap summarising significant immune cell differences.

**Method:**
1. **For each comparison** (HP vs Ctrl, HP vs IPF, IPF vs Ctrl):
   - Load merged xCell enrichment scores and metadata.
   - For each cell type: **Mann-Whitney U test** (Wilcoxon rank-sum, two-sided).
   - Compute fold change, BH-corrected FDR.
   - Generate individual boxplots for all cell types in a multi-page PDF.
2. **Combined Bubble Heatmap:**
   - X-axis: Significant cell types (ordered by number of comparisons they are significant in).
   - Y-axis: Comparisons (HP vs Ctrl, HP vs IPF, IPF vs Ctrl).
   - Bubble size: FDR significance category (< 0.0001, < 0.001, < 0.01, < 0.05).
   - Bubble colour: log₂ Fold Change (blue = downregulated, red = upregulated).
   - Non-significant cell types and entries are excluded.
3. **Repeated for ComBat-corrected enrichment scores** — same pipeline applied to batch-corrected xCell scores from Step 3.

**Key Outputs:**
- `results/xCell_Stats_Plots/` — Statistics CSVs, per-cell-type boxplot PDFs, bubble heatmap (raw scores)
- `results/xCell_Stats_Plots_Corrected/` — Same outputs for batch-corrected scores

---

## 7. Results Summary

### Discovered Gene Signatures

| Comparison | # Genes | Gene Panel | Key Finding |
|:-----------|:--------|:-----------|:------------|
| **cHP vs Control** | 5 | ADRA1A, BTLA, BTNL9, PTPRB, AFF3 | Immune receptor and adhesion molecules distinguish cHP from controls |
| **cHP vs IPF** | 3 | RNF208, CALHM6, NKX6-2 | Minimal 3-gene panel differentiates these clinically similar diseases |
| **IPF vs Control** | 6 | KLRG2, ADRA1A, ITLN2, BAAT, EPB41L5, FHL2 | Diverse functional set including metabolic (BAAT) and immune (KLRG2) genes |

### Model Performance

Results from the Ridge Regression + GLM pipeline (specific AUC values available in the per-condition CSV files under `results/lasso/`):

- **Gene sweep plots** show optimal performance plateaus, justifying the selected `optimal_k` values.
- **Combinatorial GLM** results identify top gene combinations with high AUC on both training and test sets.
- **Independent validation** on GSE184316, GSE150910-Biopsy, and GSE150910-Explant confirms cross-platform and cross-tissue generalisability.

### Immune Microenvironment

The xCell analysis and bubble heatmaps reveal:
- Differentially enriched immune cell types across conditions (detailed in `results/xCell_Stats_Plots/Stats_*.csv`).
- Batch correction of enrichment scores confirms findings are robust to platform effects.

---

## 8. Directory Structure

```
meta_analysis/
│
├── README.md                          ← This file
├── xCell.Rproj                        ← RStudio project file
│
├── 01_batch_correction.R              ← Step 1: Merge + ComBat batch correction
├── 02_tpm.R                           ← Step 2: TPM normalisation for xCell
├── 03_xCell.R                         ← Step 3: xCell score merging + batch correction
├── 03b_subset_common_genes.R          ← Step 3b: Subset to metaDEG genes
├── 04_correlation_heatmaps.R          ← Step 4: Sample & gene correlation heatmaps
├── 05_mRMR.ipynb                      ← Step 5: mRMR feature selection (Python)
├── 06_lasso_refactored.R              ← Step 6: Ridge + GLM modelling pipeline
├── 07_boxplot.R                       ← Step 7: Gene boxplots + statistics
├── 08_independent_AUC_ROC.R           ← Step 8: Independent validation
├── pairwise_comparison.R              ← xCell pairwise statistics + bubble heatmap
├── xCell_analysis.Rmd                 ← Exploratory xCell t-test notebook
│
├── cHP_Ctrl_Common_genes_same_direction_expr.csv   ← MetaDEGs: cHP vs Control
├── cHP_IPF_Common_genes_same_direction_expr.csv    ← MetaDEGs: cHP vs IPF
├── IPF_Ctrl_Common_genes_same_direction_expr.csv   ← MetaDEGs: IPF vs Control
├── Signatures.xlsx                                 ← Gene signature reference
│
├── data/
│   ├── GSE184316/                     ← Raw counts + sample info (GSE184316)
│   └── GSE150910/                     ← Raw counts + sample info (GSE150910)
│                                        Including biopsy/explant metadata
│
├── train_data/                        ← Training splits + mRMR gene lists
│   ├── cHP_Ctrl_train_data.csv
│   ├── cHP_Ctrl_top_mrmr_genes.txt    (221 genes)
│   ├── cHP_IPF_train_data.csv
│   ├── cHP_IPF_top_mrmr_genes.txt     (157 genes)
│   ├── IPF_Ctrl_train_data.csv
│   └── IPF_Ctrl_top_mrmr_genes.txt    (221 genes)
│
├── test_data/                         ← Test splits
│   ├── cHP_Ctrl_test_data.csv
│   ├── cHP_IPF_test_data.csv
│   └── IPF_Ctrl_test_data.csv
│
└── results/
    ├── PCA_plots/                     ← PCA before/after ComBat (expression data)
    ├── Raw_merged_data/               ← Raw merged count matrices
    ├── Corrected_matrix/              ← Batch-corrected expression matrices
    ├── Subsetted_corrected_matrix/    ← Subsetted to metaDEG genes only
    ├── Correlation_Heatmaps/          ← Sample & gene correlation heatmaps
    ├── mRMR_plots/                    ← Cumulative mRMR feature importance plots
    ├── lasso/                         ← Ridge/GLM modelling results
    │   ├── cHP_Ctrl/                  ← Gene sweep, ROCs, combination results
    │   ├── cHP_IPF/
    │   └── IPF_Ctrl/
    ├── boxplots/                      ← Gene expression boxplots + statistics
    │   ├── cHP_Ctrl/
    │   ├── cHP_IPF/
    │   └── IPF_Ctrl/
    ├── independent_AUC_ROC/           ← Independent validation AUC results
    │   ├── cHP_Ctrl/
    │   ├── cHP_IPF/
    │   └── IPF_Ctrl/
    ├── xCell_input/                   ← TPM matrices for xCell upload
    ├── xCell_results/                 ← xCell scores, PCA, corrected scores
    ├── xCell_Stats_Plots/             ← Pairwise xCell stats + bubble heatmap
    └── xCell_Stats_Plots_Corrected/   ← Same for batch-corrected xCell scores
```

---

## 9. How to Reproduce

### Prerequisites

1. **R** (≥ 4.2) with the following packages installed:
   - CRAN: `tidyverse`, `ggplot2`, `ggpubr`, `glmnet`, `pROC`, `reshape2`, `pheatmap`, `RColorBrewer`, `gridExtra`, `scales`
   - Bioconductor: `DESeq2`, `sva`, `biomaRt`
2. **Python** (≥ 3.8) with: `pandas`, `numpy`, `scikit-learn`, `mrmr-selection` (or equivalent mRMR library)
3. **xCell** web tool access: [https://xcell.ucsf.edu/](https://xcell.ucsf.edu/)

### Step-by-Step Reproduction

```bash
# Clone or download the project
# Ensure all data files are present in ./data/

# Open the project in RStudio (xCell.Rproj)
```

#### Step 1: Batch Correction
```r
source("01_batch_correction.R")
# Outputs: results/PCA_plots/, results/Corrected_matrix/, results/Raw_merged_data/
```

#### Step 2: TPM Normalisation
```r
source("02_tpm.R")
# Outputs: results/xCell_input/ (6 TPM files)
# NOTE: Requires internet connection for Ensembl biomaRt queries
```

#### Step 3a: xCell Deconvolution (External)
1. Upload each TPM file from `results/xCell_input/` to [xCell web portal](https://xcell.ucsf.edu/).
2. Download results and place in appropriate subdirectories under `results/xCell_results/`.

#### Step 3b: xCell Score Analysis
```r
source("03_xCell.R")
# Outputs: results/xCell_results/ (merged scores, PCA, corrected scores)
```

#### Step 3c: Subset Expression Matrices
```r
source("03b_subset_common_genes.R")
# Outputs: results/Subsetted_corrected_matrix/
```

#### Step 4: Correlation Heatmaps
```r
source("04_correlation_heatmaps.R")
# Outputs: results/Correlation_Heatmaps/
```

#### Step 5: mRMR Feature Selection
```bash
# Open and run 05_mRMR.ipynb in Jupyter Notebook / JupyterLab
# This produces:
#   train_data/*_train_data.csv
#   train_data/*_top_mrmr_genes.txt
#   test_data/*_test_data.csv
#   results/mRMR_plots/
```

#### Step 6: Ridge + GLM Modelling
```r
source("06_lasso_refactored.R")
# INTERACTIVE: The script will pause for each condition and prompt you
#   to enter optimal_k (the number of top mRMR genes to carry forward).
# Review the AUC/Accuracy sweep plot and gene_sweep_performance.csv
#   before entering your choice.
#
# Recommended values based on original analysis:
#   cHP_Ctrl: optimal_k = 5
#   cHP_IPF:  optimal_k = 3
#   IPF_Ctrl: optimal_k = 6
#
# Outputs: results/lasso/
```

#### Step 7: Boxplots
```r
source("07_boxplot.R")
# Outputs: results/boxplots/
```

#### Step 8: Independent Validation
```r
source("08_independent_AUC_ROC.R")
# Outputs: results/independent_AUC_ROC/, logFC heatmaps in results/boxplots/
```

#### Pairwise xCell Statistics
```r
source("pairwise_comparison.R")
# Outputs: results/xCell_Stats_Plots/, results/xCell_Stats_Plots_Corrected/
```

### Important Notes

- **Script order matters**: Scripts are numbered sequentially and depend on outputs from earlier steps.
- **Internet required**: Step 2 (`02_tpm.R`) requires internet access for Ensembl queries via `biomaRt`.
- **xCell is external**: Step 3 requires manual upload/download from the xCell web portal.
- **Interactive step**: Step 6 (`06_lasso_refactored.R`) requires user input to select `optimal_k` for each condition.
- **Reproducibility**: Random seeds are set (`set.seed(123)`) for the Ridge regression cross-validation.
- **Working directory**: All scripts assume the project root (`meta_analysis/`) as the working directory.

---

## 10. Software & Dependencies

### R Packages

| Package | Purpose | Source |
|:--------|:--------|:-------|
| `tidyverse` | Data wrangling (dplyr, tidyr, readr, etc.) | CRAN |
| `DESeq2` | DESeqDataSet creation, VST transformation | Bioconductor |
| `sva` | ComBat batch correction | Bioconductor |
| `biomaRt` | Gene annotation & transcript lengths from Ensembl | Bioconductor |
| `ggplot2` | Visualisation (PCA, boxplots, bubble plots) | CRAN |
| `ggpubr` | Statistical comparison brackets on plots | CRAN |
| `glmnet` | Ridge/LASSO regression (cv.glmnet) | CRAN |
| `pROC` | ROC curves, AUC, 95% CI calculation | CRAN |
| `pheatmap` | Heatmap generation (correlation, logFC) | CRAN |
| `RColorBrewer` | Colour palettes | CRAN |
| `gridExtra` | Multi-panel plot arrangement | CRAN |
| `reshape2` | Data reshaping | CRAN |
| `scales` | Axis scaling for ggplot2 | CRAN |

### Python Packages

| Package | Purpose |
|:--------|:--------|
| `pandas` | Data manipulation |
| `numpy` | Numerical operations |
| `scikit-learn` | Train/test splitting, preprocessing |
| `mrmr-selection` (or equivalent) | mRMR feature selection algorithm |

### External Tools

| Tool | Purpose | URL |
|:-----|:--------|:----|
| xCell | Cell type enrichment analysis from gene expression | [https://xcell.ucsf.edu/](https://xcell.ucsf.edu/) |

---

## 11. Key Output Files

### Statistical Results

| File | Content |
|:-----|:--------|
| `results/lasso/*/gene_sweep_performance.csv` | AUC & Accuracy for each gene count (k=2..N) |
| `results/lasso/*/per_gene_AUC_train.csv` | Per-gene AUC on training set |
| `results/lasso/*/per_gene_AUC_test.csv` | Per-gene AUC on test set |
| `results/lasso/*/GLM_combinations_AUC_results.csv` | All gene combination AUCs |
| `results/boxplots/*/<condition>_gene_statistics.csv` | Mean expression, p-values, padj, log₂FC per gene |
| `results/independent_AUC_ROC/*/<condition>_single_gene_AUC.csv` | Independent validation AUC per gene |
| `results/independent_AUC_ROC/*/<condition>_combination_AUC.csv` | Independent validation AUC per combination |
| `results/xCell_Stats_Plots/Stats_*.csv` | xCell cell-type statistics (p-values, FDR, fold change) |

### Key Figures

| File | Content |
|:-----|:--------|
| `results/PCA_plots/Merged_*_PCA_*.png` | PCA before/after batch correction |
| `results/Correlation_Heatmaps/*.png` | Sample and gene correlation heatmaps |
| `results/lasso/*/LASSO_AUC_and_Accuracy_plot.png` | Gene sweep performance curve |
| `results/lasso/*/Predictor_genes_roc_*_panel.png` | Multi-gene ROC panel plots |
| `results/boxplots/*/*/<Gene>_*.png` | Gene expression boxplots with significance |
| `results/boxplots/*/<condition>_logFC_heatmap.png` | Cross-dataset log₂FC heatmaps |
| `results/xCell_Stats_Plots/BubbleHeatmap_AllComparisons.png` | Combined immune cell bubble heatmap |

---

## 12. Author

**Abhijit Saha**

---

## 13. License

This project is for academic research purposes. Please contact the author for usage permissions.
