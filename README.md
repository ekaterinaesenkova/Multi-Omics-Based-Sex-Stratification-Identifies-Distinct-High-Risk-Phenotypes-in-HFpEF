# Sex-stratified network integration of inflammatory proteins and clinical data in HFpEF

Analysis code for:

> Esenkova EE, Koeck T, Rapp S, Bauer KI, Zeid S, Wild PS, Casiraghi E, Araldi E.
> **A shared inflammatory high-risk phenotype in HFpEF is linked to diabetes and insulin resistance in women and to renal dysfunction in men: sex-stratified network integration with external validation.**
> *Cardiovascular Diabetology* (submitted, 2026).

## Overview

The pipeline integrates 92 circulating inflammatory proteins (Olink Target 96 Inflammation) and 49 clinical variables in patients with heart failure with preserved ejection fraction (HFpEF). The analysis is done **separately in women and men**. It:

1. builds sex-specific patient similarity networks with **Similarity Network Fusion (SNF)** and identifies clusters by **spectral clustering**;
2. relates the clusters to **worsening heart failure (WHF)**, diabetes, insulin resistance (HOMA-IR) and glucose-lowering therapy;
3. tests sex differences formally (cluster × sex, diabetes × sex and protein × sex interactions, plus a pooled-SNF comparator);
4. validates the cluster structure in an independent cohort by **classifier transfer (XGBoost)** and **de novo clustering**;
5. characterises cluster-discriminating features (SHAP values, protein comparisons) and whole-blood transcriptomic correlates (Reactome GSEA).

| | Discovery | Validation |
|---|---|---|
| Cohort | MyoVasc (NCT04064450) | Gutenberg Health Study (GHS) |
| HFpEF participants with proteomics | 698 (319 women, 379 men) | 246 (134 women, 112 men) |

## Data availability

**No individual-level data are included in this repository.** MyoVasc and GHS data are under controlled access according to the participants' informed consent. Requests for data access can be addressed to the principal investigators of the studies (see the *Availability of data and materials* section of the article).

The data-loading functions at the top of the script (`read_SQL_MyoVasc_BL()`, `read_SQL_A6_data()`, `baseline13()`) are internal to the University Medical Center Mainz data infrastructure and are not part of this repository. To run the pipeline on other data, replace Block 0 with your own import that returns one row per participant, keyed by a participant ID column (`v11_sid01`), with the variables listed in `clin_vars` and `prot_vars`.

## Repository contents

```
Sex-specific_HFpEF.R   complete analysis pipeline (analyses, tables, figures)
README.md              this file
```

## Requirements

Analyses were run in **R 4.2.1**. CRAN packages:

```r
install.packages(c("SNFtool", "survival", "survminer", "caret", "missRanger", "xgboost",
                   "dplyr", "stringr", "janitor", "readxl", "Hmisc",
                   "ggplot2", "patchwork", "ggrepel"))
```

Bioconductor packages for the transcriptomic analysis:

```r
BiocManager::install(c("DESeq2", "ReactomePA", "enrichplot", "org.Hs.eg.db"))
```

Random seeds are fixed (`set.seed(256)`; the classifier cross-validation repetition *r* uses `set.seed(r)`).

## Pipeline

Run the script from top to bottom. Each block saves its results (`.rds` / `.csv`) to `project_dir`, so later blocks can be re-run on their own.

| Block | Content | Manuscript | Output |
|---|---|---|---|
| 0 | Data import, protein harmonisation, derived variables (LV mass index, TyG index) | Methods | – |
| 1 | **Discovery:** HFpEF selection; per sex, missing-data filtering (>30% per variable, >50% per participant), random-forest imputation (missRanger), near-zero-variance and correlation filtering (\|r\| ≥ 0.9), z-standardisation, SNF (K = 20, σ = 0.5, t = 20) and spectral clustering (number of clusters by eigengap and rotation cost, range 2–6); Cox models for WHF | Results 1, Fig. 1, Table 1 | `WP1_discovery.rds` |
| 2 | Cardiac characteristics by cluster; stability (200 × 80% subsamples); SNF hyperparameter grid (K = 10–30, σ = 0.3–0.8) | Results 1 | `WP2_doc_stability.rds`, `WP2_cardiac_by_cluster_MyoVasc.csv` |
| 3 / 3b | **GHS validation:** Olink harmonisation check; XGBoost classifier (high-risk = positive class; 20 × stratified 10-fold CV); transfer to GHS after within-cohort imputation and scaling; de novo SNF in GHS matched to MyoVasc clusters by profile correlation (no outcome information); Cox models | Results 3, Fig. 2, Suppl. Table 4 | `WP4_GHS_validation.rds` |
| 4 / 4b | Diabetes, HOMA-IR, TyG, glucose-lowering therapy; one-at-a-time sensitivity adjustment (eGFR, diabetes, BMI, atrial fibrillation, RAS inhibitors, glucose-lowering therapy, HOMA-IR); clustering without glycaemic markers; fixed-effect meta-analysis of both cohorts | Results 2 | `WP5_diabetes.rds` |
| 5 | Cluster × sex and diabetes × sex interactions (ratio of HRs); pooled-SNF comparator; protein × sex interactions with Benjamini–Hochberg correction; protein concordance between sexes | Results 4, Fig. 7 | `WP6_protein_sex_interactions.csv` |
| 7 | Reproductive and hormonal history in women (menarche, pregnancies, menopause, oral contraceptives, hormone replacement therapy), including item non-response | Results 2, Suppl. Table 5 | `SuppTable5_reproductive_hormonal_factors.csv` |
| Tables & figures | Baseline tables; cumulative-incidence plots; SHAP importance; protein box plots; Reactome GSEA dot plots; protein concordance | Table 1, Figs. 1–7, Suppl. Tables 1–5 | `Figure*.png/.tiff`, `SuppTable*.csv` |

### Transcriptomic analysis

Whole-blood 3′ mRNA-seq data (QuantSeq) were variance-stabilised (DESeq2). For each gene, a linear model `expression ~ cluster + lymphocytes + neutrophils + monocytes` was fitted, separately by sex. Genes were ranked by the cluster t-statistic, and preranked GSEA was run against Reactome (`gsePathway`; gene-set size 10–500; Benjamini–Hochberg-adjusted p < 0.05). The cluster coefficient is a difference on the variance-stabilised scale, not a conventional log2 fold change.

## Revision history

- **2026 (resubmission):** revised pipeline. Changes:
  - cluster labels are matched to participants by ID throughout;
  - random-forest imputation (missRanger) within sex and data layer;
  - feature filtering before network construction;
  - stability and hyperparameter analyses;
  - de novo replication in GHS;
  - formal sex-interaction analyses and pooled comparator;
  - diabetes and insulin-resistance analyses;
  - reproductive and hormonal factors;
  - corrected RNA-seq effect-size reporting.

  All results in the article are based on this version.
- **2025:** initial version (superseded).

## Citation

If you use this code, please cite the article above (citation details will be added on publication).

## Contact

Ekaterina E. Esenkova, Preventive Cardiology, University Medical Center of the Johannes Gutenberg University Mainz, Germany.
Corresponding author: Elisa Araldi (araldiel@uni-mainz.de).


