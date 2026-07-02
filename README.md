# Sex-Specific Subphenotyping of HFpEF

Multi-omic (proteomics + clinical phenomics) patient stratification in heart
failure with preserved ejection fraction (HFpEF), analyzed separately for
women and men, with external validation in an independent cohort.

## Pipeline

1. **Data integration** — proteomic and clinical/phenomic feature panels are
   imputed, filtered (near-zero-variance and highly correlated features
   removed), and combined via **Similarity Network Fusion (SNF)**.
2. **Clustering** — spectral clustering on the fused network defines
   candidate patient subgroups (2–10 clusters evaluated).
3. **Survival validation** — Cox proportional hazards models and
   Kaplan-Meier curves test whether the discovered clusters differ in
   worsening-heart-failure risk.
4. **Baseline characterization** — clinical/lab tables comparing clusters.
5. **Out-of-cohort validation** — an XGBoost classifier is trained on the
   discovery-cohort cluster labels and used to project cluster membership
   onto an independent validation cohort, whose survival outcomes are then
   checked against the same clustering.
6. **Transcriptomics** — RNA-seq counts for the discovered clusters are
   tested for differential expression, followed by KEGG/Reactome gene set
   enrichment analysis (GSEA).

Steps 1–5 are run once for the female subgroup and once for the male
subgroup.

## Data availability

This repository contains **code only**. The analysis uses patient-level
clinical, proteomic, and RNA-seq data from the MYOVASC and Gutenberg Health
Study (GHS) cohorts, which are not publicly available due to participant
privacy. The data-loading calls near the top of the script
(`read_SQL_MyoVasc_BL`, `read_SQL_A6_data`, etc.) depend on internal
database access and institution-specific helper scripts and will not run
outside that environment.

This script is shared to document the analysis methodology for
reproducibility/review purposes, not as a runnable end-to-end pipeline.

## Configuration

All machine- and institution-specific paths have been factored out into
environment variables at the top of the script, rather than hardcoded:

| Variable | Purpose |
|---|---|
| `PROJECT_DIR` | Working directory for the project |
| `MYOVASC_SQL_DIR` | Location of MYOVASC data-loading scripts |
| `GHS_SQL_DIR` | Location of GHS data-loading scripts |
| `SHARED_FUNCTIONS_DIR` | Location of shared helper functions (e.g. `baseline13.r`) |
| `PROJECT_DATA_DIR` | Location of local data extracts (defaults to `<PROJECT_DIR>/data`) |
| `PROJECT_OUTPUT_DIR` | Where results (RDS/RData) are written (defaults to `<PROJECT_DIR>/output`) |
| `DB_USER` | Database username for internal SQL-loading functions |

Set these via a local `.Renviron` file (which you should **not** commit)
before running, e.g.:

```
PROJECT_DIR=/path/to/project
MYOVASC_SQL_DIR=/path/to/myovasc/r_sql
GHS_SQL_DIR=/path/to/ghs/r_sql
SHARED_FUNCTIONS_DIR=/path/to/shared/r_functions
DB_USER=your_username
```

## Requirements

R packages: `SNFtool`, `survival`, `survminer`, `survcomp`, `flexclust`,
`ggplot2`, `Rtsne`, `rstatix`, `stringr`, `janitor`, `ggsurvfit`, `readxl`,
`caret`, `Hmisc`, `xgboost`, `pROC`, `PRROC`, `dplyr`, `clusterProfiler`,
`org.Hs.eg.db`, `ReactomePA`, `enrichplot`.
