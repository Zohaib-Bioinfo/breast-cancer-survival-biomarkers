# ML-Derived Biomarkers for Breast Cancer Survival Prognosis

Independent prognostic validation of machine learning–derived transcriptomic biomarkers in the METABRIC cohort, extending a GSE45827-trained XGBoost classifier to clinical survival outcomes.

## Overview

Most published breast cancer subtype classifiers are evaluated purely on classification accuracy, with no test of whether the genes they identify actually predict clinical outcomes. This project addresses that gap by:

1. Taking the 46 SHAP-selected biomarker genes from a GSE45827-trained XGBoost subtype classifier (see [breast-cancer-subtype-classification](https://github.com/Zohaib-Bioinfo/breast-cancer-subtype-classification))
2. Validating their **independent prognostic value** in the METABRIC cohort (n=1,608) using multivariate Cox regression — adjusting for age, lymph node status, Nottingham Prognostic Index (NPI), and PAM50 molecular subtype

## Key Findings

- **CDCA5** (HR = 1.144, p = 0.023): Higher expression independently associated with worse overall survival
- **CMC2** (HR = 0.919, p = 0.049): Higher expression independently associated with better overall survival
- Adding 46 ML-derived biomarker genes to clinical variables improved concordance index from **0.657 → 0.679 (+0.022)**
- CDCA5 high vs low expressers: Log-rank p = 1.69e-09 with sustained survival separation over 300 months

## Dataset

| Dataset | Source | Samples | Use |
|---|---|---|---|
| GSE45827 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45827) | 130 tumor samples | XGBoost classifier training + SHAP gene selection |
| METABRIC | [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) | 1,608 patients | Independent survival validation cohort |

## Methods

### Phase 1 — Biomarker Gene Discovery (GSE45827)
- Affymetrix GPL570 microarray data retrieved via `GEOparse`
- Probe IDs mapped to HGNC symbols via `mygene`
- 300 features selected by mutual information (`SelectKBest`)
- XGBoost multi-class classifier trained on 4 subtypes: Basal, HER2, Luminal A, Luminal B
- 5-fold stratified cross-validation accuracy: **94.6% ± 3.1%**
- SHAP TreeExplainer used to rank feature importance → top 46 genes retained

### Phase 2 — METABRIC Data Acquisition
- Clinical survival data (OS, RFS) downloaded from cBioPortal
- Gene expression z-scores for 46 biomarker genes retrieved via cBioPortal REST API
- 1,608 patients retained after filtering to 4 PAM50 subtypes and removing missing OS data

### Phase 3 — Survival Analysis
- **Kaplan-Meier** curves by molecular subtype (multivariate log-rank test)
- **Multivariate Cox Proportional Hazards** regression (`lifelines`, penalizer=0.1):
  - Model 1 (baseline): Age + lymph nodes + NPI + subtype
  - Model 2 (novel): Model 1 + 46 ML biomarker genes
- **Concordance Index (C-index)** comparison between models
- **CDCA5 stratified KM analysis**: median split, log-rank test

## Results

### Subtype Survival Separation
Kaplan-Meier curves confirm expected PAM50 survival hierarchy in METABRIC
(Log-rank p = 1.70e-09). Basal subtype shows characteristic crossing hazard pattern — poor early survival with relative stabilization beyond 250 months.

### Independently Prognostic Genes

| Gene | Hazard Ratio | 95% CI | p-value | Direction |
|---|---|---|---|---|
| CDCA5 | 1.144 | 1.019–1.284 | 0.023 | Higher expression → worse survival |
| CMC2 | 0.919 | 0.845–1.000 | 0.049 | Higher expression → better survival |

Both significant after adjusting for age, lymph nodes, NPI, and PAM50 subtype.

### Model Comparison

| Model | C-index |
|---|---|
| Clinical only | 0.6567 |
| Clinical + ML genes | 0.6788 |
| **Improvement** | **+0.0221** |

## Figures

| Figure | Description |
|---|---|
| `figures/KM_survival_by_subtype.png` | Kaplan-Meier overall survival by PAM50 subtype |
| `figures/forest_plot_biomarkers.png` | Forest plot — Cox hazard ratios for all 46 biomarker genes |
| `figures/KM_CDCA5_expression.png` | Survival by CDCA5 expression level (high vs low) |
| `figures/cindex_comparison.png` | C-index improvement: clinical vs clinical + ML genes |

## Repository Structure

```
breast-cancer-survival-biomarkers/
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
├── notebooks/
│   └── PrognosticValidationMLBiomarkers.ipynb
├── figures/
│   ├── KM_survival_by_subtype.png
│   ├── forest_plot_biomarkers.png
│   ├── KM_CDCA5_expression.png
│   └── cindex_comparison.png
└── results/
    └── metabric_survival_analysis.csv
```

## Reproduction

```bash
git clone https://github.com/Zohaib-Bioinfo/breast-cancer-survival-biomarkers.git
cd breast-cancer-survival-biomarkers
pip install -r requirements.txt
```

Open `notebooks/PrognosticValidationMLBiomarkers.ipynb` in Google Colab or Jupyter.

**Data note:** Clinical data (`data_clinical_patient.txt`) must be downloaded manually from [cBioPortal METABRIC](https://www.cbioportal.org/study/summary?id=brca_metabric) due to file size. Gene expression values are retrieved automatically via the cBioPortal REST API in the notebook.

## Related Project

This project extends: [breast-cancer-subtype-classification](https://github.com/Zohaib-Bioinfo/breast-cancer-subtype-classification) — GSE45827 XGBoost classifier (94.6% CV accuracy) that generated the biomarker gene panel used here.

## Tools & Technologies

Python, lifelines, scikit-learn, XGBoost, SHAP, pandas, matplotlib, cBioPortal REST API, GEOparse, mygene

## Limitations

- Cox model uses L2 penalization (penalizer=0.1) to handle 46 simultaneous gene predictors — results should be considered exploratory
- Expression platforms differ between GSE45827 (Affymetrix microarray) and METABRIC (Illumina HT-12 microarray) — cross-platform normalization not applied
- No experimental (wet-lab) validation of CDCA5 or CMC2 prognostic roles
- Median split for CDCA5 KM analysis is a post-hoc visualization — not a pre-specified cutpoint

## Author

Muhammad Zohaib — BS Bioinformatics, University of Agriculture Faisalabad (UAF)

## License

MIT License — see `LICENSE`
