# Breast Cancer Survival Biomarkers — METABRIC Prognostic Validation

Prognostic validation of an ML-derived transcriptomic biomarker panel in the METABRIC cohort (n=1,608), testing whether genes identified via SHAP-guided XGBoost classification carry independent prognostic value beyond PAM50 molecular subtype and established clinical variables.

> **Companion repository:** The biomarker panel validated here — along with subsequent pipeline robustness testing (classifier comparison, feature selection stability, external cohort classification validation, and extended pathway enrichment) — is documented in [breast-cancer-subtype-classification](https://github.com/Zohaib-Bioinfo/breast-cancer-subtype-classification). That repository establishes *how* the 46-gene panel was discovered and *why* it should be trusted; this repository establishes whether it carries genuine prognostic value.

## Key Findings

- **CDCA5** (Cell Division Cycle Associated 5 / Sororin): HR=1.144, 95% CI 1.019-1.284, p=0.023
- **CMC2** (COX Assembly Mitochondrial Protein 2): HR=0.919, 95% CI 0.845-1.000, p=0.049
- Both findings confirmed in a proportional hazards-corrected sensitivity analysis (CDCA5: HR=1.123, p=0.049; CMC2: HR=0.916, p=0.039)
- The 46-gene panel improved the concordance index from 0.657 (clinical-only) to 0.679
- Neither gene survived Benjamini-Hochberg FDR correction across all 46 genes tested; both are reported as exploratory candidates
- Companion pipeline robustness testing (see linked repository) shows CDCA5 supported by convergent evidence across four independent axes, while CMC2's supporting evidence is comparatively thinner — the two genes warrant differentiated confidence, not equal treatment

## Repository Structure

```
breast-cancer-survival-biomarkers/
├── README.md
├── notebooks/
│   └── 01_metabric_survival_validation.ipynb   # Cox regression, KM curves, PH testing, sensitivity analysis, C-index
├── data/
│   ├──  metabric_survival_analysis.csv          # Processed METABRIC cohort (1,608 patients x 73 columns)
|   └──  data_clinical_patient.txt
|   └──  data_clinical_sample.txt
├── results/                            # Cohort characteristics, Cox results, Schoenfeld PH test results
├── figures/                            # Figs 1-4 as referenced in the manuscript
                                  
```

## Data Availability

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

## Methods Summary

Clinical and expression data (46 biomarker gene z-scores, Illumina HT-12 v3) were retrieved from METABRIC via the cBioPortal REST API (study ID: `brca_metabric`). Analysis used the `lifelines` Python library: Kaplan-Meier estimation and multivariate log-rank testing for subtype-stratified survival, and L2-penalized multivariate Cox proportional hazards regression (two nested models: clinical-only baseline vs. clinical + 46-gene panel) adjusting for age at diagnosis, lymph node count, Nottingham Prognostic Index, and PAM50 subtype. The proportional hazards assumption was formally tested via Schoenfeld residuals (Grambsch-Therneau framework); seven covariates violated PH and were addressed in a pre-specified sensitivity analysis (stratification for clinical violators, exclusion for gene violators).

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

## Reproduction

```bash
git clone https://github.com/Zohaib-Bioinfo/breast-cancer-survival-biomarkers.git
cd breast-cancer-survival-biomarkers
pip install -r requirements.txt
```

Open `notebooks/PrognosticValidationMLBiomarkers.ipynb` in Google Colab or Jupyter.

**Data note:** Clinical data (`data_clinical_patient.txt`) must be downloaded manually from [cBioPortal METABRIC](https://www.cbioportal.org/study/summary?id=brca_metabric) due to file size. Gene expression values are retrieved automatically via the cBioPortal REST API in the notebook.

## Limitations

- Cox model uses L2 penalization (penalizer=0.1) to handle 46 simultaneous gene predictors — results should be considered exploratory
- Expression platforms differ between GSE45827 (Affymetrix microarray) and METABRIC (Illumina HT-12 microarray) — cross-platform normalization not applied
- No experimental (wet-lab) validation of CDCA5 or CMC2 prognostic roles
- Median split for CDCA5 KM analysis is a post-hoc visualization — not a pre-specified cutpoint

## Reproducing This Work

The notebook is self-contained and runs in Google Colab or any Python 3.12 environment with `lifelines`, `pandas`, `numpy`, `requests`, `matplotlib`, and `seaborn` installed. The processed METABRIC CSV in `data/` allows full reproduction without re-querying the cBioPortal API, though the notebook includes the original API retrieval code for transparency.

## Citation

If you use this analysis or its outputs, please cite the associated manuscript (full text and supporting information available in the [companion repository](https://github.com/Zohaib-Bioinfo/breast-cancer-subtype-classification)).

## Author

Muhammad Zohaib — BS Bioinformatics, University of Agriculture Faisalabad (UAF)

## License

MIT License — see `LICENSE`
