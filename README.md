# Breast Cancer Survival Biomarkers — METABRIC Prognostic Validation
### (Repo 2 of 2 — see `breast-cancer-subtype-classification` for Repo 1)

Prognostic validation of an ML-derived transcriptomic biomarker panel in the METABRIC cohort (n=1,608), testing whether genes identified via SHAP-guided XGBoost classification carry independent prognostic value beyond PAM50 molecular subtype and established clinical variables.

> **Companion repository:** The biomarker panel validated here — along with subsequent pipeline robustness testing (classifier comparison, feature selection stability, external cohort classification validation, and extended pathway enrichment) — is documented in [breast-cancer-subtype-classification](https://github.com/Zohaib-Bioinfo/breast-cancer-subtype-classification). That repository establishes *how* the 43-gene panel was discovered and *why* it should be trusted; this repository establishes whether it carries genuine prognostic value.

This repo answers one question: do the biomarker genes identified by the
GSE45827 subtype classifier (Repo 1) carry **independent prognostic
value** — i.e., do they predict overall survival in an entirely separate
cohort (METABRIC), beyond what's already explained by standard clinical
variables (age, lymph node status, NPI, molecular subtype)?

> ⚠️ **Reproducibility risk — this is the last of three manual
> copy-paste hand-offs.** `01_metabric_prognostic_validation.ipynb`
> Cell 5's `biomarker_genes` list is a **manually retyped copy** of
> `reannotate_probes_R.ipynb`'s output (`model_B_top50_reannotated.csv`)
> from Repo 1. There is no automated file read connecting the two repos.
> A typo here silently changes which genes get tested for prognostic
> value, with no error raised. Two earlier hand-offs exist inside Repo 1
> itself — see that repo's README for the full chain.

---

## What this notebook does (`01_metabric_prognostic_validation.ipynb`)

Single notebook, run top to bottom, in four phases:

**Phase 0 — Load & clean METABRIC clinical data** (Cells 0–3)
- Uploads `data_clinical_patient.txt` + `data_clinical_sample.txt`
  (cBioPortal METABRIC clinical files)
- Filters to the 4 subtypes matching the GSE45827 classifier
  (`CLAUDIN_SUBTYPE`: LumA/LumB/Her2/Basal → Luminal A/Luminal B/Her2/Basal)
- Encodes `OS_STATUS`/`RFS_STATUS` to binary event indicators, drops rows
  missing survival data
- **Result: 1,608-patient cohort** (Luminal A 700, Luminal B 475, Her2 224, Basal 209)

**Phase 1 — Pull biomarker gene expression via cBioPortal API** (Cells 4–7)
- `biomarker_genes` = the 43-symbol panel pasted in from Repo 1's
  `reannotate_probes_R.ipynb`
- Queries `brca_metabric_mrna_median_all_sample_Zscores` (z-scored
  microarray expression) — **42 of 43 genes found** in cBioPortal
- Pivots long-format API response to a patient × gene wide matrix, joins
  to the clinical cohort

**Phase 2 — Kaplan-Meier + Cox survival modeling** (Cells 8–14)
- KM curves by subtype
- **Model 1 (clinical only):** age, lymph nodes positive, NPI, subtype
  dummies → **C-index 0.6567**
- **Model 2 (clinical + genes):** same + biomarker genes (**41 genes**
  enter the design matrix — 1 further drop from the 42 found, likely a
  collinearity/matrix-construction issue worth double-checking if you
  rerun this) → **C-index 0.6752** (+0.0185 over clinical-only)
- Forest plot of gene hazard ratios
- Lead-gene KM stratification: dynamically picks whichever gene has the
  lowest p-value among those significant at p < 0.05 in the just-fit
  model, rather than a hardcoded gene — this run picked **CDCA5**
  (log-rank p = 1.69×10⁻⁹)

**Phase 3 — Sensitivity analysis & multiple-testing correction** (Cells 15–22)
- Schoenfeld residual test flags **ESR1** as violating the proportional-
  hazards assumption; sensitivity model refit with ESR1 removed and
  stratified on age/NPI groups instead → **C-index drops to 0.6097**
  (worth noting: markedly lower than the main model, since it no longer
  uses AGE_AT_DIAGNOSIS/NPI as continuous predictors)
- **Genes with raw p < 0.05 in the main model: CDCA5 (p=0.0155),
  IL23A (p=0.0193), AQP5 (p=0.0228).** CMC2 does **not** reach
  significance in either model (main p=0.353, sensitivity p=0.503)
- **After Benjamini-Hochberg FDR correction, none of the 41 genes reach
  q < 0.05** — CDCA5/IL23A/AQP5 all sit at q≈0.312. This is the
  notebook's own printed output (Cell 21), not a judgment call — worth
  stating plainly in the manuscript rather than reporting only the raw
  p-values, since raw-vs-FDR framing changes how strong the prognostic
  claim can be.

---

## Key numbers (for manuscript cross-checking)

| Metric | Value |
|---|---|
| Cohort size | 1,608 patients |
| Genes requested / found in cBioPortal / used in Cox model | 43 / 42 / 41 |
| C-index, clinical only | 0.6567 |
| C-index, clinical + genes | 0.6752 |
| C-index, sensitivity model (ESR1 removed, PH-stratified) | 0.6097 |
| Raw p < 0.05 | CDCA5, IL23A, AQP5 |
| FDR q < 0.05 | **none** (lowest q = 0.3116) |
| Lead gene, KM stratification | CDCA5 (log-rank p = 1.69×10⁻⁹) |
| PH-violating covariate | ESR1 |

---

## Known rough edges

1. **Manual gene-list paste-in from Repo 1** (see warning banner above)
   — the single biggest reproducibility risk in this notebook.

2. **Checkpoint-reload pattern.** Cell 15 reloads `merged_clean` from a
   saved `metabric_survival_analysis.csv` and rebuilds `biomarker_genes`
   from scratch mid-notebook, rather than relying on in-memory state from
   Cells 0–9. This is safe as written (the retyped list in Cell 15
   matches Cell 5), but it's a second copy of the same 43-gene list in
   one notebook — if one gets edited and not the other, they'll silently
   diverge. Worth consolidating into a single source-of-truth variable
   or a small helper function.

3. **43 → 42 → 41 attrition isn't fully explained in-notebook.** The
   42→41 drop between "found in cBioPortal" and "genes in Cox model" — 
   likely one gene is being excluded when `gene_cols` is constructed —
   isn't logged anywhere. Worth adding a `print(set(found_genes) -
   set(gene_cols))`-style check if this notebook is rerun, so the
   missing gene is identified by name rather than just by count.

4. **Raw-p vs. FDR-q framing.** Cells 11–14 (forest plot, C-index bar
   chart) report the model built on raw p < 0.05 significance; Cell 21's
   FDR correction (run afterward) shows none of those genes survive
   correction. Both are legitimate to report, but the manuscript should
   be explicit about which standard is being used where, rather than
   citing "significant" without qualification.

---

## Environment

- Google Colab (Python)
- Packages: `requests` (cBioPortal API), `pandas`, `numpy`, `lifelines`
  (KM, Cox, Schoenfeld PH test, log-rank test), `matplotlib`,
  `statsmodels` (Benjamini-Hochberg FDR correction)
- Data: `data_clinical_patient.txt` / `data_clinical_sample.txt` from
  cBioPortal's `brca_metabric` study (upload manually or fetch via API)

## Upstream dependency

Gene panel origin: `breast-cancer-subtype-classification` repo →
`03_feature_stability_analysis.ipynb` → `reannotate_probes_R.ipynb` →
`model_B_top50_reannotated.csv` → manually pasted into Cell 5 here.
See that repo's README for the full derivation chain and internal-CV
accuracy (93.08% ± 5.10%) the panel is built on top of.

## Data Availability

| Dataset | Source | Samples | Use |
|---|---|---|---|
| GSE45827 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45827) | 130 tumor samples | Upstream discovery cohort — classifier training + SHAP gene selection (see companion repository) |
| METABRIC | [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) | 1,608 patients | Independent survival validation cohort (this repository) |

## Repository Structure

```
breast-cancer-survival-biomarkers/
├── README.md
├── requirements.txt
├── .gitignore
├── notebooks/
│   └── 01_metabric_survival_validation.ipynb   # Cox regression, KM curves, PH testing, sensitivity analysis, C-index
├── data/
│   ├── metabric_survival_analysis.csv          # Processed METABRIC cohort (1,608 patients x 73 columns)
|   ├──  data_clinical_patient.txt
|   └──  data_clinical_sample.txt
├── results/                                     # Cohort characteristics, Cox results, Schoenfeld PH test results
└── figures/                                     # Figs 1-4 as referenced in the manuscript
```

**Note on raw clinical files:** `data_clinical_patient.txt` and `data_clinical_sample.txt` are large raw cBioPortal exports. Downloaded directly from [cBioPortal METABRIC](https://www.cbioportal.org/study/summary?id=brca_metabric) ; `metabric_survival_analysis.csv` above already contains the fully processed result and is sufficient for reproducing every table and figure in the manuscript without re-downloading anything.

## Results

**Subtype survival separation:** Kaplan-Meier curves confirm the expected PAM50 survival hierarchy in METABRIC (log-rank p=1.70×10⁻⁹). The Basal subtype shows a characteristic crossing-hazard pattern — poor early survival with relative stabilization beyond 250 months.

**Independently prognostic genes:**

| Gene | Hazard Ratio | 95% CI | p-value | Direction |
|---|---|---|---|---|
| CDCA5 | 1.144 | 1.019-1.284 | 0.023 | Higher expression → worse survival |
| CMC2 | 0.919 | 0.845-1.000 | 0.049 | Higher expression → better survival |

Both significant after adjusting for age, lymph nodes, NPI, and PAM50 subtype.

**Model comparison:**

| Model | C-index |
|---|---|
| Clinical only | 0.657 |
| Clinical + ML genes | 0.679 |
| Improvement | +0.022 |

## Reproducing This Work

```
git clone https://github.com/Zohaib-Bioinfo/breast-cancer-survival-biomarkers.git
cd breast-cancer-survival-biomarkers
pip install -r requirements.txt
```

Open `notebooks/01_metabric_survival_validation.ipynb` in Google Colab, Kaggle, or any Python 3.12 environment with `lifelines`, `pandas`, `numpy`, `requests`, `matplotlib`, and `seaborn` installed. The processed METABRIC CSV in `data/` allows full reproduction without re-querying the cBioPortal API, though the notebook includes the original API retrieval code for transparency. If you want to reproduce the raw-data acquisition step itself rather than start from the processed CSV, download `data_clinical_patient.txt` and `data_clinical_sample.txt` manually from cBioPortal (see Data Availability above) — these are not included in the repository due to file size.

## Limitations

- Cox model uses L2 penalization (penalizer=0.1) to handle 46 simultaneous gene predictors — results should be considered exploratory
- Expression platforms differ between GSE45827 (Affymetrix microarray) and METABRIC (Illumina HT-12 microarray) — cross-platform normalization not applied
- No experimental (wet-lab) validation of CDCA5 or CMC2 prognostic roles
- Median split for CDCA5 KM analysis is a post-hoc visualization, not a pre-specified cutpoint
- Neither CDCA5 nor CMC2 survived Benjamini-Hochberg FDR correction across all 46 genes tested; see companion repository for convergent robustness evidence differentiating confidence between the two

## Citation

If you use this analysis or its outputs, please cite the associated manuscript (full text and supporting information available in the [companion repository](https://github.com/Zohaib-Bioinfo/breast-cancer-subtype-classification)).

## Author

Muhammad Zohaib — BS Bioinformatics, Department of Computer Science, University of Agriculture Faisalabad (UAF)

## License

MIT License — see `LICENSE`
