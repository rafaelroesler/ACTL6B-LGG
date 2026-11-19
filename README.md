# ACTL6B expression and neuronal differentiation in lower-grade glioma

This repository organizes the original R analysis scripts associated with the ACTL6B TCGA-LGG manuscript. The files have been **renamed, but their code has not been altered**; see `SOURCE_MANIFEST.csv` for provenance and SHA-256 checksums.

## Scripts and workflow

| Script | Analysis |
|---|---|
| `01_figures_1_2_workspace.R` | Recreates/exports original Figure 1–2 objects from a validated workspace snapshot |
| `02_cox_initial_models.R` | Initial Cox survival models |
| `03_cox_molecular_adjusted.R` | Age-, grade-, molecular-class-adjusted Cox model |
| `04_cox_stratified_sensitivity.R` | Cox sensitivity model stratified by grade and molecular class |
| `05_cox_forest_plot.R` | Figure 3 forest plot, using previously transcribed model estimates |
| `06_coexpression_four_heatmaps.R` | ACTL6B co-expression and targeted heatmaps |
| `07_idh_stratified_coexpression.R` | IDH-mutant versus IDH-wildtype co-expression; standalone input-loading workflow |
| `08_baf_coexpression.R` | BAF chromatin-remodeling gene expression and ACTL6A/ACTL6B relationships |
| `09_go_bp_gsea.R` | GO Biological Process GSEA |
| `10_random_forest.R` | Supervised Random Forest regression with repeated cross-validation |
| `11_xcell.R` | xCell enrichment analysis |
| `12_estimate_scores.R` | ESTIMATE stromal/immune signature scoring |
| `13_estimate_figure.R` | ACTL6B ESTIMATE association figure from score outputs |

## Requirements and reproducibility limitations

- Original analyses were conducted in R/RStudio. Several scripts expect the existing RStudio session to contain `expr_lgg`, `analysis_data`, and sometimes `lgg_subtypes`; these are **not included** here.
- Scripts use hard-coded Windows paths under `C:/ACTL6B_TCGA_LGG`; configure the same folder structure before sourcing, or edit paths for your environment.
- `01_figures_1_2_workspace.R` and `02_cox_initial_models.R` / `03_cox_molecular_adjusted.R` require `C:/ACTL6B_TCGA_LGG/ACTL6B_Figures1_2_verified_workspace.rds` (not supplied).
- `04_cox_stratified_sensitivity.R` requires the patient-level CSV produced by `03_cox_molecular_adjusted.R`.
- `13_estimate_figure.R` requires score files produced by `12_estimate_scores.R`.
- `07_idh_stratified_coexpression.R` can load Xena expression data, but requires molecular subtype annotation from the local CSV or TCGAbiolinks.
- `05_cox_forest_plot.R` uses estimates transcribed from a prior model output rather than reading fresh Cox coefficients; cross-check against the model before publication.
- Package installation and input data acquisition are not automated consistently across scripts. Consult the package checks at the beginning of each script.
- The scripts have been inspected and organized, **not rerun or independently validated** in this environment. The figures and numerical results must be checked against the manuscript before public release.
- No clinical patient-level files or large expression datasets are distributed in this archive. Add source accessions, data availability, package versions, and licensing details before publishing.

## Suggested run order

1. Restore the validated TCGA-LGG data and RStudio objects.
2. Run `03_cox_molecular_adjusted.R` before `04_cox_stratified_sensitivity.R`.
3. Run `12_estimate_scores.R` before `13_estimate_figure.R`.
4. Run other scripts individually after checking their input requirements.

`results/` and `figures/` are placeholders for reviewed outputs. The repository does not claim to include every analysis performed in the manuscript.
