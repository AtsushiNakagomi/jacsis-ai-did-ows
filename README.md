# JACSIS generative-AI and subjective wellbeing — four-wave panel analysis

R analysis script for:

> Nakagomi, A., Akutsu, Y., Yasuoka, M. & Tabuchi, T. (2026). Generative AI use and health and wellbeing in Japanese adults: a longitudinal outcome-wide difference-in-differences analysis. *[Journal]*. DOI: [will be added on acceptance].

This repository contains the analysis code and variable codebook used to generate the results, tables, and forest plots reported in the manuscript. The study uses four consecutive waves of the Japan COVID-19 and Society Internet Survey (JACSIS, 2022 / 2023 / 2024 / 2025) to estimate the effect of generative-AI use initiation on a 13-item panel of subjective-wellbeing outcomes drawn from the VanderWeele Flourishing Index, K6 distress, LSNS-6 social network, UCLA-3 loneliness, and time-use.

## Repository contents

| File | Description |
|---|---|
| `analysis.R` | Main analysis script. Constructs all exposures, confounders, and outcomes; runs the four-tier analytic pipeline (Main DiD, Sens 1 DiD, Sens 2 ANCOVA, Sens 3 ANCOVA); fits exploratory factor analysis and VIF diagnostics; produces all manuscript tables (Table 1-3, Sup Tables 1A/1B/2/3/4) and per-tier forest plots. |
| `CODEBOOK.md` | Mapping between JACSIS survey items and analytic variables used in the script. |
| `LICENSE` | MIT License. |
| `README.md` | This file. |

## Data availability

The JACSIS data are not publicly available due to ethical restrictions on participant privacy. De-identified data may be shared upon reasonable request, subject to approval by the JACSIS steering committee and relevant ethics review boards. See the manuscript Data Availability section for details.

`analysis.R` expects three panel CSV files at `data/`. These must be supplied by approved researchers; column names must match those documented in `CODEBOOK.md`. Each file is the wide-format merge of the listed JACSIS waves on `Monitor_ID`, with each variable suffixed by the wave year (e.g., `Q37S1_2025`, `AGE_2024`).

| File | Waves merged | Tier(s) that read it |
|---|---|---|
| `data/jacsis_2022_2025.csv`         | 2022 + 2025               | Sens 3|
| `data/jacsis_2023_2024_2025.csv`    | 2023 + 2024 + 2025        | Main, Sens 2|
| `data/jacsis_2022_2023_2025.csv`    | 2022 + 2023 + 2025        | Sens 1|

## Analytic design

The pipeline fits 13 wellbeing outcomes × 5 AI-use exposures (`use` = binary initiation indicator; `intensity` = mean of nine purpose-frequency items; `creative` / `daily` / `social` = three purpose-composite scores fitted jointly) across four tiers:

| Tier | Model | Cohort | Identification | Outcome window |
|---|---|---|---|---|
| **Main**   | DiD (first-differences) | `Q37S1_2025 ∈ {5, 6}` vs `== 1` | parallel-trends test on 2023 → 2024 | ΔY 2024 → 2025 |
| **Sens 1** | DiD                     | `Q37S1_2025 ∈ {4, 5, 6}` vs `== 1` | parallel-trends test on 2022 → 2023 | ΔY 2023 → 2025 |
| **Sens 2** | ANCOVA                  | same cohort as Main | conditional independence given `Y_2024` + covariates | `Y_2025` given `Y_2024` |
| **Sens 3** | ANCOVA                  | `Q37S1_2025 == 1` vs `∈ {3, 4, 5, 6}` (code 2 excluded) | conditional independence given `Y_2022` + covariates | `Y_2025` given `Y_2022` |

All four tiers adjust for nine baseline confounders (sex, age in 5-year bands, education, employment, income, marital, living alone, smartphone time, PC/tablet time) plus the levels of the 12 other outcomes at the tier's baseline wave. Standard errors are OLS with HC1 robust covariance. Multiple testing is controlled by Benjamini–Hochberg FDR within each exposure family (13 tests per family). Coefficients are reported on both the raw outcome scale and standardised by the baseline-wave outcome SD (`estimate_oz = estimate / sd_pre`).

See `CODEBOOK.md` for full variable definitions, outcome constructs, cross-wave Q-code crosswalks, and pre-/post-/baseline conventions.

## Requirements

- R version 4.5.x
- R packages: `readr`, `dplyr`, `tidyr`, `purrr`, `sandwich`, `lmtest`, `psych`, `car`, `ggplot2`

Install missing packages from CRAN before running:

```r
install.packages(c("readr","dplyr","tidyr","purrr","sandwich",
                   "lmtest","psych","car","ggplot2"))
```

## Usage

1. Place the three panel CSV files at the paths listed under **Data availability**.
2. From the repository root, run:
   ```
   Rscript analysis.R
   ```
3. Outputs are written to `output/`.

## Outputs

After successful execution, `output/` will contain:

```
output/
└── manuscript/
    ├── table1_main.{csv,md}                              ← Table 1
    ├── table2_main_results.{csv,md}                      ← Table 2 (with E-values)
    ├── table3_cross_tier_ladder.{csv,md}                 ← Table 3
    ├── sup_table_1a_main_purposes_anyuse.{csv,md}        ← Supplementary Table 1A
    ├── sup_table_1b_main_purposes_daily.{csv,md}         ← Supplementary Table 1B
    ├── sup_table_2_sens1_characteristics.{csv,md}        ← Supplementary Table 2
    ├── sup_table_3_sens3_characteristics.{csv,md}        ← Supplementary Table 3
    ├── sup_table_4_main_raw_levels_changes.{csv,md}      ← Supplementary Table 4
    ├── diagnostics/
    │   ├── efa_loadings.csv                              ← EFA factor loadings
    │   ├── efa_fit.csv                                   ← TLI, RMSEA + 90% CI, BIC, cumulative variance, Cronbach's alphas
    │   └── vif.csv                                       ← Variance inflation factors (Main DiD design matrix)
    └── tables_figures/
        ├── main/
        │   ├── results_main_did.csv                      ← 65 rows = 13 outcomes × 5 exposures (Main DiD; with E-values)
        │   ├── results_main_did_usersonly.csv            ← Users-only A + B sensitivity (52 rows)
        │   ├── results_main_pretrends.csv                ← Pre-trends test (ΔY 2023 → 2024)
        │   ├── forest_main_did.png                       ← Forest plot (raw beta)
        │   ├── forest_main_did_oz.png                    ← Forest plot (outcome-standardised beta)
        │   ├── forest_main_did_usersonly{,_oz}.png       ← Users-only A + B forest plots
        │   └── forest_main_pretrends{,_oz}.png           ← Pre-trends forest plots
        ├── sens1/                                        ← Sens 1 DiD parallel set
        ├── sens2/                                        ← Sens 2 ANCOVA parallel set
        └── sens3/                                        ← Sens 3 ANCOVA parallel set
```

## Reproducibility

- Random seed: `set.seed(20260512)` is set at the top of the script.
- All statistical estimates use closed-form formulas (no Monte-Carlo step), so results are deterministic given fixed input data.
- The exploratory factor analysis uses `psych::fa()` with the maximum-likelihood extraction and oblimin rotation; this is also deterministic for a fixed input matrix.

## Citation

If you use or adapt this code, please cite:

> Nakagomi, A., Akutsu, Y., Yasuoka, M. & Tabuchi, T. (2026). Generative AI use and health and wellbeing in Japanese adults: a longitudinal outcome-wide difference-in-differences analysis. *[Journal]*. DOI: [will be added on acceptance].

## License

The code in this repository is released under the MIT License. See `LICENSE` for details.

## Contact

Atsushi Nakagomi — a.nakagomi@chiba-u.jp
Center for Preventive Medicine, Chiba University
ORCID: [0000-0002-3908-696X](https://orcid.org/0000-0002-3908-696X)
