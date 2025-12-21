# postpartum-depression-Risk-Prediction-ML
Integrated postpartum depression study (Bangladesh): burden &amp; social gradients, TMLE causal effects of modifiable determinants, ML risk prediction (elastic net + Super Learner), fairness/subgroup audits, and policy simulations for risk-targeted prevention. Reproducible R pipeline + manuscript tables/figures.

# Postpartum Depression in Bangladesh — Causal Effects, Risk Prediction, Fairness, and Policy Simulations

Reproducible R workflow for a full end-to-end study of postpartum depression (PPD) risk in Bangladesh, spanning: (1) cohort description and data quality, (2) burden and social patterning, (3) causal effects of modifiable determinants using Targeted Maximum Likelihood Estimation (TMLE), (4) machine-learning risk prediction (elastic net + Super Learner) with calibration, (5) subgroup/fairness performance and explainability (SHAP), and (6) sensitivity analyses and policy-relevant risk-targeting simulations.

---

## Project goals (what this repo delivers)

**Primary objectives**
1. **Quantify PPD burden** using EPDS-defined high-risk status.
2. **Characterize social and clinical gradients** in PPD prevalence.
3. **Estimate causal effects** of modifiable exposures on PPD risk (TMLE; contrasted with GLM).
4. **Develop and evaluate prediction models** with cross-validation and calibration.
5. **Assess subgroup performance and fairness gaps** across key strata.
6. **Simulate policy targeting** (prevented cases under interventions applied to highest predicted-risk groups).

**Primary outcome**
- **PPD high-risk status** defined as **EPDS ≥ 13** (binary), with EPDS total score also analyzed descriptively.

---

## Data source

This study uses a publicly available, de-identified dataset:

Raisa, J. F., & Kaiser, M. S. (2025). *Data for Postpartum Depression Prediction in Bangladesh* (Version 2) [Data set]. Mendeley Data. https://doi.org/10.17632/4nznnrk8cg.2

Local copies used in this repo:
- `data/PPD_dataset_v2.csv`
- `data/PPD_Data Dictionary_v2.csv`

---

## Repository structure

.
├── data/
│ ├── PPD_dataset_v2.csv
│ └── PPD_Data Dictionary_v2.csv
├── scripts/
│ ├── Section 1 – Participant Characteristics & Data Quality.R
│ ├── Section 2 – Burden & Social Patterning of Postpartum Depression.R
│ ├── Section 3 – Causal Effects of Modifiable Risk Factors.R
│ ├── Section 4 – Development and Performance of Prediction Models.R
│ ├── Section 5 – Fairness, Subgroup Performance, and Explainability.R
│ └── Section 6 – Sensitivity Analyses and Policy-Relevant Simulations.R
└── output/
├── section1_descriptives/
│ ├── figures/
│ │ ├── figure1_epds_distribution_by_ppd.png
│ │ └── figure2_age_density_by_ppd.png
│ └── tables/
│ ├── table1_baseline_characteristics.html
│ └── missingness_summary_section1.csv
├── section2_causal/
│ ├── figures/
│ │ └── figure_ppd_prevalence_gradients.png
│ └── tables/
│ ├── table2_ppd_prevalence_by_social_patterning.html
│ ├── overall_ppd_prevalence_section2.csv
│ └── ppd_prevalence_by_group_section2.csv
├── section3_prediction/
│ ├── figures/
│ │ └── figure2_tmle_forest_relative_risks.png
│ └── tables/
│ ├── table3_tmle_vs_glm_causal_effects.html
│ └── tmle_vs_glm_causal_results_section3.csv
├── section4_prediction/
│ ├── figures/
│ │ └── figure3_roc_and_calibration_curves.png
│ └── tables/
│ ├── table4_prediction_performance.html
│ ├── elastic_net_full_coefficients.csv
│ └── superlearner_full_weights.csv
├── section5_fairness/
│ ├── figures/
│ │ ├── figure4_shap_global_importance.png
│ │ └── figure5_subgroup_fairness_sensitivity.png
│ └── tables/
│ ├── table5_subgroup_performance.html
│ ├── subgroup_performance_section5.csv
│ ├── fairness_gaps_section5.csv
│ ├── global_importance_shap_section5.csv
│ └── global_importance_shap_or_coef_section5.csv
└── section6_sensitivity_policy/
├── figures/
│ └── figure6_risk_targeting_policy_simulation.png
└── tables/
├── table6_sensitivity_and_policy_summary.html
├── sensitivity_and_policy_summary_section6.csv
└── policy_simulation_curve_points_section6.csv

markdown
Copy code

---

## Outputs (what to look at first)

### Figures (key visuals)
- **Figure 1:** EPDS score distribution by PPD status  
  `output/section1_descriptives/figures/figure1_epds_distribution_by_ppd.png`
- **Figure 2:** Maternal age distribution by PPD status  
  `output/section1_descriptives/figures/figure2_age_density_by_ppd.png`
- **Figure 3:** Social/clinical gradients in PPD prevalence  
  `output/section2_causal/figures/figure_ppd_prevalence_gradients.png`
- **Figure 4:** TMLE causal effects (relative risks; log scale)  
  `output/section3_prediction/figures/figure2_tmle_forest_relative_risks.png`
- **Figure 5:** Prediction model ROC + calibration curves  
  `output/section4_prediction/figures/figure3_roc_and_calibration_curves.png`
- **Figure 6:** Global feature importance (SHAP)  
  `output/section5_fairness/figures/figure4_shap_global_importance.png`
- **Figure 7:** Subgroup sensitivity (TPR at a fixed threshold)  
  `output/section5_fairness/figures/figure5_subgroup_fairness_sensitivity.png`
- **Figure 8:** Policy simulations: prevented cases vs targeted proportion  
  `output/section6_sensitivity_policy/figures/figure6_risk_targeting_policy_simulation.png`

### Tables (manuscript-ready HTML)
Open HTML tables in your browser:
- **Table 1:** Baseline characteristics  
  `output/section1_descriptives/tables/table1_baseline_characteristics.html`
- **Table 2:** PPD prevalence by social patterning  
  `output/section2_causal/tables/table2_ppd_prevalence_by_social_patterning.html`
- **Table 3:** TMLE vs GLM causal effects  
  `output/section3_prediction/tables/table3_tmle_vs_glm_causal_effects.html`
- **Table 4:** Prediction performance  
  `output/section4_prediction/tables/table4_prediction_performance.html`
- **Table 5:** Subgroup performance and fairness gaps  
  `output/section5_fairness/tables/table5_subgroup_performance.html`
- **Table 6:** Sensitivity analyses + policy summary  
  `output/section6_sensitivity_policy/tables/table6_sensitivity_and_policy_summary.html`

---

## Methods overview (high-level)

### Descriptives and burden
- EPDS total score summarized and visualized by outcome group.
- PPD prevalence reported overall and stratified by key sociodemographic/clinical factors.

### Causal inference (TMLE)
- TMLE estimates the causal effect of modifiable exposures on PPD risk with flexible nuisance models.
- Effects are reported as **relative risks** with uncertainty estimates; GLM estimates are included for comparison.

### Prediction modeling
- **Elastic net** (regularized logistic regression) for sparse, stable prediction.
- **Super Learner** ensemble for flexible, data-adaptive risk prediction.
- Performance assessed with discrimination (ROC/AUROC) and calibration plots.

### Fairness & subgroup performance
- Stratified performance metrics across predefined subgroups.
- Fairness gaps summarized and visualized.
- **SHAP** used for global model interpretability and feature importance.

### Policy simulations
- Translate causal and predictive results into decision metrics:
  “If we target the highest predicted-risk X% of women, what fraction of cases might be prevented under plausible interventions?”

---

## Quickstart (reproduce the full pipeline)

### 1) Requirements
- R (recommended: ≥ 4.2)
- RStudio (optional but recommended)

### 2) Install packages
Open R and install core packages used across sections:
```r
install.packages(c(
  "tidyverse","data.table","readr","stringr","forcats",
  "tableone","janitor","skimr",
  "glmnet","pROC","yardstick","caret",
  "SuperLearner",
  "ggplot2","patchwork","scales",
  "rmarkdown","knitr"
))
Note: Some scripts may require additional packages (e.g., TMLE-related or SHAP tooling). If a script errors with “package not found,” install the missing package and re-run.

3) Run scripts in order (recommended)
From the project root:

scripts/Section 1 – Participant Characteristics & Data Quality.R

scripts/Section 2 – Burden & Social Patterning of Postpartum Depression.R

scripts/Section 3 – Causal Effects of Modifiable Risk Factors.R

scripts/Section 4 – Development and Performance of Prediction Models.R

scripts/Section 5 – Fairness, Subgroup Performance, and Explainability.R

scripts/Section 6 – Sensitivity Analyses and Policy-Relevant Simulations.R

Each script writes results into its corresponding output/section*_.../figures and output/section*_.../tables folders.

Reproducibility notes
Ensure your working directory is the repo root (so relative paths resolve).

If you clone this repo, keep folder names unchanged.

HTML tables are designed to be manuscript-ready. If you need Word/PDF tables, render HTML through your preferred workflow or convert using RMarkdown.


Authors and contributions (CRediT)
Sunday A. Adetunji (Lead Author): Conceptualization; Methodology (statistical and ML); Software; Formal analysis; Data curation; Visualization; Writing—original draft (all sections); Writing—review & editing; Supervision; Project administration.

Tosin F. Alabi: Conceptualization (study design); Writing—review & editing (discussion and interpretation).

Tomiisin B. Balogun: Conceptualization (study design; policy translation); Writing—original draft (policy framing); Writing—review & editing (policy-focused revisions).

Citation
Cite the dataset
Raisa, J. F., & Kaiser, M. S. (2025). Data for Postpartum Depression Prediction in Bangladesh (Version 2) [Data set]. Mendeley Data. https://doi.org/10.17632/4nznnrk8cg.2



Contact
Corresponding author: Sunday A. Adetunji, MD
Email: adetunjs@oregonstate.edu
ORCID: https://orcid.org/0000-0001-9321-9957

Abbreviations
AUROC (area under the receiver operating characteristic curve); CI (confidence interval); CRediT (Contributor Roles Taxonomy); EPDS (Edinburgh Postnatal Depression Scale); GLM (generalized linear model); IPV (intimate partner violence); ML (machine learning); PHQ-9 (Patient Health Questionnaire–9); PPD (postpartum depression); RR (relative risk); ROC (receiver operating characteristic); SHAP (Shapley additive explanations); TMLE (targeted maximum likelihood estimation).

