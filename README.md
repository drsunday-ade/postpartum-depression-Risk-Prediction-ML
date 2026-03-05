# A parsimonious machine-learning model for postpartum depression risk from prenatal PHQ-2 symptom trajectories: internal validation using EPDS ≥13

Reproducible **R** analysis pipeline and manuscript-ready outputs (tables/figures) for a **prediction-model development study with internal validation** using a publicly available, de-identified postpartum cohort from **Bangladesh** (n=800; births within 24 months). The primary endpoint is **high depressive symptom burden (EPDS ≥13)**, and the core prognostic feature is a prespecified **prenatal PHQ-2 trajectory** (before pregnancy × during pregnancy: NN, NP, PN, PP).

---

## Manuscript
**Adetunji SA**, Balogun T, Oyewusi R.  
**A parsimonious machine-learning model for postpartum depression risk from prenatal PHQ-2 symptom trajectories: internal validation using EPDS ≥13.** *Manuscript.*

### Authors & affiliations
- **Sunday A. Adetunji, MD, MPH** (1,2)  
  1. College of Health, Oregon State University, Corvallis, OR, USA  
  2. College of Health, Obafemi Awolowo University, Ile-Ife, Nigeria  
- **Tomiisin Balogun, BSc** (3)  
  3. Joseph Ayo Babalola University, Arakeji, Osun State, Nigeria  
- **Rhoda Oyewusi, RN, RM, PON** (4)  
  4. University of Lagos, Department of Nursing/Midwifery, Lagos, Nigeria  

**Corresponding author:** Sunday A. Adetunji — adetunjs@oregonstate.edu — ORCID: 0000-0001-9321-9957

---

## Data
This repository does **not** redistribute the dataset.

**Public dataset:** Raisa JF, Kaiser MS (2025). *Data for Postpartum Depression Prediction in Bangladesh*. Mendeley Data, V2. doi:10.17632/4nznnrk8cg.2

### Expected local path
After downloading from Mendeley Data, place the CSV here:

data/PPD_dataset_v2.csv


---

## Methods (high-level)
- **Outcome:** EPDS-high = **EPDS ≥13** (primary); sensitivity analyses for EPDS ≥12 and EPDS ≥11
- **Core predictor:** prenatal **PHQ-2 trajectory** (NN, NP, PN, PP)
- **Primary models:** nested **elastic-net penalized logistic regression**
  - Model A: trajectory only  
  - Model B: + demographics/household  
  - Model C: + pregnancy context  
  - Model D: + psychosocial/postpartum context (as available)
- **Benchmark (secondary):** XGBoost (Model D) with probability recalibration (where implemented)
- **Internal validation:** stratified bootstrap optimism correction
- **Performance:** AUC, Brier score, calibration intercept/slope, and decision-curve analysis (0.05–0.60 thresholds)

---

## Reproducibility: quick start
1) Clone:
```bash
git clone https://github.com/drsunday-ade/postpartum-depression-Risk-Prediction-ML.git
cd postpartum-depression-Risk-Prediction-ML

Add the dataset:

data/PPD_dataset_v2.csv

Restore/install dependencies (choose what matches the repo setup):

If renv.lock is present:

renv::restore()

Otherwise, install required packages listed in the project scripts.

Run the analysis:

Open the R project and run the main pipeline script located in the repository (see the analysis/ or src/ directory).

The pipeline generates manuscript-ready tables/figures in the outputs directory (see below).

Tip: If you want a single-command execution, add an explicit entrypoint (e.g., analysis/run_pipeline.R) and call it via Rscript.

Outputs (typical)

outputs/tables/ — descriptive tables, performance tables (AUC/Brier/calibration), subgroup summaries

outputs/figures/ — EPDS distribution, PHQ-2 trajectory risk gradient, ROC, calibration, decision curves

manuscript/ — manuscript text and supplementary materials (if included)

Ethics

This study is a secondary analysis of a publicly available, de-identified dataset. The dataset creators report that original data collection obtained ethical approval and informed consent. No participant contact occurred in this study.

Funding & competing interests

Funding: None.

Competing interests: None declared.

How to cite
Our paper

Adetunji SA, Balogun T, Oyewusi R. A parsimonious machine-learning model for postpartum depression risk from prenatal PHQ-2 symptom trajectories: internal validation using EPDS ≥13. Manuscript.

Data source

Raisa JF, Kaiser MS. Data for Postpartum Depression Prediction in Bangladesh. Mendeley Data, V2. doi:10.17632/4nznnrk8cg.2

Contact

Sunday A. Adetunji, MD, MPH
adetunjs@oregonstate.edu
 | sundayadetunjisa@gmail.com

ORCID: https://orcid.org/0000-0001-9321-9957
