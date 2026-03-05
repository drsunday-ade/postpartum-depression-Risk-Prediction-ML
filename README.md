# Postpartum Depression Risk Prediction (Bangladesh) — Reproducible ML Pipeline

This repository contains the reproducible analysis code for a prediction-model development study with internal validation using a publicly available, de-identified postpartum dataset from Bangladesh (women with a birth within the prior 24 months). The primary objective is to estimate postpartum depression (PPD) risk using a parsimonious machine-learning approach derived from prenatal PHQ-2 symptom trajectories, with outcome defined as EPDS ≥13.

## How to cite

### Our paper (this repository)
Adetunji SA, Balogun T, Oyewusi RO. **A parsimonious machine-learning model for postpartum depression risk from prenatal PHQ-2 symptom trajectories: internal validation using EPDS ≥13.** *Manuscript*.

**BibTeX**
```bibtex
@article{adetunji_ppd_parsimonious_ml,
  title   = {A parsimonious machine-learning model for postpartum depression risk from prenatal PHQ-2 symptom trajectories: internal validation using EPDS \\ge 13},
  author  = {Adetunji, Sunday A. and Balogun, Tomiisin and Oyewusi, Rhoda O.},
  journal = {Manuscript},
  note    = {Code: https://github.com/drsunday-ade/postpartum-depression-Risk-Prediction-ML}
}
Data source

Raisa JF, Kaiser MS. Data for Postpartum Depression Prediction in Bangladesh. Mendeley Data, V2. doi:10.17632/4nznnrk8cg.2

Note: This repository does not redistribute the dataset. Please obtain it directly from Mendeley Data.

Data

Download the dataset from Mendeley Data (doi:10.17632/4nznnrk8cg.2).

Place the file here:

data/PPD_dataset_v2.csv
Quick start
git clone https://github.com/drsunday-ade/postpartum-depression-Risk-Prediction-ML.git
cd postpartum-depression-Risk-Prediction-ML

Then run the main analysis entry point in this repository (e.g., analysis/ or src/).
If you provide an environment file (recommended), use one of the following and run the pipeline:

R (renv)

renv::restore()

Python (requirements.txt)

python -m venv .venv
source .venv/bin/activate  # (Windows: .venv\Scripts\activate)
pip install -r requirements.txt
Methods (high-level)

Prediction-model development with internal validation

Feature engineering including prenatal PHQ-2 symptom trajectory variables (as specified in the derived-variable documentation/scripts)

Model training and evaluation focused on:

Discrimination (e.g., AUC)

Calibration (e.g., calibration curve / calibration-in-the-large where implemented)

Transparent reporting of preprocessing and derived variables

Outputs

Typical outputs include:

Performance metrics (discrimination and calibration)

Model objects (where saved)

Reproducible tables/figures for manuscript reporting

Derived-variable specification (if included in the repo)

Repository structure (recommended)

data/ — local data (not tracked; place PPD_dataset_v2.csv here)

src/ or R/ — reusable functions

analysis/ — main pipeline scripts/notebooks

outputs/ — tables, figures, and logs

docs/ — derived-variable specification and reporting notes

Ethics and transparency

This work is a secondary analysis of a publicly available, de-identified dataset. The dataset creators report that original data collection obtained ethical approval and informed consent. No participant contact occurred in this study.

Funding and competing interests

Funding: None.

Competing interests: None declared.

License

Add a LICENSE file (e.g., MIT) to specify reuse terms.

Contact

Sunday A. Adetunji, MD, MPH
Email: sundayadetunjisa@gmail.com

ORCID: https://orcid.org/0000-0001-9321-9957

Acknowledgments

We thank Raisa Jasiya Fairiz and Kaiser M Shamim for making the Bangladesh postpartum depression prediction dataset publicly available through Mendeley Data.
