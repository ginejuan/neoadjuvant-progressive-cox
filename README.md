# Confounding by indication in the evaluation of neoadjuvant chemotherapy

[![DOI](https://zenodo.org/badge/1293370566.svg)](https://doi.org/10.5281/zenodo.21259046)

Reproducible analysis code and de-identified data for the manuscript
*"Confounding by indication in the evaluation of neoadjuvant chemotherapy:
a sequential Cox analysis of an eight-year single-centre breast cancer cohort"*
(Clinical and Translational Oncology, under review).

Single-centre retrospective cohort of 461 women with breast cancer
(Puerto Real University Hospital, 2016–2018; median follow-up 7.3 years,
reverse Kaplan–Meier).

## Contents

```
code/reproducible_analysis.py      Reproduces every result in the manuscript
code/data_recovery_extractors.py   Source-document extraction functions (documentation)
data/dataset_deidentified.csv      De-identified dataset (N = 461; no personal data)
requirements.txt                   Python dependencies
```

## Reproduce the results

```bash
pip install -r requirements.txt
cd code
python reproducible_analysis.py      # reads ../data/dataset_deidentified.csv
```

The script regenerates the descriptive statistics, Kaplan–Meier estimates,
multivariable and sequential (progressive) Cox models, proportional-hazards
tests and time-split analyses, baseline standardised mean differences, the
power and sensitivity analyses (post-hoc power, stage-IV exclusion,
cause-specific analysis, Ki67 model, HER2-negative restriction), the
revision analyses (6/12-month landmark analysis for radiotherapy,
common-sample progressive models, multiply-imputed full adjustment,
time-dependent stability of the NACT estimate, time origin at diagnosis)
and the E-values.

## Data

`dataset_deidentified.csv` contains one row per patient with survival times,
event indicators and model covariates only — no identifiers, names or dates.
It includes `her2pos` (HER2-positive surrogate subtype) and `delay_dx_days`
(biopsy-to-surgery interval in days, for the time-origin sensitivity
analysis). Covariate values recovered from source documents (grade, Ki67,
age, radiotherapy) are already integrated (see manuscript Supplementary
Figure S1).

## Citation

Fernández Alba JJ, et al. Confounding by indication in the evaluation of
neoadjuvant chemotherapy: a sequential Cox analysis of an eight-year
single-centre breast cancer cohort. *Clinical and Translational Oncology*
(under review). Code and data archived at Zenodo:
https://doi.org/10.5281/zenodo.21259046

## License

Code released under the MIT License (see `LICENSE`); the de-identified
dataset is released under CC BY 4.0 (see `data/LICENSE`).
