# Progressive Cox models expose confounding by indication in neoadjuvant chemotherapy

Reproducible analysis code and de-identified data for the manuscript
*"Progressive Cox models expose confounding by indication in neoadjuvant
chemotherapy: an eight-year single-centre breast cancer cohort"*
(Clinical and Translational Oncology).

Single-centre retrospective cohort of 461 women with breast cancer
(Puerto Real University Hospital, 2016–2018; median follow-up 8.0 years).

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
multivariable and progressive Cox models (Tables 2–5), proportional-hazards
tests (Table 6), baseline standardised mean differences (Table 4), the power
and sensitivity analyses (Section 3.8: post-hoc power, stage-IV exclusion,
cause-specific/competing-risk analysis, Ki67 model) and the E-values.

> Note: `reproducible_analysis.py` reads `data/dataset_deidentified.csv`. If run
> from `code/`, adjust the path (`../data/dataset_deidentified.csv`) or run from
> the repository root.

## Data

`dataset_deidentified.csv` contains one row per patient with survival times,
event indicators and model covariates only — no identifiers, names or dates.
Covariate values recovered from source documents (grade, Ki67, age,
radiotherapy) are already integrated (see manuscript Supplementary Figure S1).

## Citation

Fernández Alba JJ, et al. Progressive Cox models expose confounding by
indication in neoadjuvant chemotherapy: an eight-year single-centre breast
cancer cohort. *Clinical and Translational Oncology* (under review).

## License

Code released under the MIT License (see `LICENSE`).
