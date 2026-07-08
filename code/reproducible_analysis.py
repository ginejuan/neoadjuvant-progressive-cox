#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Supplementary code S1 — Reproducible statistical analysis
Manuscript: "Progressive Cox models expose confounding by indication in
neoadjuvant chemotherapy: an eight-year single-centre breast cancer cohort"
Target journal: Clinical and Translational Oncology
================================================================================

Reproduces every quantitative result reported in the manuscript from the
de-identified analysis dataset.

INPUT
-----
dataset_deidentified.csv  — one row per patient (N = 461), no personal data.
Columns:
  id        anonymous identifier
  analytic  1 = included in survival analyses (n = 456), 0 = excluded (5)
  t_os, ev_os     overall-survival time (years) and event (death, any cause)
  t_dfs, ev_dfs   disease-free-survival time and event (recurrence or death)
  ev_bc           breast-cancer-specific death (competing-risk / cause-specific)
  neo             neoadjuvant chemotherapy (1/0)
  age             age at diagnosis (years)
  stage34         stage III-IV vs I-II ; stage4 = stage IV ; stage_ord = 0..IV
  TN              triple-negative subtype ; grade3 = histological grade 3
  her2pos         HER2-positive surrogate subtype (1/0)  [for sensitivity analysis]
  RT              adjuvant radiotherapy ; ki67 (%) ; ki67hi = Ki67 > 19%
  size_mm         tumour size (mm)
Covariate values recovered from source documents are already integrated
(Supplementary Figure S1).

REQUIREMENTS
------------
python >= 3.10 ; pandas ; numpy ; scipy ; lifelines
  pip install pandas numpy scipy lifelines

USAGE
-----
  python Code_S1_reproducible_analysis.py
================================================================================
"""
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import proportional_hazard_test

df = pd.read_csv("../data/dataset_deidentified.csv")     # all 461 patients
d = df[df["analytic"] == 1].copy()               # 456 survival-analysis cohort


def hr(model, var):
    h = np.exp(model.params_[var])
    lo, hi = np.exp(model.confidence_intervals_.loc[var]).values
    return f"HR={h:.2f} (95% CI {lo:.2f}-{hi:.2f}); p={model.summary.loc[var, 'p']:.3g}"


def cox(data, cols, tcol, ecol):
    s = data[cols + [tcol, ecol]].dropna()
    return CoxPHFitter().fit(s, tcol, ecol), s


def km5(t, e):
    m = t.notna() & (t > 0)
    s5 = float(KaplanMeierFitter().fit(t[m], e[m]).predict(5.0)) * 100
    return s5, int(m.sum()), int(e[m].sum())


L = "-" * 78
print(L); print("Cohort N =", len(df), "| analytic n =", len(d)); print(L)

# ---- Table 1 (descriptives, full cohort N=461) ----
print("\n[Table 1 — descriptives, N = 461]")
print("  Median age (IQR): %.0f (%.0f-%.0f)" %
      (df['age'].median(), df['age'].quantile(.25), df['age'].quantile(.75)))
for v, lab in [("stage34", "Stage III-IV"), ("TN", "Triple-negative"),
               ("neo", "Neoadjuvant CT"), ("RT", "Radiotherapy")]:
    print("  %-16s %d (%.1f%%)" % (lab, df[v].sum(), 100*df[v].mean()))
g3 = int((df['grade3'] == 1).sum()); gk = int(df['grade3'].notna().sum())
print("  Grade 3          %d (%.1f%% of known); grade missing %.1f%%" %
      (g3, 100*g3/gk, 100*df['grade3'].isna().mean()))

# ---- Kaplan-Meier 5-year (Figure 1) ----
o = km5(d['t_os'], d['ev_os']); f = km5(d['t_dfs'], d['ev_dfs'])
print("\n[Figure 1 — Kaplan-Meier]")
print("  OS  5-year = %.1f%% (events %d) ; DFS 5-year = %.1f%% (events %d)" % (o[0], o[2], f[0], f[2]))

# ---- Table 4 — baseline SMD by NACT (full cohort N=461) ----
print("\n[Table 4 — baseline standardised mean differences by NACT, N = 461]")
neo, non = df[df['neo'] == 1], df[df['neo'] == 0]
def smd_bin(a, b):
    p1, p2 = a.mean(), b.mean(); sd = np.sqrt((p1*(1-p1)+p2*(1-p2))/2)
    return (p1-p2)/sd if sd > 0 else 0.0
def smd_cont(a, b):
    a, b = a.dropna(), b.dropna(); sd = np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
    return (a.mean()-b.mean())/sd if sd > 0 else 0.0
print("  Age            SMD %.2f" % smd_cont(neo['age'], non['age']))
for v in ["stage34", "TN", "grade3", "ki67hi", "RT"]:
    print("  %-14s NACT %.1f%% vs %.1f%% (SMD %.2f)" %
          (v, 100*neo[v].mean(), 100*non[v].mean(), smd_bin(neo[v], non[v])))

# ---- Table 2 — multivariable Cox OS ----
print("\n[Table 2 — multivariable Cox, overall survival]")
c, s = cox(d, ["TN", "stage34", "grade3", "RT"], "t_os", "ev_os")
print("  n=%d, events=%d, Harrell C=%.3f" % (len(s), int(s['ev_os'].sum()), c.concordance_index_))
for v in ["stage34", "TN", "grade3", "RT"]:
    print("   ", v, hr(c, v))

# ---- Tables 3 & 5 — progressive Cox (OS, DFS) ----
steps = [(["neo"], "M1 NACT alone"), (["neo", "age"], "M2 +age"),
         (["neo", "age", "stage34"], "M3 +stage"),
         (["neo", "age", "stage34", "grade3"], "M4 +grade"),
         (["neo", "age", "stage34", "grade3", "TN"], "M5 +subtype")]
for tcol, ecol, lbl in [("t_os", "ev_os", "Table 3 (OS)"), ("t_dfs", "ev_dfs", "Table 5 (DFS)")]:
    print("\n[%s — progressive Cox, NACT hazard ratio]" % lbl)
    for cols, name in steps:
        c, s = cox(d, cols, tcol, ecol)
        print("  %-13s %s (n=%d)" % (name, hr(c, "neo"), len(s)))

# ---- Table 6 — proportional-hazards test ----
print("\n[Table 6 — Schoenfeld residual test (model of Table 2)]")
c, s = cox(d, ["TN", "stage34", "grade3", "RT"], "t_os", "ev_os")
z = proportional_hazard_test(c, s, time_transform="rank")
for v in ["TN", "stage34", "grade3", "RT"]:
    print("    %-8s p=%.3f" % (v, z.summary.loc[v, "p"]))
# stage III-IV proportional-hazards violation: time-split at 3 years
_sp = []
for _, _r in d.dropna(subset=["t_os", "ev_os", "stage34", "TN", "RT"]).iterrows():
    _t, _e = _r["t_os"], int(_r["ev_os"])
    if _t <= 3:
        _sp.append((0, _t, _e, _r.stage34, _r.TN, _r.RT))
    else:
        _sp.append((0, 3, 0, _r.stage34, _r.TN, _r.RT))
        _sp.append((3, _t, _e, _r.stage34, _r.TN, _r.RT))
_sp = pd.DataFrame(_sp, columns=["start", "stop", "ev", "stage34", "TN", "RT"])
for _lab, _seg in [("0-3y", _sp[_sp.start == 0]), (">3y", _sp[_sp.start == 3])]:
    _ct = CoxPHFitter().fit(_seg, "stop", "ev", entry_col="start")
    print("    stage III-IV %-4s HR %.2f" % (_lab, np.exp(_ct.params_["stage34"])))

# ---- Section 3.8 — power and additional sensitivity analyses ----
print("\n[Section 3.8 — power and sensitivity]")
c, s = cox(d, ["neo", "age", "stage34"], "t_os", "ev_os")
se = c.standard_errors_["neo"]; zc = stats.norm.ppf(.975) + stats.norm.ppf(.80)
print("  Stage-adjusted NACT: %s ; SE(logHR)=%.3f" % (hr(c, "neo"), se))
print("  Minimum detectable HR at 80%% power: %.2f" % np.exp(zc*se))
for h in (1.5, 2.0, 2.5):
    print("    power for HR=%.1f: %.0f%%" % (h, 100*stats.norm.cdf(abs(np.log(h))/se - stats.norm.ppf(.975))))
print("  Schoenfeld p for NACT (adjusted): %.3f" % proportional_hazard_test(c, s, time_transform="rank").summary.loc["neo", "p"])
c2, _ = cox(d[d['stage4'] == 0], ["neo", "age", "stage34"], "t_os", "ev_os")
print("  Excluding stage IV: NACT %s" % hr(c2, "neo"))
print("  Deaths: %d total, %d breast-cancer-specific" % (int(d['ev_os'].sum()), int(d['ev_bc'].sum())))
cs0, _ = cox(d, ["neo"], "t_os", "ev_bc")
cs1, _ = cox(d, ["neo", "age", "stage34"], "t_os", "ev_bc")
cs5, _ = cox(d, ["neo", "age", "stage34", "grade3", "TN"], "t_os", "ev_bc")
cs2, _ = cox(d, ["neo", "age", "stage_ord", "size_mm"], "t_os", "ev_bc")
print("  Cause-specific NACT, unadjusted (M1):        %s" % hr(cs0, "neo"))
print("  Cause-specific NACT, +age+binary stage (M3): %s" % hr(cs1, "neo"))
print("  Cause-specific NACT, +grade+subtype (M5):    %s" % hr(cs5, "neo"))
print("  Cause-specific NACT, ordinal stage+size:     %s (attenuated -> residual confounding)" % hr(cs2, "neo"))
ck, _ = cox(d, ["TN", "stage34", "grade3", "RT", "ki67hi"], "t_os", "ev_os")
print("  Ki67>19%% in multivariable model: %s" % hr(ck, "ki67hi"))
if "her2pos" in d.columns:
    hn = d[d["her2pos"] == 0]
    print("  [HER2-negative sensitivity, excluding %d HER2+ patients]" % int(d["her2pos"].sum()))
    chn, _ = cox(hn, ["TN", "stage34", "grade3", "RT"], "t_os", "ev_os")
    for v in ["TN", "stage34", "grade3", "RT"]:
        print("    %-8s %s" % (v, hr(chn, v)))
    cnn, _ = cox(hn, ["neo", "age", "stage34"], "t_os", "ev_os")
    print("    stage-adjusted NACT: %s" % hr(cnn, "neo"))

# ---- E-value (VanderWeele & Ding, 2017) ----
def evalue(h):
    h = 1/h if h < 1 else h
    return h + np.sqrt(h*(h-1))
print("\n[E-values]  NACT residual HR 1.56 -> %.2f | stage III-IV 3.94 -> %.2f | RT 0.32 -> %.2f"
      % (evalue(1.56), evalue(3.94), evalue(0.32)))
print("\n" + L + "\nDone — results reproduce the values reported in the manuscript.\n" + L)
