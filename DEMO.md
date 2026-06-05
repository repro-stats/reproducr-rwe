# DEMO — reproducr-rwe walkthrough

This document walks through the complete `reproducr` pipeline applied to a
Real World Evidence study comparing two treatments for cardiovascular outcomes
in patients with type 2 diabetes.

---

## Study design

**Research question:** Does treatment A reduce the risk of MACE compared to
treatment B in type 2 diabetes patients, as observed in routine clinical
practice?

**Dataset:** Simulated EHR cohort — 6,000 patients, new-user active-comparator
design, 36 months of follow-up.

**Confounders addressed:** Age, sex, BMI, HbA1c, eGFR, cardiovascular history,
smoking, statin use, index year.

**Methods:** Propensity score matching (primary), IPTW (sensitivity analysis),
Cox proportional hazards, Kaplan-Meier.

---

## Step 1 — Audit the analysis script

```r
library(reproducr)

report <- audit_script("analysis.R", renv = TRUE)
print(report)
```

```
-- reproducr audit report [2026-06-04 09:14] --

  Files scanned:     1
  Packages found:    5
  Calls detected:    31
  R version:         4.4.2
  Platform:          aarch64-apple-darwin20
  Versions from:     renv.lock

  Next step: risks <- risk_score(report)
```

The audit detects 31 qualified `pkg::fn` calls across 5 packages — all
resolved from `renv.lock` for stable version reporting.

---

## Step 2 — Score for risk

```r
risks <- risk_score(report)
print(risks)
```

```
-- reproducr risk score --

  HIGH:      0
  MEDIUM:    0
  LOW:       4

[LOW]  base::sort    locale_check — sort() output is locale-sensitive
[LOW]  base::format  locale_check — format() output is locale-sensitive
...
```

No high or medium risks — the locked environment and explicit `set.seed()`
throughout the pipeline produce a clean audit. The LOW locale flags are
informational — the analysis uses `Sys.setlocale("LC_COLLATE", "C")` for
any sorting operations to ensure cross-platform consistency.

---

## Step 3 — Run the analysis

```r
source("analysis.R")
```

```
Cohort: 6000 patients (2847 treated, 3153 control)

--- Table 1: Pre-matching baseline characteristics ---
...

Matching: 2731 treated matched to 2731 control

--- Balance after matching ---
Max SMD after matching: 0.047 (threshold: 0.10)

--- Primary results: Hazard ratios for MACE ---
Unadjusted:    HR = 0.761 (95% CI: 0.694-0.835, p < 0.001)
Multivariable: HR = 0.812 (95% CI: 0.739-0.892, p < 0.001)
PS-matched:    HR = 0.798 (95% CI: 0.724-0.880, p < 0.001)
IPTW (ATE):    HR = 0.805 (95% CI: 0.733-0.885, p < 0.001)

12-month event-free survival: control 82.4%, treatment 87.1%
36-month event-free survival: control 61.2%, treatment 68.9%

24 outputs ready for certification.
```

Consistent treatment effect across all four estimators — HR approximately
0.80, indicating a ~20% relative reduction in MACE risk with treatment A.

---

## Step 4 — Certify the outputs

```r
certify(
  outputs = OUTPUTS,
  tag     = "submission-v1",
  script  = "analysis.R"
)
```

```
reproducr: certified 24 output(s) [2026-06-04] under tag 'submission-v1'
```

24 outputs hashed and stored in `.reproducr.rds`.

---

## Step 5 — Check for drift

Six months later, after a package upgrade:

```r
source("analysis.R")

check_drift(OUTPUTS, against = "submission-v1")
```

```
-- reproducr drift check vs 'submission-v1' --

  Verdict  : ALL OUTPUTS MATCH
  OK       : 24
  Drifted  : 0
```

All 24 outputs match — the locked `renv` environment and explicit seeds
guarantee numerical reproducibility across runs.

---

## Step 6 — Generate the academic report

```r
repro_report(
  report, risks,
  format      = "md",
  style       = "academic",
  output_file = "reproducibility_report.md"
)
```

The academic report produces a ready-to-paste methods paragraph:

> *All analyses were conducted in R (version 4.4.2) on Ubuntu 24.04.
> Propensity score matching was performed using MatchIt (v4.5.5).
> Inverse probability of treatment weighting used WeightIt (v1.3.0).
> Survival analyses used the survival package (v3.7-0).
> The analysis environment was managed with renv (v1.1.5).
> Reproducibility auditing (reproducr v0.1.1) identified no high- or
> medium-severity risks. The full audit report and certification records
> are available as supplementary materials.*

---

## Step 7 — Badge the README

```r
repro_badge(report, risks, output = "README")
```

The README badge updates automatically on every CI run, reflecting the current
reproducibility status of the analysis.

---

## Why this matters for RWE

RWE analyses are particularly vulnerable to silent breaking changes because:

1. **Long timelines** — study design, analysis, and publication may span years
2. **Active packages** — `MatchIt`, `WeightIt`, `cobalt` are under active
   development with frequent releases
3. **Regulatory scrutiny** — HTA submissions and regulatory RWE packages
   require documented, reproducible methodology
4. **Replication** — external validation requires exact numerical reproducibility

`reproducr` makes these risks visible and provides a verifiable audit trail
from analysis to submission.
