# reproducr-rwe <a href="https://repro-stats.github.io/reproducr/"><img src="https://raw.githubusercontent.com/repro-stats/reproducr/main/man/figures/logo.svg" align="right" height="120" alt="reproducr website" /></a>

<!-- badges: start -->
[![reproducibility](https://img.shields.io/badge/reproducibility-at%20risk-red)](https://repro-stats.github.io/reproducr/)
[![R-CMD-check](https://github.com/repro-stats/reproducr-rwe/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/repro-stats/reproducr-rwe/actions/workflows/R-CMD-check.yml)
<!-- badges: end -->

> End-to-end `reproducr` pipeline for a simulated Real World Evidence study —
> propensity score matching, IPTW sensitivity analysis, `renv` environment
> locking, academic-style audit report.

| | |
|---|---|
| **Domain** | Real World Evidence / pharmacoepidemiology |
| **Dataset** | Simulated EHR cohort (n = 6,000 patients) |
| **Analysis** | PS matching, IPTW, Cox PH, Kaplan-Meier |
| **Environment** | renv — locked environment |
| **Report style** | academic |
| **Outputs certified** | 22 |
| **Audience** | Epidemiologists, HEOR analysts, pharmacoepidemiologists |

See the full walkthrough: [DEMO.md](DEMO.md)

---

## Study overview

**Research question:** Does treatment A reduce the risk of major adverse
cardiovascular events (MACE) compared to treatment B in patients with type 2
diabetes, as observed in routine clinical practice?

**Design:** Retrospective cohort study using simulated electronic health
records. New-user, active-comparator design with 36 months of follow-up.

**Key results (simulated data):**

| Analysis | HR | 95% CI | p-value |
|---|---|---|---|
| Unadjusted | — | — | — |
| Multivariable adjusted | — | — | — |
| PS-matched (primary) | — | — | — |
| IPTW (sensitivity) | — | — | — |

*Results populated by CI on each run — see `.reproducr.rds` for certified values.*

---

## What this demonstrates

**Why RWE analyses are particularly vulnerable to silent breaking changes:**

RWE pipelines are long, involve many packages at active development stages,
and are often run years apart — at study design, during analysis, and at
publication. Each gap is an opportunity for package behaviour to shift silently.

Specific risks in this pipeline:

- **`MatchIt`** — default matching method and caliper behaviour have evolved
  across versions; matched pairs can differ silently
- **`WeightIt`** — estimand defaults and weight trimming behaviour have changed
- **`cobalt`** — balance statistic calculations updated in recent versions
- **`survival`** — `survfit()` output structure changed between major versions

**Tier 1 — Scan & score**

`audit_script()` scans `analysis.R`, resolves versions from `renv.lock`, and
`risk_score()` flags any known breaking changes, missing seeds, or
locale-sensitive operations across all packages in the pipeline.

**Tier 2 — Baseline & drift**

`certify()` hashes 22 key outputs — cohort counts, balance statistics,
hazard ratios with confidence intervals, p-values, and Kaplan-Meier landmarks.
`check_drift()` compares against the previous certified run on every push.

**Tier 3 — Report & export**

`repro_report()` generates an academic-style methods paragraph suitable for
inclusion in an HEOR manuscript or HTA submission. `repro_badge()` updates the
README badge automatically.

---

## Reproducibility design choices

**`set.seed(2026L)` at the top of `analysis.R`** — every stochastic call in
the pipeline (cohort simulation, matching) is covered by a single seed set at
the script entry point.

**`renv` for environment locking** — `renv.lock` records the exact version of
every package. `reproducr` then verifies that those locked versions do not
contain known silent breaking changes.

**Qualified calls throughout** — all function calls use `pkg::fn` notation.
This makes every dependency explicit and allows `audit_script()` to detect
them. It also makes the code self-documenting.

---

## Running locally

```r
# Clone the repo, then:
renv::restore()        # restore the locked environment
source("analysis.R")

library(reproducr)
report <- audit_script("analysis.R", renv = TRUE)
risks  <- risk_score(report)
print(risks)
```

---

## CI/CD

Two workflows run on every push and weekly:

| Workflow | Purpose |
|---|---|
| `R-CMD-check.yml` | Restore renv, run structural tests |
| `reproducr-audit.yml` | Audit, certify, detect drift, update badge |

---

## Part of the reproducr gallery

| Example | Domain | renv | Report style |
|---|---|---|---|
| [reproducr-ecology](https://github.com/repro-stats/reproducr-ecology) | Ecology / penguins | No | minimal |
| [reproducr-clinical](https://github.com/repro-stats/reproducr-clinical) | Clinical trials / oncology | Yes | pharma |
| **reproducr-rwe** (this repo) | Real world evidence | Yes | academic |
