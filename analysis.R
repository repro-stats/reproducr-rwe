# =============================================================================
# Real World Evidence Analysis — reproducr-rwe
#
# Study:   Comparative effectiveness of treatment A vs treatment B
#          in patients with type 2 diabetes from simulated EHR data
#
# Design:  Retrospective cohort study using propensity score matching
#          Primary endpoint: time to first major adverse cardiovascular
#          event (MACE) over 36 months of follow-up
#
# Author:  Ndoh Penn
# Repo:    https://github.com/repro-stats/reproducr-rwe
# =============================================================================

set.seed(2026L)

# ---- 0. Dependencies --------------------------------------------------------

library(survival)
library(MatchIt)
library(cobalt)

# ---- 1. Simulate EHR cohort -------------------------------------------------

n <- 6000L

# Patient characteristics
age        <- stats::rnorm(n, mean = 62, sd = 11)
female     <- stats::rbinom(n, 1L, 0.48)
bmi        <- stats::rnorm(n, mean = 29.5, sd = 5.2)
hba1c      <- stats::rnorm(n, mean = 7.8, sd = 1.4)
egfr       <- stats::rnorm(n, mean = 74, sd = 18)
cv_history <- stats::rbinom(n, 1L, 0.22)
smoking    <- stats::rbinom(n, 1L, 0.18)
statin_use <- stats::rbinom(n, 1L, 0.61)
index_year <- base::sample(2018L:2022L, n, replace = TRUE)

# Treatment assignment — confounded by age, HbA1c, CV history
log_odds_treat <- -3.8 +
  0.025 * age +
  0.18  * hba1c +
  0.52  * cv_history +
  -0.14 * (egfr / 10) +
  0.08  * bmi

ps_true  <- 1 / (1 + exp(-log_odds_treat))
treat    <- stats::rbinom(n, 1L, ps_true)

# Time-to-MACE outcome (Weibull)
baseline_hazard <- exp(
  -4.8 +
  0.04  * age +
  0.35  * cv_history +
  0.18  * hba1c +
  -0.22 * treat +          # true treatment effect: HR ~ 0.80
  -0.12 * statin_use +
  0.09  * smoking
)

time_to_event <- stats::rweibull(n, shape = 1.4, scale = 1 / baseline_hazard)
follow_up     <- stats::runif(n, min = 0.5, max = 36)
event_time    <- pmin(time_to_event, follow_up)
event         <- as.integer(time_to_event <= follow_up)

# Assemble cohort
cohort <- data.frame(
  patient_id  = seq_len(n),
  treat       = treat,
  age         = round(age, 1),
  female      = female,
  bmi         = round(bmi, 1),
  hba1c       = round(hba1c, 1),
  egfr        = round(pmax(egfr, 15), 1),
  cv_history  = cv_history,
  smoking     = smoking,
  statin_use  = statin_use,
  index_year  = index_year,
  time        = round(event_time, 2),
  event       = event,
  stringsAsFactors = FALSE
)

cat(sprintf("Cohort: %d patients (%d treated, %d control)\n",
            nrow(cohort), sum(cohort$treat), sum(!cohort$treat)))

# ---- 2. Baseline characteristics --------------------------------------------

vars <- c("age", "female", "bmi", "hba1c", "egfr",
          "cv_history", "smoking", "statin_use")

# Table 1 — baseline characteristics by treatment group (base R)
cat("\n--- Table 1: Pre-matching baseline characteristics ---\n")
tab1 <- do.call(rbind, lapply(vars, function(v) {
  x0 <- cohort[[v]][cohort$treat == 0]
  x1 <- cohort[[v]][cohort$treat == 1]
  if (length(unique(cohort[[v]])) == 2L) {
    # Binary variable — report proportion
    data.frame(
      variable = v,
      control  = sprintf("%.1f%%", mean(x0) * 100),
      treated  = sprintf("%.1f%%", mean(x1) * 100),
      smd      = sprintf("%.3f",
                   (mean(x1) - mean(x0)) /
                   sqrt((stats::var(x1) + stats::var(x0)) / 2))
    )
  } else {
    # Continuous variable — report mean (SD)
    data.frame(
      variable = v,
      control  = sprintf("%.1f (%.1f)", mean(x0), stats::sd(x0)),
      treated  = sprintf("%.1f (%.1f)", mean(x1), stats::sd(x1)),
      smd      = sprintf("%.3f",
                   (mean(x1) - mean(x0)) /
                   sqrt((stats::var(x1) + stats::var(x0)) / 2))
    )
  }
}))
print(tab1, row.names = FALSE)

# ---- 3. Propensity score estimation -----------------------------------------

ps_formula <- treat ~ age + female + bmi + hba1c + egfr +
              cv_history + smoking + statin_use + index_year

ps_model <- stats::glm(ps_formula, data = cohort, family = stats::binomial())

cohort$ps <- stats::predict(ps_model, type = "response")

# C-statistic (concordance) for PS model
ps_concordance <- survival::concordance(
  stats::glm(ps_formula, data = cohort, family = stats::binomial())
)$concordance
cat(sprintf("\nPS model C-statistic: %.3f\n", ps_concordance))

# ---- 4. Propensity score matching -------------------------------------------

match_out <- MatchIt::matchit(
  formula = ps_formula,
  data    = cohort,
  method  = "nearest",
  ratio   = 1L,
  replace = FALSE,
  caliper = 0.2 * stats::sd(log(cohort$ps / (1 - cohort$ps)))
)

cat(sprintf("\nMatching: %d treated matched to %d control\n",
            sum(match_out$weights[cohort$treat == 1] > 0),
            sum(match_out$weights[cohort$treat == 0] > 0)))

matched_data <- MatchIt::match.data(match_out)

# ---- 5. Balance assessment --------------------------------------------------

bal <- cobalt::bal.tab(
  match_out,
  stats    = c("mean.diffs", "variance.ratios"),
  thresholds = c(m = 0.1)
)

cat("\n--- Balance after matching ---\n")
print(bal)

# Max standardised mean difference after matching
smd_after <- max(abs(bal$Balance$Diff.Adj), na.rm = TRUE)
cat(sprintf("Max SMD after matching: %.3f (threshold: 0.10)\n", smd_after))

# ---- 6. Outcome analysis — Cox proportional hazards -------------------------

cox_unadj <- survival::coxph(
  survival::Surv(time, event) ~ treat,
  data = cohort
)

cox_adj <- survival::coxph(
  survival::Surv(time, event) ~ treat + age + female + bmi +
    hba1c + egfr + cv_history + smoking + statin_use,
  data = cohort
)

cox_matched <- survival::coxph(
  survival::Surv(time, event) ~ treat,
  data    = matched_data,
  weights = matched_data$weights
)

# Extract primary results
hr_unadj   <- exp(stats::coef(cox_unadj)["treat"])
ci_unadj   <- exp(stats::confint(cox_unadj)["treat", ])
p_unadj    <- summary(cox_unadj)$coefficients["treat", "Pr(>|z|)"]

hr_adj     <- exp(stats::coef(cox_adj)["treat"])
ci_adj     <- exp(stats::confint(cox_adj)["treat", ])
p_adj      <- summary(cox_adj)$coefficients["treat", "Pr(>|z|)"]

hr_matched <- exp(stats::coef(cox_matched)["treat"])
ci_matched <- exp(stats::confint(cox_matched)["treat", ])
p_matched  <- summary(cox_matched)$coefficients["treat", "Pr(>|z|)"]

cat("\n--- Primary results: Hazard ratios for MACE ---\n")
cat(sprintf("Unadjusted:      HR = %.3f (95%% CI: %.3f-%.3f, p = %.4f)\n",
            hr_unadj, ci_unadj[1], ci_unadj[2], p_unadj))
cat(sprintf("Multivariable:   HR = %.3f (95%% CI: %.3f-%.3f, p = %.4f)\n",
            hr_adj, ci_adj[1], ci_adj[2], p_adj))
cat(sprintf("PS-matched:      HR = %.3f (95%% CI: %.3f-%.3f, p = %.4f)\n",
            hr_matched, ci_matched[1], ci_matched[2], p_matched))

# ---- 7. IPTW sensitivity analysis (manual — base R only) -------------------
#
# Stabilised IPTW weights for Average Treatment Effect (ATE):
#   treated:   w = P(A=1) / PS
#   control:   w = P(A=0) / (1 - PS)
# Weights trimmed at 99th percentile to reduce influence of extreme values.

p_treat     <- mean(cohort$treat)
iptw_treat  <- p_treat / cohort$ps
iptw_ctrl   <- (1 - p_treat) / (1 - cohort$ps)
cohort$iptw <- ifelse(cohort$treat == 1L, iptw_treat, iptw_ctrl)

# Trim at 99th percentile
w99         <- stats::quantile(cohort$iptw, 0.99)
cohort$iptw <- pmin(cohort$iptw, w99)

cat(sprintf("IPTW weight range: %.2f - %.2f (trimmed at 99th pctile: %.2f)\n",
            min(cohort$iptw), max(cohort$iptw), w99))

cox_iptw <- survival::coxph(
  survival::Surv(time, event) ~ treat,
  data    = cohort,
  weights = cohort$iptw
)

hr_iptw <- exp(stats::coef(cox_iptw)["treat"])
ci_iptw <- exp(stats::confint(cox_iptw)["treat", ])
p_iptw  <- summary(cox_iptw)$coefficients["treat", "Pr(>|z|)"]

cat(sprintf("IPTW (ATE):      HR = %.3f (95%% CI: %.3f-%.3f, p = %.4f)\n",
            hr_iptw, ci_iptw[1], ci_iptw[2], p_iptw))

# ---- 8. Kaplan-Meier curves (matched cohort) --------------------------------

km_fit <- survival::survfit(
  survival::Surv(time, event) ~ treat,
  data = matched_data
)

# 12- and 24-month event-free survival
km_summary_12 <- base::summary(km_fit, times = 12)
km_summary_24 <- base::summary(km_fit, times = 24)

# Safe extraction — survfit may not reach all time points in small cohorts
.safe_surv <- function(km_sum, idx) {
  s <- km_sum$surv
  if (is.null(s) || length(s) < idx) NA_real_ else s[[idx]]
}

km_12_ctrl  <- .safe_surv(km_summary_12, 1L)
km_12_treat <- .safe_surv(km_summary_12, 2L)
km_24_ctrl  <- .safe_surv(km_summary_24, 1L)
km_24_treat <- .safe_surv(km_summary_24, 2L)

cat(sprintf("\n12-month event-free survival: control %.1f%%, treatment %.1f%%\n",
            km_12_ctrl * 100, km_12_treat * 100))
cat(sprintf("24-month event-free survival: control %.1f%%, treatment %.1f%%\n",
            km_24_ctrl * 100, km_24_treat * 100))

# ---- 9. Collect outputs for certification -----------------------------------

OUTPUTS <- list(
  # Cohort
  n_total         = nrow(cohort),
  n_treated       = sum(cohort$treat),
  n_control       = sum(!cohort$treat),
  event_rate      = round(mean(cohort$event), 4),

  # Balance
  smd_max_after   = round(smd_after, 4),
  n_matched_pairs = sum(match_out$weights[cohort$treat == 1] > 0),

  # Primary — PS-matched
  hr_matched      = round(hr_matched, 4),
  ci_matched_lo   = round(ci_matched[1], 4),
  ci_matched_hi   = round(ci_matched[2], 4),
  p_matched       = round(p_matched, 6),

  # Multivariable adjusted
  hr_adj          = round(hr_adj, 4),
  ci_adj_lo       = round(ci_adj[1], 4),
  ci_adj_hi       = round(ci_adj[2], 4),
  p_adj           = round(p_adj, 6),

  # IPTW sensitivity
  hr_iptw         = round(hr_iptw, 4),
  ci_iptw_lo      = round(ci_iptw[1], 4),
  ci_iptw_hi      = round(ci_iptw[2], 4),
  p_iptw          = round(p_iptw, 6),

  # KM landmarks
  km_12_ctrl      = round(km_12_ctrl,  4),
  km_12_treat     = round(km_12_treat, 4),
  km_24_ctrl      = round(km_24_ctrl,  4),
  km_24_treat     = round(km_24_treat, 4)
)

cat(sprintf("\n%d outputs ready for certification.\n", length(OUTPUTS)))