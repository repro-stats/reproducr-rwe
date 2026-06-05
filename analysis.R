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
library(WeightIt)
library(tableone)

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
  -7.5 +
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

# ---- 2. Baseline characteristics (Table 1) ----------------------------------

vars       <- c("age", "female", "bmi", "hba1c", "egfr",
                "cv_history", "smoking", "statin_use")
factor_vars <- c("female", "cv_history", "smoking", "statin_use")

tab1_pre <- tableone::CreateTableOne(
  vars       = vars,
  strata     = "treat",
  data       = cohort,
  factorVars = factor_vars
)

cat("\n--- Table 1: Pre-matching baseline characteristics ---\n")
print(tab1_pre, smd = TRUE)

# ---- 3. Propensity score estimation -----------------------------------------

ps_formula <- treat ~ age + female + bmi + hba1c + egfr +
              cv_history + smoking + statin_use + index_year

ps_model <- stats::glm(ps_formula, data = cohort, family = stats::binomial())
cohort$ps <- stats::predict(ps_model, type = "response")

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
  stats      = c("mean.diffs", "variance.ratios"),
  thresholds = c(m = 0.1)
)

cat("\n--- Balance after matching ---\n")
print(bal)

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

# ---- 7. IPTW sensitivity analysis (WeightIt) --------------------------------

wt_out <- WeightIt::weightit(
  formula  = ps_formula,
  data     = cohort,
  method   = "ps",
  estimand = "ATE"
)

cohort$iptw <- wt_out$weights

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

km_summary_6  <- base::summary(km_fit, times = 6)
km_summary_12 <- base::summary(km_fit, times = 12)

.safe_surv <- function(km_sum, idx) {
  s <- km_sum$surv
  if (is.null(s) || length(s) < idx) NA_real_ else s[[idx]]
}

km_6_ctrl   <- .safe_surv(km_summary_6,  1L)
km_6_treat  <- .safe_surv(km_summary_6,  2L)
km_12_ctrl  <- .safe_surv(km_summary_12, 1L)
km_12_treat <- .safe_surv(km_summary_12, 2L)

cat(sprintf("\n6-month event-free survival:  control %.1f%%, treatment %.1f%%\n",
            km_6_ctrl * 100, km_6_treat * 100))
cat(sprintf("12-month event-free survival: control %.1f%%, treatment %.1f%%\n",
            km_12_ctrl * 100, km_12_treat * 100))

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

  # IPTW sensitivity (WeightIt)
  hr_iptw         = round(hr_iptw, 4),
  ci_iptw_lo      = round(ci_iptw[1], 4),
  ci_iptw_hi      = round(ci_iptw[2], 4),
  p_iptw          = round(p_iptw, 6),

  # KM landmarks
  km_6_ctrl       = round(km_6_ctrl,   4),
  km_6_treat      = round(km_6_treat,  4),
  km_12_ctrl      = round(km_12_ctrl,  4),
  km_12_treat     = round(km_12_treat, 4)
)

cat(sprintf("\n%d outputs ready for certification.\n", length(OUTPUTS)))