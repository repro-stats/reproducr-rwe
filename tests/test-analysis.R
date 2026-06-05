# tests/test-analysis.R
# Structural tests for the RWE analysis pipeline
# These verify the analysis produces the expected outputs without asserting
# specific numerical values (which vary by platform and package version).

source("analysis.R", local = TRUE)

stopifnot(
  "Cohort size is correct"           = OUTPUTS$n_total == 6000L,
  "Treated and control sum to total" = OUTPUTS$n_treated + OUTPUTS$n_control == 6000L,
  "Event rate is between 0 and 1"    = OUTPUTS$event_rate > 0 & OUTPUTS$event_rate < 1,
  "Max SMD is below 0.15"            = OUTPUTS$smd_max_after < 0.15,
  "Primary HR is positive"           = OUTPUTS$hr_matched > 0,
  "CI lower < HR < CI upper"         = OUTPUTS$ci_matched_lo < OUTPUTS$hr_matched &
                                       OUTPUTS$hr_matched < OUTPUTS$ci_matched_hi,
  "p-value is between 0 and 1"       = OUTPUTS$p_matched >= 0 & OUTPUTS$p_matched <= 1,
  "IPTW HR is positive"              = OUTPUTS$hr_iptw > 0,
  "KM estimates are between 0 and 1" = all(c(OUTPUTS$km_12_ctrl, OUTPUTS$km_12_treat,
                                              OUTPUTS$km_24_ctrl, OUTPUTS$km_24_treat)
                                            >= 0 & c(OUTPUTS$km_12_ctrl, OUTPUTS$km_12_treat,
                                                     OUTPUTS$km_24_ctrl, OUTPUTS$km_24_treat)
                                            <= 1),
  "24 outputs produced"              = length(OUTPUTS) == 24L
)

cat("All tests passed.\n")