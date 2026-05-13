# =============================================================================
# Four-wave (JACSIS 2022 / 2023 / 2024 / 2025) prospective study of
# generative-AI initiation and subjective wellbeing.
#
# Analysis script for:
#   Nakagomi, A., Akutsu, Y., Yasuoka, M. & Tabuchi, T. (2026). 
#   Title: Generative AI use and health and wellbeing in Japanese adults: a longitudinal outcome-wide difference-in-differences analysis.
#   [Journal]. DOI: [will be added on acceptance].
#
# Author:     Atsushi Nakagomi (a.nakagomi@chiba-u.jp)
# License:    MIT (see LICENSE file)
#
# DATA AVAILABILITY
# -----------------
# JACSIS data are not publicly available due to ethical restrictions on
# participant privacy. De-identified data may be shared upon reasonable
# request, subject to approval by the JACSIS steering committee and the
# relevant ethics review boards. See manuscript Data Availability section
# for details. The three panel CSV files expected by this script (see paths
# below) must be supplied by approved researchers; column names should match
# those documented in CODEBOOK.md.
#
# REPRODUCIBILITY
# ---------------
# - R version 4.5.3
# - Required packages: readr, dplyr, tidyr, purrr, sandwich, lmtest,
#                       psych, car, ggplot2
# - Random seed: set.seed(20260512) is set below.
#
# USAGE
# -----
#   1. Place the three panel CSV files at the paths defined below
#      (DATA_DIR/PANEL_*). Column names must match CODEBOOK.md.
#   2. From the repository root: `Rscript analysis.R`
#   3. Outputs are written to the directory defined by OUTPUT_DIR.
#
# PIPELINE
# --------
# Four analytic tiers, all on the same 13 wellbeing outcomes x 5 AI-use
# exposures:
#   Main    DiD     (codes 5+6 vs 1; outcome change 2024 -> 2025)
#   Sens 1  DiD     (codes 4+5+6 vs 1; outcome change 2023 -> 2025)
#   Sens 2  ANCOVA  (same cohort as Main; Y_2025 ~ exposure + Y_2024 + cov)
#   Sens 3  ANCOVA  (full sample; Y_2025 ~ exposure + Y_2022 + cov)
# Followed by diagnostics (EFA on the AI-purpose items + Cronbach's alpha +
# VIF) and publication tables (Table 1-3, Supplementary Tables 1A/1B/2/3/4).
# =============================================================================

set.seed(20260512)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_DIR        <- "data"        # directory containing the three panel CSVs
PANEL_2W        <- file.path(DATA_DIR, "jacsis_2022_2025.csv")
PANEL_3W_MAIN   <- file.path(DATA_DIR, "jacsis_2023_2024_2025.csv")
PANEL_3W_SENS1  <- file.path(DATA_DIR, "jacsis_2022_2023_2025.csv")

OUTPUT_DIR      <- "output"      # all outputs go under here
MS_DIR    <- file.path(OUTPUT_DIR, "manuscript")
TF_DIR    <- file.path(MS_DIR,     "tables_figures")
DG_DIR    <- file.path(MS_DIR,     "diagnostics")
TF_MAIN   <- file.path(TF_DIR,     "main")
TF_SENS1  <- file.path(TF_DIR,     "sens1")
TF_SENS2  <- file.path(TF_DIR,     "sens2")
TF_SENS3  <- file.path(TF_DIR,     "sens3")
for (d in c(MS_DIR, TF_DIR, DG_DIR, TF_MAIN, TF_SENS1, TF_SENS2, TF_SENS3))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Aliases used by the per-tier sections below.
INPUT_FILE  <- PANEL_2W        # 2-wave panel (Sens 3)
INPUT_3WAVE <- PANEL_3W_MAIN   # 3-wave panel for Main + Sens 2
INPUT_S1    <- PANEL_3W_SENS1  # 3-wave panel for Sens 1
INPUT_S3    <- PANEL_2W        # 2-wave panel for Sens 3 (= INPUT_FILE)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(sandwich); library(lmtest); library(psych); library(car)
})

cat(sprintf("[setup] OUTPUT_DIR=%s\n", OUTPUT_DIR))
cat(sprintf("[setup] PANEL_2W       = %s\n", PANEL_2W))
cat(sprintf("[setup] PANEL_3W_MAIN  = %s\n", PANEL_3W_MAIN))
cat(sprintf("[setup] PANEL_3W_SENS1 = %s\n", PANEL_3W_SENS1))

# ---------------------------------------------------------------------------
# H E L P E R   F U N C T I O N S
# ---------------------------------------------------------------------------

as_num <- function(x) suppressWarnings(as.numeric(x))

row_mean <- function(data, cols, na_rm = TRUE) {
  sub <- as.data.frame(data[, cols, drop = FALSE])
  m <- matrix(as_num(unlist(sub, use.names = FALSE)),
              nrow = nrow(sub), ncol = length(cols))
  rowMeans(m, na.rm = na_rm)
}

row_sum_tx <- function(data, cols, fn = identity, na_rm = TRUE) {
  sub <- as.data.frame(data[, cols, drop = FALSE])
  m <- matrix(as_num(unlist(sub, use.names = FALSE)),
              nrow = nrow(sub), ncol = length(cols))
  m2 <- fn(m)
  if (!is.matrix(m2)) {
    m2 <- matrix(as.numeric(m2), nrow = nrow(m), ncol = ncol(m))
  }
  rowSums(m2, na.rm = na_rm)
}

chronic_binary <- function(x) {
  v <- as_num(x)
  out <- as.integer(v %in% c(3, 4, 5))
  out[is.na(out)] <- 0L
  out
}

# Robust SE (HC1) coefficient extractor
robust_coef <- function(model, term) {
  vc <- tryCatch(sandwich::vcovHC(model, type = "HC1"), error = function(e) NULL)
  co <- coef(model)
  if (is.null(vc) || !term %in% rownames(vc) || !term %in% names(co)) {
    sm <- summary(model)$coefficients
    if (!term %in% rownames(sm)) {
      return(c(estimate = NA_real_, se = NA_real_, z = NA_real_,
               p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
    }
    est <- sm[term, 1]; se <- sm[term, 2]
    z <- est / se; p <- 2 * pnorm(-abs(z))
    return(c(estimate = est, se = se, z = z, p = p,
             ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se))
  }
  est <- co[[term]]; se <- sqrt(vc[term, term])
  z <- est / se; p <- 2 * pnorm(-abs(z))
  c(estimate = est, se = se, z = z, p = p,
    ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se)
}

# car::vif wrapper
compute_vif_safe <- function(model) {
  tryCatch({
    v <- car::vif(model)
    if (is.matrix(v)) v <- v[, ncol(v)]
    out <- as.numeric(v)
    names(out) <- if (is.null(names(v))) seq_along(out) else names(v)
    out
  }, error = function(e) NA_real_)
}

drop_constant_cols <- function(d) {
  keep <- vapply(d, function(x) length(unique(x[is.finite(x)])) >= 2L, logical(1))
  d[, keep, drop = FALSE]
}

# Outcome-z: append columns based on baseline-wave SD per outcome.
swap_to_oz <- function(d) {
  d$estimate <- d$estimate_oz
  d$se       <- d$se_oz
  d$ci_lo    <- d$ci_lo_oz
  d$ci_hi    <- d$ci_hi_oz
  d
}

add_outcome_z <- function(results_df, sd_pre, outcome_col = "outcome") {
  ids <- results_df[[outcome_col]]
  s   <- sd_pre[ids]
  s[!is.finite(s) | s <= 0] <- NA_real_
  results_df$sd_pre      <- unname(s)
  results_df$estimate_oz <- results_df$estimate / s
  results_df$se_oz       <- results_df$se       / s
  results_df$ci_lo_oz    <- results_df$ci_lo    / s
  results_df$ci_hi_oz    <- results_df$ci_hi    / s
  results_df
}

# Recode functions
k6_recode  <- function(M) pmax(0, pmin(4, 5 - M))
ucla_recode <- function(M) pmax(0, pmin(3, 4 - M))

# E-value — VanderWeele & Ding 2017. For a standardized
# effect size d on a continuous outcome, RR ≈ exp(0.91 × |d|) (Chinn 2000),
# then E = RR + sqrt(RR × (RR − 1)). E-value at CI = E-value evaluated at the
# CI bound closer to the null; if the CI crosses zero, the data are
# consistent with no effect and E-value at CI = 1 (no confounding needed).
.e_value_from_d <- function(d) {
  if (!is.finite(d)) return(NA_real_)
  rr <- exp(0.91 * abs(d))
  rr + sqrt(rr * (rr - 1))
}
add_e_values_oz <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  est <- df$estimate_oz
  lo  <- df$ci_lo_oz
  hi  <- df$ci_hi_oz
  df$e_value <- vapply(est, .e_value_from_d, numeric(1))
  ci_excludes_null <- is.finite(lo) & is.finite(hi) & (sign(lo) == sign(hi)) &
                      lo != 0 & hi != 0
  bound_closer <- ifelse(abs(lo) < abs(hi), lo, hi)
  e_ci <- vapply(bound_closer, .e_value_from_d, numeric(1))
  df$e_value_ci <- ifelse(ci_excludes_null, e_ci, 1.0)
  df
}

# Time-use 2025 re-bin.
# 2025 splits the 4-5h, 6-7h, 8-9h, 10-11h ranges into individual hours and
# adds 4 levels; わからない moves from position 12 → 16. Re-bin so that the
# numeric value carries the same real-world meaning as 2022/2023/2024:
#   1-6 align (none / <30min / 30min / 1h / 2h / 3h)
#   2025 codes 7+8 (4h, 5h)   → 7  ("4-5時間")
#                9+10 (6h, 7h)  → 8  ("6-7時間")
#               11+12 (8h, 9h)  → 9  ("8-9時間")
#               13+14 (10h, 11h)→ 10 ("10-11時間")
#                15 (12h+)       → 11 ("12時間以上")
#                16 (わからない) → 12 ("わからない")
recode_timeuse_2025 <- function(x) {
  ifelse(x %in% 1:6,   x,
  ifelse(x %in% 7:8,   7L,
  ifelse(x %in% 9:10,  8L,
  ifelse(x %in% 11:12, 9L,
  ifelse(x %in% 13:14, 10L,
  ifelse(x == 15,      11L,
  ifelse(x == 16,      12L, NA_integer_)))))))
}

.timeuse_outcomes <- c("d7_sleep", "d7_sitting", "d7_walking")

# Time-use band → hours conversion. After recode_timeuse_2025
# the variable is on the 12-level band scale; band_to_hours converts it to
# numeric hours/day using midpoints for ranges, with わからない → NA.
#   1=なし(0時間)→0; 2=30分未満→0.25; 3=30分程度→0.5; 4=1時間→1; 5=2時間→2;
#   6=3時間→3; 7=4-5時間→4.5; 8=6-7時間→6.5; 9=8-9時間→8.5; 10=10-11時間→10.5;
#   11=12時間以上→12; 12=わからない→NA.
.timeuse_band_to_hours <- c(0, 0.25, 0.5, 1, 2, 3, 4.5, 6.5, 8.5, 10.5, 12, NA_real_)
band_to_hours <- function(x) {
  v <- suppressWarnings(as.integer(x))
  out <- rep(NA_real_, length(v))
  ok <- !is.na(v) & v >= 1L & v <= 12L
  out[ok] <- .timeuse_band_to_hours[v[ok]]
  out
}

apply_timeuse_hours <- function(df, vars = .timeuse_outcomes) {
  for (id in intersect(vars, names(df))) {
    df[[id]] <- band_to_hours(df[[id]])
  }
  df
}

# ---------------------------------------------------------------------------
# Outcome specs (13 — shared across all DiD/ANCOVA tiers)
# Each spec carries the per-wave column maps (c22 / c23 / c24 / c25).
# ---------------------------------------------------------------------------
sens_specs <- list(
  list(id="d1_pm_avg",     kind="continuous", aggregator="mean",
       c22=c("Q78.3_2022","Q78.4_2022"),
       c23=c("Q86.3_2023","Q86.4_2023"),
       c24=c("Q76.3_2024","Q76.4_2024"),
       c25=c("Q74.3_2025","Q74.4_2025")),
  list(id="d1_k6",         kind="sum_recode",
       c22=paste0("Q68.", 1:6, "_2022"),
       c23=paste0("Q75.", 1:6, "_2023"),
       c24=paste0("Q65.", 1:6, "_2024"),
       c25=paste0("Q65.", 1:6, "_2025"),
       recode=k6_recode),
  list(id="d2_hls_avg",    kind="continuous", aggregator="mean",
       c22=c("Q78.1_2022","Q78.2_2022"),
       c23=c("Q86.1_2023","Q86.2_2023"),
       c24=c("Q76.1_2024","Q76.2_2024"),
       c25=c("Q74.1_2025","Q74.2_2025")),
  list(id="d3_mp_avg",     kind="continuous", aggregator="mean",
       c22=c("Q78.5_2022","Q78.14_2022"),
       c23=c("Q86.5_2023","Q86.14_2023"),
       c24=c("Q76.5_2024","Q76.14_2024"),
       c25=c("Q74.5_2025","Q74.14_2025")),
  list(id="d3_ikigai",     kind="continuous",
       c22="Q78.6_2022", c23="Q86.6_2023", c24="Q76.6_2024", c25="Q74.6_2025"),
  list(id="d4_cv_avg",     kind="continuous", aggregator="mean",
       c22=c("Q78.8_2022","Q78.15_2022"),
       c23=c("Q86.8_2023","Q86.15_2023"),
       c24=c("Q76.8_2024","Q76.15_2024"),
       c25=c("Q74.8_2025","Q74.15_2025")),
  list(id="d5_rel_avg",    kind="continuous", aggregator="mean",
       c22=c("Q78.9_2022","Q78.10_2022"),
       c23=c("Q86.9_2023","Q86.10_2023"),
       c24=c("Q76.9_2024","Q76.10_2024"),
       c25=c("Q74.9_2025","Q74.10_2025")),
  list(id="d5_lsns_total", kind="continuous", aggregator="sum",
       c22=paste0("Q18.", 1:6, "_2022"),
       c23=paste0("Q18.", 1:6, "_2023"),
       c24=paste0("Q17.", 1:6, "_2024"),
       c25=paste0("Q20.", 1:6, "_2025")),
  list(id="d5_ucla",       kind="sum_recode",
       c22=paste0("Q68S1.", 1:3, "_2022"),
       c23=paste0("Q76.",   1:3, "_2023"),
       c24=paste0("Q66.",   1:3, "_2024"),
       c25=paste0("Q66.",   1:3, "_2025"),
       recode=ucla_recode),
  list(id="d6_worry_avg",  kind="continuous", aggregator="mean",
       c22=c("Q78.11_2022","Q78.12_2022"),
       c23=c("Q86.11_2023","Q86.12_2023"),
       c24=c("Q76.11_2024","Q76.12_2024"),
       c25=c("Q74.11_2025","Q74.12_2025")),
  list(id="d7_sleep",      kind="continuous",
       c22="Q28.9_2022", c23="Q31.9_2023", c24="Q28.9_2024", c25="Q32.9_2025"),
  list(id="d7_sitting",    kind="continuous",
       c22="Q28.5_2022", c23="Q31.5_2023", c24="Q28.5_2024", c25="Q32.5_2025"),
  list(id="d7_walking",    kind="continuous",
       c22="Q28.6_2022", c23="Q31.6_2023", c24="Q28.6_2024", c25="Q32.6_2025")
)
stopifnot(length(sens_specs) == 13L)

outcome_order <- c(
  "d1_pm_avg","d1_k6",
  "d2_hls_avg",
  "d3_mp_avg","d3_ikigai",
  "d4_cv_avg",
  "d5_rel_avg","d5_lsns_total","d5_ucla",
  "d6_worry_avg",
  "d7_sleep","d7_sitting","d7_walking"
)
outcome_labels <- c(
  d1_pm_avg     = "Physical & mental health",
  d1_k6         = "K6 distress (sum 0-24)",
  d2_hls_avg    = "Happiness & life satisfaction",
  d3_mp_avg     = "Meaning & purpose",
  d3_ikigai     = "Ikigai",
  d4_cv_avg     = "Character & virtue",
  d5_rel_avg    = "Relations",
  d5_lsns_total = "LSNS-6 total (0-30)",
  d5_ucla       = "UCLA-3 loneliness",
  d6_worry_avg  = "Financial & safety worry",
  d7_sleep      = "Sleep hours",
  d7_sitting    = "Sitting time",
  d7_walking    = "Walking time"
)
domain_of <- c(
  rep("1. Mental & Physical Health", 2),
  rep("2. Happiness & Life Satisfaction", 1),
  rep("3. Meaning & Purpose", 2),
  rep("4. Character & Virtue", 1),
  rep("5. Social Wellbeing", 3),
  rep("6. Financial & Material Stability", 1),
  rep("7. Health Behavior", 3)
)
names(domain_of) <- outcome_order

# ---------------------------------------------------------------------------
# Outcome materialiser: build a single Y_<wave> column for one spec on dataset ds.
# ---------------------------------------------------------------------------
build_y_for_wave <- function(ds, spec, wave) {
  cols <- spec[[wave]]
  if (spec$kind == "binary") {
    chronic_binary(ds[[cols]])
  } else if (spec$kind == "sum_recode") {
    row_sum_tx(ds, cols, fn = spec$recode, na_rm = TRUE)
  } else if (spec$kind == "binary_lt") {
    s <- if (!is.null(spec$recode))
           row_sum_tx(ds, cols, fn = spec$recode, na_rm = TRUE)
         else row_sum_tx(ds, cols, na_rm = TRUE)
    as.integer(s < spec$threshold)
  } else if (spec$kind == "binary_ge") {
    s <- if (!is.null(spec$recode))
           row_sum_tx(ds, cols, fn = spec$recode, na_rm = TRUE)
         else row_sum_tx(ds, cols, na_rm = TRUE)
    as.integer(s >= spec$threshold)
  } else {
    if (length(cols) == 1L) {
      as_num(ds[[cols]])
    } else if (!is.null(spec$aggregator) && spec$aggregator == "sum") {
      row_sum_tx(ds, cols, na_rm = TRUE)
    } else {
      row_mean(ds, cols, na_rm = TRUE)
    }
  }
}

# Materialise Y for all specs at one wave key. Returns data.frame.
# Time-use values are returned in HOURS/DAY: for 2025, first
# applies recode_timeuse_2025 (16→12 bands), then band_to_hours (band → hours).
# For 2022/2023/2024 they're already on the 12-band scale → directly band_to_hours.
# わからない (band 12) → NA, so rows naturally drop from time-use regressions
# via the `keep` filter in fit_did_one. Covariate adjustment via fill_zero_pre
# converts those NAs to 0 (acceptable: small misclassification of わからない as
# "0 hours" for OTHER outcomes' regressions; main analysis is on the focal Y).
materialise_outcomes <- function(ds, wave_key) {
  out <- data.frame(row.names = seq_len(nrow(ds)))
  for (sp in sens_specs) {
    v <- build_y_for_wave(ds, sp, wave_key)
    if (sp$id %in% .timeuse_outcomes) {
      if (wave_key == "c25") v <- recode_timeuse_2025(v)
      v <- band_to_hours(v)
    }
    out[[sp$id]] <- v
  }
  out
}

# Pre-outcome NA → 0 (consistent with locked pipeline).
fill_zero_pre <- function(d) {
  for (cn in names(d)) {
    v <- d[[cn]]; v[!is.finite(v)] <- 0; d[[cn]] <- v
  }
  d
}

# ---------------------------------------------------------------------------
# Confounder derivation (parametric on baseline wave suffix).
# Requires `ds` to be available in the calling environment (per-tier scope).
# ---------------------------------------------------------------------------
derive_confounders <- function(ds, suffix) {
  sex_col   <- paste0("SEX_",     suffix)
  age_col   <- paste0("AGE_",     suffix)
  edu_col   <- paste0("Q21.1_",   suffix)   # NB: 2025 education is Q25.1 (not used here).
  # Bug #1 fix: 2023 employment is at Q6.1 (Q5.1_2023 is a
  # vaccine-attitudes question). Per-wave switch:
  emp_col   <- switch(suffix,
                      "2022" = "Q5.1_2022",
                      "2023" = "Q6.1_2023",
                      "2024" = "Q5.1_2024",
                      "2025" = "Q5.1_2025",
                      paste0("Q5.1_", suffix))
  # Bug #2 / #3 fix: income Q-code shifts across waves.
  # 2022 = Q87.1; 2023 = Q90.1 (Q87.1_2023 is ACE); 2024 = Q80.1 (no fallback
  # needed — Q80.1_2024 is the actual 2024 income question); 2025 = Q77.1.
  inc_col   <- switch(suffix,
                      "2022" = "Q87.1_2022",
                      "2023" = "Q90.1_2023",
                      "2024" = "Q80.1_2024",
                      "2025" = "Q77.1_2025",
                      paste0("Q87.1_", suffix))
  mar_col   <- paste0("Q2_",      suffix)
  liv_col   <- paste0("Q1.1_",    suffix)
  sp_col    <- if (suffix == "2022") "Q28.13_2022"
               else if (suffix == "2023") "Q31.13_2023"
               else if (suffix == "2024") "Q28.13_2024" else NA_character_
  pc_col    <- if (suffix == "2022") "Q28.14_2022"
               else if (suffix == "2023") "Q31.14_2023"
               else if (suffix == "2024") "Q28.14_2024" else NA_character_
  out <- list()
  out$is_female <- as.integer(as_num(ds[[sex_col]]) == 2L)
  age_v <- as_num(ds[[age_col]])
  age_bands <- list(c(18,24), c(25,29), c(30,34), c(35,39), c(40,44),
                    c(45,49), c(50,54), c(55,59), c(60,64))
  for (b in age_bands) {
    out[[sprintf("age_%d_%d", b[1], b[2])]] <-
      as.integer(age_v >= b[1] & age_v <= b[2])
  }
  # Education: drop the Unknown dummy; combine Unknown + NA
  # into the "less than university" reference (very small Unknown share — < 1%
  # in real data). Two non-reference dummies: University (Q21.1 ∈ {6,7,8}) and
  # Graduate (== 9).
  edu <- as_num(ds[[edu_col]])
  out$edu_univ <- as.integer(edu %in% 6:8)
  out$edu_grad <- as.integer(edu == 9)
  emp <- as_num(ds[[emp_col]])
  out$emp_exec    <- as.integer(emp == 1)
  out$emp_self    <- as.integer(emp %in% 2:4)
  out$emp_nonreg  <- as.integer(emp %in% 7:11)
  out$emp_student <- as.integer(emp %in% 12:13)
  out$emp_notwork <- as.integer(emp %in% 14:16 | is.na(emp))
  # Income: replace quantile cuts (uneven on discrete ordinal)
  # with fixed yen-band cuts. <2M = codes 1-4 (0-199万円); 2-6M = codes 5-8
  # (200-599万円); 6-10M = codes 9-12 (600-999万円); 10M+ = codes 13-18
  # (1,000+ 万円); Unknown = codes 19, 20, NA (答えたくない/分からない).
  # Reference = <2M.
  inc <- as_num(ds[[inc_col]])
  inc_band <- ifelse(inc %in% 1:4,   1L,
              ifelse(inc %in% 5:8,   2L,
              ifelse(inc %in% 9:12,  3L,
              ifelse(inc %in% 13:18, 4L, NA_integer_))))
  out$inc_2_6m     <- as.integer(!is.na(inc_band) & inc_band == 2L)
  out$inc_6_10m    <- as.integer(!is.na(inc_band) & inc_band == 3L)
  out$inc_10m_plus <- as.integer(!is.na(inc_band) & inc_band == 4L)
  out$inc_unknown  <- as.integer(is.na(inc_band))
  out$married      <- as.integer(as_num(ds[[mar_col]]) %in% 1:3)
  out$living_alone <- as.integer(!is.na(as_num(ds[[liv_col]])) &
                                 as_num(ds[[liv_col]]) == 1L)
  # Smartphone / PC time: convert band → hours.
  # わからない (band 12) → NA → fill_zero_pre will set to 0 hours later.
  sp_raw <- if (!is.na(sp_col)) as_num(ds[[sp_col]]) else rep(NA_real_, nrow(ds))
  sp_h   <- band_to_hours(sp_raw)
  sp_h[!is.finite(sp_h)] <- 0
  out$smartphone <- sp_h
  pc_raw <- if (!is.na(pc_col)) as_num(ds[[pc_col]]) else rep(NA_real_, nrow(ds))
  pc_h   <- band_to_hours(pc_raw)
  pc_h[!is.finite(pc_h)] <- 0
  out$pc_tablet <- pc_h
  as.data.frame(out, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Per-outcome model fitter (DiD or ANCOVA) for single-exposure (use, intensity).
# Returns one row of summary output.
# ---------------------------------------------------------------------------
fit_one <- function(spec, y_pre_df, y_post_df, conf_df, pre_levels_df,
                    family_label, exposure_vec, exposure_label,
                    mode = c("did", "ancova"),
                    return_vif = FALSE) {
  mode <- match.arg(mode)
  out_id <- spec$id
  y_pre  <- y_pre_df[[out_id]]
  y_post <- y_post_df[[out_id]]
  keep   <- is.finite(y_pre) & is.finite(y_post)
  n_used <- sum(keep)
  is_bin <- spec$kind %in% c("binary", "binary_lt", "binary_ge")
  blank  <- function() {
    base <- data.frame(
      outcome = out_id, exposure = exposure_label,
      estimate = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, n = n_used,
      kind = if (is_bin) "binary" else "continuous",
      family = family_label,
      stringsAsFactors = FALSE
    )
    if (return_vif) list(row = base, vif = NULL) else base
  }
  if (n_used < 30L) return(blank())
  other_pre <- pre_levels_df[, setdiff(names(pre_levels_df), out_id),
                             drop = FALSE]
  if (mode == "did") {
    d_fit <- data.frame(
      y    = as.numeric(y_post - y_pre),
      expo = exposure_vec,
      conf_df,
      other_pre,
      check.names = FALSE
    )[keep, , drop = FALSE]
  } else {
    d_fit <- data.frame(
      y     = as.numeric(y_post),
      y_pre = as.numeric(y_pre),
      expo  = exposure_vec,
      conf_df,
      other_pre,
      check.names = FALSE
    )[keep, , drop = FALSE]
  }
  preds <- drop_constant_cols(d_fit[, setdiff(names(d_fit), "y"), drop = FALSE])
  if (!"expo" %in% names(preds)) return(blank())
  if (length(unique(d_fit$y)) < 2L) return(blank())
  d_fit2 <- cbind(y = d_fit$y, preds)
  fit <- tryCatch(lm(y ~ ., data = d_fit2), error = function(e) NULL)
  if (is.null(fit)) return(blank())
  cf <- robust_coef(fit, "expo")
  row <- data.frame(
    outcome = out_id, exposure = exposure_label,
    estimate = unname(cf["estimate"]), se = unname(cf["se"]),
    z = unname(cf["z"]), p = unname(cf["p"]),
    ci_lo = unname(cf["ci_lo"]), ci_hi = unname(cf["ci_hi"]),
    n = n_used,
    kind = if (is_bin) "binary" else "continuous",
    family = family_label,
    stringsAsFactors = FALSE
  )
  if (return_vif) {
    vifs <- compute_vif_safe(fit)
    return(list(row = row, vif = vifs))
  }
  row
}

# ---------------------------------------------------------------------------
# Per-outcome Pattern B fitter (3 purposes, joint).
# ---------------------------------------------------------------------------
fit_b_one <- function(spec, y_pre_df, y_post_df, conf_df, pre_levels_df,
                      family_label, F1, F2, F3,
                      mode = c("did", "ancova"),
                      return_vif = FALSE) {
  mode <- match.arg(mode)
  out_id <- spec$id
  y_pre  <- y_pre_df[[out_id]]
  y_post <- y_post_df[[out_id]]
  keep   <- is.finite(y_pre) & is.finite(y_post)
  n_used <- sum(keep)
  is_bin <- spec$kind %in% c("binary", "binary_lt", "binary_ge")
  blank_rows <- function() {
    rows <- do.call(rbind, lapply(c("creative","daily","social"), function(t)
      data.frame(outcome = out_id, exposure = t,
        estimate = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_, n = n_used,
        kind = if (is_bin) "binary" else "continuous",
        family = family_label,
        stringsAsFactors = FALSE)))
    if (return_vif) list(rows = rows, vif = NULL) else rows
  }
  if (n_used < 30L) return(blank_rows())
  other_pre <- pre_levels_df[, setdiff(names(pre_levels_df), out_id),
                             drop = FALSE]
  if (mode == "did") {
    d_fit <- data.frame(
      y        = as.numeric(y_post - y_pre),
      creative = F1,
      daily    = F2,
      social   = F3,
      conf_df,
      other_pre,
      check.names = FALSE
    )[keep, , drop = FALSE]
  } else {
    d_fit <- data.frame(
      y        = as.numeric(y_post),
      y_pre    = as.numeric(y_pre),
      creative = F1,
      daily    = F2,
      social   = F3,
      conf_df,
      other_pre,
      check.names = FALSE
    )[keep, , drop = FALSE]
  }
  preds  <- drop_constant_cols(d_fit[, setdiff(names(d_fit), "y"), drop = FALSE])
  needed <- c("creative", "daily", "social")
  if (!all(needed %in% names(preds))) return(blank_rows())
  if (length(unique(d_fit$y)) < 2L)   return(blank_rows())
  d_fit2 <- cbind(y = d_fit$y, preds)
  fit <- tryCatch(lm(y ~ ., data = d_fit2), error = function(e) NULL)
  if (is.null(fit)) return(blank_rows())
  rows <- do.call(rbind, lapply(needed, function(t) {
    cf <- robust_coef(fit, t)
    data.frame(outcome = out_id, exposure = t,
      estimate = unname(cf["estimate"]), se = unname(cf["se"]),
      z = unname(cf["z"]), p = unname(cf["p"]),
      ci_lo = unname(cf["ci_lo"]), ci_hi = unname(cf["ci_hi"]),
      n = n_used,
      kind = if (is_bin) "binary" else "continuous",
      family = family_label,
      stringsAsFactors = FALSE)
  }))
  if (return_vif) {
    vifs <- compute_vif_safe(fit)
    return(list(rows = rows, vif = vifs))
  }
  rows
}

# ---------------------------------------------------------------------------
# Convenience: run the full 5-exposure × 13-outcome family for one tier and
# return one combined data frame (65 rows) with `family` column distinguishing
# A0 / A / B. BH-FDR applied per-exposure (13 tests each) and outcome-z added.
# ---------------------------------------------------------------------------
run_tier_family <- function(ds, y_pre_df, y_post_df, conf_df, pre_levels_df,
                            sd_pre, mode = c("did", "ancova"),
                            family_phase = "did",
                            exposure_use = ds$treat,
                            return_vif = FALSE) {
  mode <- match.arg(mode)
  # Pattern A0 (use)
  a0_rows <- lapply(sens_specs, function(sp)
    fit_one(sp, y_pre_df, y_post_df, conf_df, pre_levels_df,
            family_label = paste0("A0_", family_phase),
            exposure_vec = exposure_use, exposure_label = "use",
            mode = mode))
  results_a0 <- do.call(rbind, a0_rows)

  # Pattern A (intensity) — gather VIF on this fit (used for the diagnostics
  # CSV when return_vif = TRUE)
  a_results <- lapply(sens_specs, function(sp)
    fit_one(sp, y_pre_df, y_post_df, conf_df, pre_levels_df,
            family_label = paste0("A_", family_phase),
            exposure_vec = ds$intensity, exposure_label = "intensity",
            mode = mode, return_vif = return_vif))
  if (return_vif) {
    results_a <- do.call(rbind, lapply(a_results, `[[`, "row"))
    vif_records <- lapply(seq_along(a_results), function(i) {
      v <- a_results[[i]]$vif
      if (is.null(v) || all(is.na(v))) return(NULL)
      data.frame(outcome = sens_specs[[i]]$id, pattern = "A",
                 term = names(v), vif = as.numeric(v),
                 stringsAsFactors = FALSE)
    })
    vif_a <- do.call(rbind, vif_records)
  } else {
    results_a <- do.call(rbind, a_results)
    vif_a <- NULL
  }

  # Pattern B (3 purposes joint)
  b_results <- lapply(sens_specs, function(sp)
    fit_b_one(sp, y_pre_df, y_post_df, conf_df, pre_levels_df,
              family_label = paste0("B_", family_phase),
              F1 = ds$creative_2025, F2 = ds$daily_2025, F3 = ds$social_2025,
              mode = mode, return_vif = return_vif))
  if (return_vif) {
    results_b <- do.call(rbind, lapply(b_results, `[[`, "rows"))
    vif_b_recs <- lapply(seq_along(b_results), function(i) {
      v <- b_results[[i]]$vif
      if (is.null(v) || all(is.na(v))) return(NULL)
      data.frame(outcome = sens_specs[[i]]$id, pattern = "B",
                 term = names(v), vif = as.numeric(v),
                 stringsAsFactors = FALSE)
    })
    vif_b <- do.call(rbind, vif_b_recs)
  } else {
    results_b <- do.call(rbind, b_results)
    vif_b <- NULL
  }

  combined <- rbind(results_a0, results_a, results_b)
  # BH-FDR per-exposure (13 tests each)
  combined$p_bh   <- NA_real_
  combined$bh_sig <- NA
  for (ex in c("use", "intensity", "creative", "daily", "social")) {
    ix <- combined$exposure == ex
    if (any(ix)) {
      combined$p_bh[ix]   <- p.adjust(combined$p[ix], method = "BH")
      combined$bh_sig[ix] <- combined$p_bh[ix] < 0.10
    }
  }
  combined <- add_outcome_z(combined, sd_pre)
  vif_combined <- if (return_vif) {
    rbind(vif_a, vif_b)
  } else NULL
  list(results = combined, vif = vif_combined)
}

# ---------------------------------------------------------------------------
# 4-panel (intensity / creative / daily / social) forest plot for the
# users-only A+B sensitivity (; intensity added). Mirrors the
# existing forest_panels calls but narrows panel_levels to the four
# user-side exposures — Pattern A0 (use, binary D) is degenerate among users
# and excluded.
# ---------------------------------------------------------------------------
forest_users_only_ab <- function(res, n_users, title_tag, x_lab,
                                 tf_dir, file_prefix,
                                 oz_xlab) {
  if (is.null(res) || nrow(res) == 0L) return(invisible(NULL))
  panels4 <- c("intensity", "creative", "daily", "social")
  p_raw <- forest_panels(prep_5col(res), panel_levels = panels4) +
    labs(x = x_lab)
  p_oz  <- forest_panels(swap_to_oz(prep_5col(res)), panel_levels = panels4) +
    labs(x = oz_xlab)
  ggsave(file.path(tf_dir, paste0(file_prefix, ".png")),
         p_raw, width = 13, height = 11, dpi = 150)
  ggsave(file.path(tf_dir, paste0(file_prefix, "_oz.png")),
         p_oz,  width = 13, height = 11, dpi = 150)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Users-only A + B sensitivity (; intensity added per
# ): same Pattern A (intensity) + Pattern B (3 purposes joint)
# regressions as in run_tier_family, but restricted to AI users
# (ds$treat == 1L). β captures within-user dose-response, free of the
# user-vs-non-user contrast that the full-sample patterns bundle in.
# `sd_pre` should be computed on the user subset at the tier's pre-wave so
# the outcome-z standardisation reflects within-user dispersion only.
# ---------------------------------------------------------------------------
run_users_only_ab <- function(ds, y_pre_df, y_post_df, conf_df, pre_levels_df,
                              sd_pre, mode = c("did", "ancova"),
                              family_phase) {
  mode <- match.arg(mode)
  is_user <- which(ds$treat == 1L)
  if (length(is_user) < 30L) {
    cat(sprintf("[users-only A+B] skipped (%s): n_users=%d (<30)\n",
                family_phase, length(is_user)))
    return(NULL)
  }
  ds_u <- ds[is_user, , drop = FALSE]
  yp_u <- y_pre_df[is_user, , drop = FALSE]
  yo_u <- y_post_df[is_user, , drop = FALSE]
  cf_u <- conf_df[is_user, , drop = FALSE]
  pl_u <- pre_levels_df[is_user, , drop = FALSE]
  # Pattern A (intensity)
  rows_a <- lapply(sens_specs, function(sp)
    fit_one(sp, y_pre_df = yp_u, y_post_df = yo_u,
            conf_df = cf_u, pre_levels_df = pl_u,
            family_label  = paste0("A_", family_phase, "_usersonly"),
            exposure_vec  = ds_u$intensity,
            exposure_label = "intensity",
            mode = mode))
  res_a <- do.call(rbind, rows_a)
  # Pattern B (3 purposes joint)
  rows_b <- lapply(sens_specs, function(sp)
    fit_b_one(sp, y_pre_df = yp_u, y_post_df = yo_u,
              conf_df = cf_u, pre_levels_df = pl_u,
              family_label = paste0("B_", family_phase, "_usersonly"),
              F1 = ds_u$creative_2025,
              F2 = ds_u$daily_2025,
              F3 = ds_u$social_2025,
              mode = mode))
  res_b <- do.call(rbind, rows_b)
  out <- rbind(res_a, res_b)
  out$p_bh   <- NA_real_
  out$bh_sig <- NA
  for (ex in c("intensity", "creative", "daily", "social")) {
    ix <- out$exposure == ex
    if (any(ix)) {
      out$p_bh[ix]   <- p.adjust(out$p[ix], method = "BH")
      out$bh_sig[ix] <- out$p_bh[ix] < 0.10
    }
  }
  out <- add_outcome_z(out, sd_pre)
  out
}

# ---------------------------------------------------------------------------
# Forest plotting (5-column unified; raw + outcome-z variants).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))

INDENT <- "    "
y_levels_topdown <- local({
  doms <- unique(domain_of)
  out <- character(0)
  for (dom in doms) {
    out <- c(out, dom)
    ids <- names(domain_of)[domain_of == dom]
    out <- c(out, paste0(INDENT, outcome_labels[ids]))
  }
  out
})
y_levels_factor <- rev(y_levels_topdown)
header_labels <- unique(domain_of)
y_face_vec  <- ifelse(y_levels_factor %in% header_labels, "bold", "plain")
y_color_vec <- ifelse(y_levels_factor %in% header_labels, "grey15", "grey30")
y_size_vec  <- ifelse(y_levels_factor %in% header_labels, 18,    16)

make_y_data <- function(d, panel_levels = NULL) {
  d$y_label <- paste0(INDENT, outcome_labels[d$outcome])
  cols_keep <- names(d)
  build_headers_one <- function(panel_value = NULL) {
    h <- data.frame(
      outcome = paste0("__hdr_", seq_along(header_labels)),
      y_label = header_labels,
      stringsAsFactors = FALSE
    )
    for (col in setdiff(cols_keep, c("outcome", "y_label", "panel"))) {
      h[[col]] <- NA
    }
    if (!is.null(panel_value) && "panel" %in% cols_keep) {
      h$panel <- factor(panel_value, levels = panel_levels)
    }
    h[, cols_keep, drop = FALSE]
  }
  if (is.null(panel_levels)) {
    headers <- build_headers_one()
  } else {
    headers <- do.call(rbind, lapply(panel_levels, build_headers_one))
  }
  out <- rbind(d, headers)
  out$y_label <- factor(out$y_label, levels = y_levels_factor)
  out
}

forest_theme <- theme_classic(base_size = 9) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    strip.background.x = element_rect(fill = "grey95", color = NA),
    strip.text.x       = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(face  = y_face_vec,
                                color = y_color_vec,
                                size  = y_size_vec),
    axis.title.x       = element_text(size = 18),
    axis.text.x = element_text(size = 9),
    legend.text = element_text(size = 15),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.spacing = unit(0.4, "lines")
  )

# Panel strip display labels (; wrapped so the
# longer labels do not clip in the panel strip width). Internal exposure
# names map to the printable column headers used in every forest plot.
panel_display_labels <- c(
  "use"       = "Use",
  "intensity" = "Intensity",
  "creative"  = "Productivity\n/ creative",
  "daily"     = "Daily-life\n/ information",
  "social"    = "Social\n/ emotional"
)

forest_panels <- function(data, panel_levels, title = NULL, subtitle = NULL) {
  # Per : titles and subtitles are no longer rendered; legend
  # label spelled out; panel strips use display names; x-axis label is set
  # by callers via `+ labs(x = ...)`.
  d <- make_y_data(data, panel_levels = panel_levels)
  ggplot(d, aes(x = estimate, y = y_label)) +
    facet_grid(cols = vars(panel),
               labeller = labeller(panel = panel_display_labels)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0, color = "grey40", linewidth = 0.4,
                   na.rm = TRUE) +
    geom_point(aes(fill = sig), shape = 21, size = 2.4, stroke = 0.6,
               color = "black", na.rm = TRUE) +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"),
                      labels = c("FALSE" = "NS",
                                 "TRUE"  = "BH-significant (q < 0.10)"),
                      breaks = c("FALSE", "TRUE"),
                      name = NULL, na.value = NA, na.translate = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = NULL, subtitle = NULL,
         x = "Estimate beta (linear: per 1-point increase in exposure)") +
    forest_theme
}

# Prepare unified 5-column plot data from a combined results data frame.
prep_5col <- function(df_combined) {
  d <- df_combined
  d$sig <- ifelse(is.na(d$bh_sig), FALSE, d$bh_sig)
  d$panel <- factor(d$exposure,
                    levels = c("use", "intensity", "creative", "daily", "social"))
  d
}

# ---------------------------------------------------------------------------
# Per-tier exposure setup helper. Computes ds$treat, ds$intensity,
# ds$creative_2025, ds$daily_2025, ds$social_2025 in-place.
# Cohort filter is applied upstream.
# ---------------------------------------------------------------------------
setup_exposures <- function(ds) {
  ds$intensity <- row_mean(ds, paste0("Q37S3.", 1:9, "_2025"), na_rm = TRUE)
  ds$intensity[!is.finite(ds$intensity)] <- 1
  ds$creative_2025 <- row_mean(ds, paste0("Q37S3.", c(1,2,4,5), "_2025"),
                               na_rm = TRUE)
  ds$daily_2025    <- row_mean(ds, paste0("Q37S3.", c(3,6,7),   "_2025"),
                               na_rm = TRUE)
  ds$social_2025   <- row_mean(ds, paste0("Q37S3.", c(8,9),     "_2025"),
                               na_rm = TRUE)
  for (cn in c("creative_2025", "daily_2025", "social_2025")) {
    v <- ds[[cn]]; v[!is.finite(v)] <- 1; ds[[cn]] <- v
  }
  ds
}

sd_pre_for <- function(y_df) {
  out <- vapply(sens_specs, function(s) {
    v <- y_df[[s$id]]; v <- v[is.finite(v)]
    if (length(v) < 2L) NA_real_ else sd(v)
  }, numeric(1))
  names(out) <- vapply(sens_specs, `[[`, character(1), "id")
  out
}

# ============================================================================
# 7.  M A I N   A N A L Y S I S   (DiD on PANEL_3W_MAIN)
#     Cohort: Q37S1_2025 ∈ {5, 6} vs == 1.   Outcome window: ΔY 2024 → 2025.
#     Pre-trends: ΔY 2023 → 2024 (both pre-exposure for codes 5/6).
# ============================================================================

# INPUT_3WAVE set above (PANEL_3W_MAIN).

main_did_ok       <- FALSE
main_efa_loadings <- NULL
main_efa_fit      <- NULL
main_vif_df       <- NULL
if (!file.exists(INPUT_3WAVE)) {
  cat(sprintf("[main] skipped: %s not found\n", INPUT_3WAVE))
} else {
  df3 <- readr::read_csv(INPUT_3WAVE, show_col_types = FALSE)
  cat(sprintf("[main] 3-wave INPUT=%s n_raw=%d\n", INPUT_3WAVE, nrow(df3)))

  ai_start <- as_num(df3$Q37S1_2025)
  cohort   <- ifelse(ai_start %in% c(5L, 6L), "treat",
              ifelse(ai_start == 1L,         "control", NA_character_))
  keep_co  <- !is.na(cohort)
  ds <- df3[keep_co, , drop = FALSE]
  ds$treat <- as.integer(cohort[keep_co] == "treat")
  cat(sprintf("[main] cohort: treat (Q37S1 in {5,6})=%d  control (Q37S1=1)=%d  excluded=%d\n",
              sum(ds$treat == 1L), sum(ds$treat == 0L),
              sum(!keep_co)))

  ds <- setup_exposures(ds)

  # Materialize Y at 2023, 2024, 2025
  y_2023_df <- materialise_outcomes(ds, "c23")
  y_2024_df <- materialise_outcomes(ds, "c24")
  y_2025_df <- materialise_outcomes(ds, "c25")
  pre_2023_df <- fill_zero_pre(y_2023_df)
  pre_2024_df <- fill_zero_pre(y_2024_df)
  # Time-use outcomes are already in hours with わからない=NA;
  # rows with NA Y_pre or Y_post will be dropped by fit_did_one's `keep` filter.

  conf_2023_df <- derive_confounders(ds, "2023")
  conf_2024_df <- derive_confounders(ds, "2024")

  sd_pre_main_did <- sd_pre_for(y_2024_df)
  sd_pre_main_pre <- sd_pre_for(y_2023_df)

  # PRIMARY DiD (ΔY 2024 → 2025) — also captures VIF
  main_did_fam <- run_tier_family(
    ds              = ds,
    y_pre_df        = y_2024_df,
    y_post_df       = y_2025_df,
    conf_df         = conf_2024_df,
    pre_levels_df   = pre_2024_df,
    sd_pre          = sd_pre_main_did,
    mode            = "did",
    family_phase    = "did",
    exposure_use    = ds$treat,
    return_vif      = TRUE
  )
  res_main_did <- add_e_values_oz(main_did_fam$results)   # E-values main_vif_df  <- main_did_fam$vif
  write.csv(res_main_did, file.path(TF_MAIN, "results_main_did.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Users-only A + B sensitivity (; intensity added per
  # ). SD_pre is computed on the user subset only — for Main this
  # is codes 5+6 = AI initiators during 2025 — so outcome-z reflects
  # within-user dispersion at the immediate pre-wave (2024).
  is_user_main_did       <- which(ds$treat == 1L)
  sd_pre_main_did_users  <- sd_pre_for(y_2024_df[is_user_main_did, , drop = FALSE])
  res_main_did_uo <- run_users_only_ab(
    ds = ds, y_pre_df = y_2024_df, y_post_df = y_2025_df,
    conf_df = conf_2024_df, pre_levels_df = pre_2024_df,
    sd_pre = sd_pre_main_did_users, mode = "did", family_phase = "did"
  )
  if (!is.null(res_main_did_uo)) {
    res_main_did_uo <- add_e_values_oz(res_main_did_uo)
    write.csv(res_main_did_uo,
              file.path(TF_MAIN, "results_main_did_usersonly.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[main] users-only A+B (n_users=%d): sig BH q<0.10 intensity=%d creative=%d daily=%d social=%d\n",
                sum(ds$treat == 1L),
                sum(res_main_did_uo$exposure == "intensity" & res_main_did_uo$bh_sig %in% TRUE),
                sum(res_main_did_uo$exposure == "creative"  & res_main_did_uo$bh_sig %in% TRUE),
                sum(res_main_did_uo$exposure == "daily"     & res_main_did_uo$bh_sig %in% TRUE),
                sum(res_main_did_uo$exposure == "social"    & res_main_did_uo$bh_sig %in% TRUE)))
    forest_users_only_ab(res_main_did_uo, n_users = sum(ds$treat == 1L),
                         title_tag = "Main",
                         x_lab = "Delta change-score beta per +1-pt exposure (DiD 2024->2025)",
                         tf_dir = TF_MAIN,
                         file_prefix = "forest_main_did_usersonly",
                         oz_xlab = "Standardized DiD estimate (β / SD of outcome at baseline)")
    cat("[main] users-only Pattern A+B: saved 1 raw + 1 oz forest plot\n")
  }

  # PRE-TRENDS DiD (ΔY 2023 → 2024)
  res_main_pretrends <- run_tier_family(
    ds              = ds,
    y_pre_df        = y_2023_df,
    y_post_df       = y_2024_df,
    conf_df         = conf_2023_df,
    pre_levels_df   = pre_2023_df,
    sd_pre          = sd_pre_main_pre,
    mode            = "did",
    family_phase    = "pretrends",
    exposure_use    = ds$treat
  )$results
  res_main_pretrends <- add_e_values_oz(res_main_pretrends)
  write.csv(res_main_pretrends, file.path(TF_MAIN, "results_main_pretrends.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("[main] DiD rows=%d (sig BH q<0.10: use=%d intensity=%d creative=%d daily=%d social=%d)\n",
              nrow(res_main_did),
              sum(res_main_did$exposure == "use"       & res_main_did$bh_sig %in% TRUE),
              sum(res_main_did$exposure == "intensity" & res_main_did$bh_sig %in% TRUE),
              sum(res_main_did$exposure == "creative"  & res_main_did$bh_sig %in% TRUE),
              sum(res_main_did$exposure == "daily"     & res_main_did$bh_sig %in% TRUE),
              sum(res_main_did$exposure == "social"    & res_main_did$bh_sig %in% TRUE)))
  cat(sprintf("[main] pre-trends rows=%d (sig BH q<0.10: use=%d intensity=%d creative=%d daily=%d social=%d)\n",
              nrow(res_main_pretrends),
              sum(res_main_pretrends$exposure == "use"       & res_main_pretrends$bh_sig %in% TRUE),
              sum(res_main_pretrends$exposure == "intensity" & res_main_pretrends$bh_sig %in% TRUE),
              sum(res_main_pretrends$exposure == "creative"  & res_main_pretrends$bh_sig %in% TRUE),
              sum(res_main_pretrends$exposure == "daily"     & res_main_pretrends$bh_sig %in% TRUE),
              sum(res_main_pretrends$exposure == "social"    & res_main_pretrends$bh_sig %in% TRUE)))

  # Forest plots
  ns <- sum(ds$treat == 1L); nc <- sum(ds$treat == 0L)
  panels5 <- c("use", "intensity", "creative", "daily", "social")
  did_xlab <- "Delta change-score beta (use: 0/1; others: 1-pt scale)"
  oz_xlab  <- "Standardized DiD estimate (β / SD of outcome at baseline)"

  p_main_did <- forest_panels(prep_5col(res_main_did),
    panel_levels = panels5,
    title = "Main - DiD - DeltaY (2024->2025) by 5 exposures",
    subtitle = sprintf("n_treat = %d  |  n_control = %d  |  white = NS, black = BH q<0.10", ns, nc)) +
    labs(x = did_xlab)
  p_main_pre <- forest_panels(prep_5col(res_main_pretrends),
    panel_levels = panels5,
    title = "Main - pre-trends DeltaY (2023->2024) by future 5 exposures",
    subtitle = "Significant non-zero beta here flags parallel-trends violation") +
    labs(x = did_xlab)
  p_main_did_oz <- forest_panels(swap_to_oz(prep_5col(res_main_did)),
    panel_levels = panels5,
    title = "Main - DiD DeltaY (2024->2025) - outcome-standardised",
    subtitle = sprintf("n_treat = %d  |  n_control = %d", ns, nc)) +
    labs(x = oz_xlab)
  p_main_pre_oz <- forest_panels(swap_to_oz(prep_5col(res_main_pretrends)),
    panel_levels = panels5,
    title = "Main - pre-trends DeltaY (2023->2024) - outcome-standardised",
    subtitle = "Significant non-zero beta here flags parallel-trends violation") +
    labs(x = oz_xlab)

  ggsave(file.path(TF_MAIN, "forest_main_did.png"),       p_main_did,     width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_MAIN, "forest_main_pretrends.png"), p_main_pre,     width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_MAIN, "forest_main_did_oz.png"),    p_main_did_oz,  width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_MAIN, "forest_main_pretrends_oz.png"), p_main_pre_oz, width = 13, height = 11, dpi = 300)
  cat("[main] saved 2 raw + 2 oz forest plots\n")

  # Stash for §9 (Sens 2 ANCOVA reuses the same cohort-filtered ds + Y panels).
  main_ds            <- ds
  main_y_2024_df     <- y_2024_df
  main_y_2025_df     <- y_2025_df
  main_pre_2024_df   <- pre_2024_df
  main_conf_2024_df  <- conf_2024_df
  main_sd_pre_did    <- sd_pre_main_did

  # ---- EFA + Cronbach's α on Main treated arm only -------
  # Treated arm = codes 5, 6 = 2025 initiators (n = 1,940 on real data).
  ds_treat <- ds[ds$treat == 1L, , drop = FALSE]
  purpose_cols_main <- paste0("Q37S3.", 1:9, "_2025")
  efa_main <- tryCatch(
    suppressMessages(suppressWarnings(
      psych::fa(ds_treat[, purpose_cols_main], nfactors = 3,
                rotate = "oblimin", fm = "ml", scores = "regression")
    )),
    error = function(e) { cat("[main] EFA failed:", conditionMessage(e), "\n"); NULL }
  )
  purpose_labels_short_main <- c(
    "1.writing/SNS", "2.translation/summary", "3.lookup/Q&A",
    "4.image/video gen", "5.learning/study", "6.daily-life planning",
    "7.health info", "8.conversation", "9.emotional support"
  )
  if (!is.null(efa_main)) {
    Lm <- unclass(efa_main$loadings)
    rownames(Lm) <- purpose_labels_short_main[seq_len(nrow(Lm))]
    main_efa_loadings <- data.frame(
      item = rownames(Lm),
      ML1  = round(as.numeric(Lm[, 1]), 3),
      ML2  = round(as.numeric(Lm[, 2]), 3),
      ML3  = round(as.numeric(Lm[, 3]), 3),
      stringsAsFactors = FALSE
    )
    rmsea_m <- efa_main$RMSEA
    cumvar_m <- as.numeric(tail(efa_main$Vaccounted["Cumulative Var", ], 1))
    .alpha_safe <- function(items) {
      tryCatch(
        suppressMessages(suppressWarnings(
          psych::alpha(items, check.keys = FALSE, warnings = FALSE)$total$raw_alpha
        )),
        error = function(e) NA_real_
      )
    }
    a_int_m       <- .alpha_safe(ds_treat[, purpose_cols_main])
    a_creative_m  <- .alpha_safe(ds_treat[, purpose_cols_main[c(1, 2, 4, 5)]])
    a_daily_m     <- .alpha_safe(ds_treat[, purpose_cols_main[c(3, 6, 7)]])
    a_social_m    <- .alpha_safe(ds_treat[, purpose_cols_main[c(8, 9)]])
    main_efa_fit <- data.frame(
      metric = c("n_users", "TLI", "RMSEA", "RMSEA_lower", "RMSEA_upper",
                 "BIC", "cumulative_var",
                 "alpha_intensity_9items", "alpha_creative_4items",
                 "alpha_daily_3items", "alpha_social_2items"),
      value  = c(nrow(ds_treat),
                 as.numeric(efa_main$TLI),
                 as.numeric(rmsea_m[["RMSEA"]]),
                 as.numeric(rmsea_m[["lower"]]),
                 as.numeric(rmsea_m[["upper"]]),
                 as.numeric(efa_main$BIC),
                 cumvar_m,
                 a_int_m, a_creative_m, a_daily_m, a_social_m),
      stringsAsFactors = FALSE
    )
    cat(sprintf("[main] EFA (treated arm n=%d): TLI=%.3f RMSEA=%.3f [%.3f, %.3f] BIC=%.1f cumvar=%.3f\n",
                nrow(ds_treat), as.numeric(efa_main$TLI),
                as.numeric(rmsea_m[["RMSEA"]]),
                as.numeric(rmsea_m[["lower"]]),
                as.numeric(rmsea_m[["upper"]]),
                as.numeric(efa_main$BIC), cumvar_m))
    cat(sprintf("[main] Cronbach's alpha: intensity=%.2f creative=%.2f daily=%.2f social=%.2f\n",
                a_int_m, a_creative_m, a_daily_m, a_social_m))
  }

  cat("[main] [done]\n")
  main_did_ok <- TRUE
}

# ============================================================================
# 8.  S E N S   1   (DiD on PANEL_3W_SENS1)
#     Cohort: Q37S1_2025 ∈ {4, 5, 6} vs == 1. Outcome window: ΔY 2023 → 2025.
#     Pre-trends: ΔY 2022 → 2023.
# ============================================================================

# INPUT_S1 set above (PANEL_3W_SENS1).

if (!file.exists(INPUT_S1)) {
  cat(sprintf("[sens1] skipped: %s not found\n", INPUT_S1))
} else {
  df_s1 <- readr::read_csv(INPUT_S1, show_col_types = FALSE)
  cat(sprintf("[sens1] 3-wave (2022/23/25) INPUT=%s n_raw=%d\n",
              INPUT_S1, nrow(df_s1)))

  ai_start <- as_num(df_s1$Q37S1_2025)
  cohort   <- ifelse(ai_start %in% c(4L, 5L, 6L), "treat",
              ifelse(ai_start == 1L,             "control", NA_character_))
  keep_co  <- !is.na(cohort)
  ds <- df_s1[keep_co, , drop = FALSE]
  ds$treat <- as.integer(cohort[keep_co] == "treat")
  cat(sprintf("[sens1] cohort: treat (Q37S1 in {4,5,6})=%d  control (Q37S1=1)=%d  excluded=%d\n",
              sum(ds$treat == 1L), sum(ds$treat == 0L),
              sum(!keep_co)))

  ds <- setup_exposures(ds)

  y_2022_df <- materialise_outcomes(ds, "c22")
  y_2023_df <- materialise_outcomes(ds, "c23")
  y_2025_df <- materialise_outcomes(ds, "c25")
  pre_2022_df <- fill_zero_pre(y_2022_df)
  pre_2023_df <- fill_zero_pre(y_2023_df)

  conf_2022_df <- derive_confounders(ds, "2022")
  conf_2023_df <- derive_confounders(ds, "2023")

  sd_pre_s1_did <- sd_pre_for(y_2023_df)
  sd_pre_s1_pre <- sd_pre_for(y_2022_df)

  # PRIMARY DiD (ΔY 2023 → 2025)
  res_s1_did <- run_tier_family(
    ds            = ds,
    y_pre_df      = y_2023_df,
    y_post_df     = y_2025_df,
    conf_df       = conf_2023_df,
    pre_levels_df = pre_2023_df,
    sd_pre        = sd_pre_s1_did,
    mode          = "did",
    family_phase  = "did",
    exposure_use  = ds$treat
  )$results
  res_s1_did <- add_e_values_oz(res_s1_did)
  write.csv(res_s1_did, file.path(TF_SENS1, "results_sens1_did.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Users-only A + B sensitivity (; intensity added per
  # ). SD_pre on user subset only at 2023 baseline.
  is_user_s1_did       <- which(ds$treat == 1L)
  sd_pre_s1_did_users  <- sd_pre_for(y_2023_df[is_user_s1_did, , drop = FALSE])
  res_s1_did_uo <- run_users_only_ab(
    ds = ds, y_pre_df = y_2023_df, y_post_df = y_2025_df,
    conf_df = conf_2023_df, pre_levels_df = pre_2023_df,
    sd_pre = sd_pre_s1_did_users, mode = "did", family_phase = "did"
  )
  if (!is.null(res_s1_did_uo)) {
    res_s1_did_uo <- add_e_values_oz(res_s1_did_uo)
    write.csv(res_s1_did_uo,
              file.path(TF_SENS1, "results_sens1_did_usersonly.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[sens1] users-only A+B (n_users=%d): sig BH q<0.10 intensity=%d creative=%d daily=%d social=%d\n",
                sum(ds$treat == 1L),
                sum(res_s1_did_uo$exposure == "intensity" & res_s1_did_uo$bh_sig %in% TRUE),
                sum(res_s1_did_uo$exposure == "creative"  & res_s1_did_uo$bh_sig %in% TRUE),
                sum(res_s1_did_uo$exposure == "daily"     & res_s1_did_uo$bh_sig %in% TRUE),
                sum(res_s1_did_uo$exposure == "social"    & res_s1_did_uo$bh_sig %in% TRUE)))
    forest_users_only_ab(res_s1_did_uo, n_users = sum(ds$treat == 1L),
                         title_tag = "Sens 1",
                         x_lab = "Delta change-score beta per +1-pt exposure (DiD 2023->2025)",
                         tf_dir = TF_SENS1,
                         file_prefix = "forest_sens1_did_usersonly",
                         oz_xlab = "Standardized DiD estimate (β / SD of outcome at baseline)")
    cat("[sens1] users-only Pattern A+B: saved 1 raw + 1 oz forest plot\n")
  }

  # PRE-TRENDS (ΔY 2022 → 2023)
  res_s1_pretrends <- run_tier_family(
    ds            = ds,
    y_pre_df      = y_2022_df,
    y_post_df     = y_2023_df,
    conf_df       = conf_2022_df,
    pre_levels_df = pre_2022_df,
    sd_pre        = sd_pre_s1_pre,
    mode          = "did",
    family_phase  = "pretrends",
    exposure_use  = ds$treat
  )$results
  res_s1_pretrends <- add_e_values_oz(res_s1_pretrends)
  write.csv(res_s1_pretrends, file.path(TF_SENS1, "results_sens1_pretrends.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  cat(sprintf("[sens1] DiD rows=%d (sig BH q<0.10: use=%d intensity=%d creative=%d daily=%d social=%d)\n",
              nrow(res_s1_did),
              sum(res_s1_did$exposure == "use"       & res_s1_did$bh_sig %in% TRUE),
              sum(res_s1_did$exposure == "intensity" & res_s1_did$bh_sig %in% TRUE),
              sum(res_s1_did$exposure == "creative"  & res_s1_did$bh_sig %in% TRUE),
              sum(res_s1_did$exposure == "daily"     & res_s1_did$bh_sig %in% TRUE),
              sum(res_s1_did$exposure == "social"    & res_s1_did$bh_sig %in% TRUE)))
  cat(sprintf("[sens1] pre-trends rows=%d (sig BH q<0.10: use=%d intensity=%d creative=%d daily=%d social=%d)\n",
              nrow(res_s1_pretrends),
              sum(res_s1_pretrends$exposure == "use"       & res_s1_pretrends$bh_sig %in% TRUE),
              sum(res_s1_pretrends$exposure == "intensity" & res_s1_pretrends$bh_sig %in% TRUE),
              sum(res_s1_pretrends$exposure == "creative"  & res_s1_pretrends$bh_sig %in% TRUE),
              sum(res_s1_pretrends$exposure == "daily"     & res_s1_pretrends$bh_sig %in% TRUE),
              sum(res_s1_pretrends$exposure == "social"    & res_s1_pretrends$bh_sig %in% TRUE)))

  ns <- sum(ds$treat == 1L); nc <- sum(ds$treat == 0L)
  panels5 <- c("use", "intensity", "creative", "daily", "social")
  did_xlab <- "Delta change-score beta (use: 0/1; others: 1-pt scale)"
  oz_xlab  <- "Standardized DiD estimate (β / SD of outcome at baseline)"

  p_s1_did <- forest_panels(prep_5col(res_s1_did),
    panel_levels = panels5,
    title = "Sens 1 - DiD DeltaY (2023->2025) by 5 exposures",
    subtitle = sprintf("n_treat = %d  |  n_control = %d  |  white = NS, black = BH q<0.10", ns, nc)) +
    labs(x = did_xlab)
  p_s1_pre <- forest_panels(prep_5col(res_s1_pretrends),
    panel_levels = panels5,
    title = "Sens 1 - pre-trends DeltaY (2022->2023) by future 5 exposures",
    subtitle = "Significant non-zero beta here flags parallel-trends violation") +
    labs(x = did_xlab)
  p_s1_did_oz <- forest_panels(swap_to_oz(prep_5col(res_s1_did)),
    panel_levels = panels5,
    title = "Sens 1 - DiD DeltaY (2023->2025) - outcome-standardised",
    subtitle = sprintf("n_treat = %d  |  n_control = %d", ns, nc)) + labs(x = oz_xlab)
  p_s1_pre_oz <- forest_panels(swap_to_oz(prep_5col(res_s1_pretrends)),
    panel_levels = panels5,
    title = "Sens 1 - pre-trends DeltaY (2022->2023) - outcome-standardised",
    subtitle = "Significant non-zero beta here flags parallel-trends violation") + labs(x = oz_xlab)

  ggsave(file.path(TF_SENS1, "forest_sens1_did.png"),       p_s1_did,     width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_SENS1, "forest_sens1_pretrends.png"), p_s1_pre,     width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_SENS1, "forest_sens1_did_oz.png"),    p_s1_did_oz,  width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_SENS1, "forest_sens1_pretrends_oz.png"), p_s1_pre_oz, width = 13, height = 11, dpi = 300)
  cat("[sens1] saved 2 raw + 2 oz forest plots\n")

  cat("[sens1] [done]\n")
}

# ============================================================================
# 9.  S E N S   2   (ANCOVA on PANEL_3W_MAIN; same cohort as Main)
#     Y_2025 ~ exposure + Y_2024 + 9 confounders @ 2024 + 12 other @ 2024.
# ============================================================================

if (!main_did_ok) {
  cat("[sens2] skipped: Main did not run (no 3-wave panel)\n")
} else {
  ds              <- main_ds
  y_2024_df       <- main_y_2024_df
  y_2025_df       <- main_y_2025_df
  pre_2024_df     <- main_pre_2024_df
  conf_2024_df    <- main_conf_2024_df
  sd_pre_sens2    <- main_sd_pre_did   # same baseline (2024) SD as Main DiD

  res_sens2 <- run_tier_family(
    ds            = ds,
    y_pre_df      = y_2024_df,
    y_post_df     = y_2025_df,
    conf_df       = conf_2024_df,
    pre_levels_df = pre_2024_df,
    sd_pre        = sd_pre_sens2,
    mode          = "ancova",
    family_phase  = "ancova",
    exposure_use  = ds$treat
  )$results
  res_sens2 <- add_e_values_oz(res_sens2)
  write.csv(res_sens2, file.path(TF_SENS2, "results_sens2_ancova.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Users-only A + B sensitivity (; intensity added per
  # ). SD_pre on user subset only at 2024 baseline (same cohort
  # as Main).
  is_user_sens2          <- which(ds$treat == 1L)
  sd_pre_sens2_users     <- sd_pre_for(y_2024_df[is_user_sens2, , drop = FALSE])
  res_sens2_uo <- run_users_only_ab(
    ds = ds, y_pre_df = y_2024_df, y_post_df = y_2025_df,
    conf_df = conf_2024_df, pre_levels_df = pre_2024_df,
    sd_pre = sd_pre_sens2_users, mode = "ancova", family_phase = "ancova"
  )
  if (!is.null(res_sens2_uo)) {
    res_sens2_uo <- add_e_values_oz(res_sens2_uo)
    write.csv(res_sens2_uo,
              file.path(TF_SENS2, "results_sens2_ancova_usersonly.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[sens2] users-only A+B (n_users=%d): sig BH q<0.10 intensity=%d creative=%d daily=%d social=%d\n",
                sum(ds$treat == 1L),
                sum(res_sens2_uo$exposure == "intensity" & res_sens2_uo$bh_sig %in% TRUE),
                sum(res_sens2_uo$exposure == "creative"  & res_sens2_uo$bh_sig %in% TRUE),
                sum(res_sens2_uo$exposure == "daily"     & res_sens2_uo$bh_sig %in% TRUE),
                sum(res_sens2_uo$exposure == "social"    & res_sens2_uo$bh_sig %in% TRUE)))
    forest_users_only_ab(res_sens2_uo, n_users = sum(ds$treat == 1L),
                         title_tag = "Sens 2",
                         x_lab = "ANCOVA beta per +1-pt exposure (Y_2025 ~ exposure + Y_2024 + cov)",
                         tf_dir = TF_SENS2,
                         file_prefix = "forest_sens2_ancova_usersonly",
                         oz_xlab = "Standardized ANCOVA estimate (β / SD of outcome at baseline)")
    cat("[sens2] users-only Pattern A+B: saved 1 raw + 1 oz forest plot\n")
  }

  cat(sprintf("[sens2] ANCOVA rows=%d (sig BH q<0.10: use=%d intensity=%d creative=%d daily=%d social=%d)\n",
              nrow(res_sens2),
              sum(res_sens2$exposure == "use"       & res_sens2$bh_sig %in% TRUE),
              sum(res_sens2$exposure == "intensity" & res_sens2$bh_sig %in% TRUE),
              sum(res_sens2$exposure == "creative"  & res_sens2$bh_sig %in% TRUE),
              sum(res_sens2$exposure == "daily"     & res_sens2$bh_sig %in% TRUE),
              sum(res_sens2$exposure == "social"    & res_sens2$bh_sig %in% TRUE)))

  ns <- sum(ds$treat == 1L); nc <- sum(ds$treat == 0L)
  panels5 <- c("use", "intensity", "creative", "daily", "social")
  ancova_xlab <- "ANCOVA beta (Y_2025 ~ exposure + Y_2024 + covariates)"
  oz_xa       <- "Standardized ANCOVA estimate (β / SD of outcome at baseline)"

  p_sens2 <- forest_panels(prep_5col(res_sens2),
    panel_levels = panels5,
    title = "Sens 2 - ANCOVA Y_2025 ~ exposure + Y_2024 + covariates",
    subtitle = sprintf("n_treat = %d  |  n_control = %d  |  white = NS, black = BH q<0.10", ns, nc)) +
    labs(x = ancova_xlab)
  p_sens2_oz <- forest_panels(swap_to_oz(prep_5col(res_sens2)),
    panel_levels = panels5,
    title = "Sens 2 - ANCOVA - outcome-standardised",
    subtitle = sprintf("n_treat = %d  |  n_control = %d", ns, nc)) + labs(x = oz_xa)

  ggsave(file.path(TF_SENS2, "forest_sens2_ancova.png"),    p_sens2,    width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_SENS2, "forest_sens2_ancova_oz.png"), p_sens2_oz, width = 13, height = 11, dpi = 300)
  cat("[sens2] saved 1 raw + 1 oz forest plot\n")
  cat("[sens2] [done]\n")
}

# ============================================================================
# 10.  S E N S   3   (ANCOVA on PANEL_2W; full 2-wave sample)
#      Y_2025 ~ exposure + Y_2022 + 9 confounders @ 2022 + 12 other @ 2022.
# ============================================================================

INPUT_S3 <- INPUT_FILE   # PANEL_2W or PANEL_2W

# Stash globals so §11 (Diagnostics) can reuse them.
# (EFA + VIF are now computed in §7 Main on the Main treated arm and Main DiD
# design matrices respectively, ; Sens 3 no longer computes them.)
sens3_done <- FALSE

if (!file.exists(INPUT_S3)) {
  cat(sprintf("[sens3] skipped: %s not found\n", INPUT_S3))
} else {
  df_s3_raw <- readr::read_csv(INPUT_S3, show_col_types = FALSE)
  cat(sprintf("[sens3] 2-wave INPUT=%s n_raw=%d\n", INPUT_S3, nrow(df_s3_raw)))

  # Full sample: D = 1 if Q37S1_2025 in {3,4,5,6}, D = 0 if
  # == 1 (never used). Code 2 (past user, not now — quit) is now EXCLUDED for
  # consistency with Main / Sens 1 / Sens 2.
  df_s3_raw <- df_s3_raw %>%
    mutate(
      Q37S1_2025 = as_num(Q37S1_2025),
      D = case_when(
        Q37S1_2025 == 1                ~ 0L,
        Q37S1_2025 %in% c(3, 4, 5, 6)  ~ 1L,
        TRUE                           ~ NA_integer_
      ),
      AGE_2022 = as_num(AGE_2022)
    )
  n_pre <- nrow(df_s3_raw)
  ds <- df_s3_raw %>% filter(AGE_2022 >= 18 & !is.na(D))
  cat(sprintf("[sens3] sample: age>=18 & !is.na(D): n=%d -> %d\n", n_pre, nrow(ds)))

  # Q37S3 NA → 1 floor for never-users (skip-pattern)
  purpose_cols <- paste0("Q37S3.", 1:9, "_2025")
  for (c0 in purpose_cols) {
    v <- as_num(ds[[c0]]); v[is.na(v)] <- 1; ds[[c0]] <- v
  }
  ds$treat        <- as.integer(ds$D)   # alias used by run_tier_family (use exposure)
  ds <- setup_exposures(ds)

  y_2022_df <- materialise_outcomes(ds, "c22")
  y_2025_df <- materialise_outcomes(ds, "c25")
  pre_2022_df <- fill_zero_pre(y_2022_df)
  conf_2022_df <- derive_confounders(ds, "2022")

  sd_pre_sens3 <- sd_pre_for(y_2022_df)

  res_sens3 <- run_tier_family(
    ds            = ds,
    y_pre_df      = y_2022_df,
    y_post_df     = y_2025_df,
    conf_df       = conf_2022_df,
    pre_levels_df = pre_2022_df,
    sd_pre        = sd_pre_sens3,
    mode          = "ancova",
    family_phase  = "ancova",
    exposure_use  = ds$treat
  )$results
  res_sens3 <- add_e_values_oz(res_sens3)
  write.csv(res_sens3, file.path(TF_SENS3, "results_sens3_ancova.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Users-only A + B sensitivity (; intensity added per
  # ). SD_pre on user subset only at 2022 baseline.
  is_user_sens3          <- which(ds$treat == 1L)
  sd_pre_sens3_users     <- sd_pre_for(y_2022_df[is_user_sens3, , drop = FALSE])
  res_sens3_uo <- run_users_only_ab(
    ds = ds, y_pre_df = y_2022_df, y_post_df = y_2025_df,
    conf_df = conf_2022_df, pre_levels_df = pre_2022_df,
    sd_pre = sd_pre_sens3_users, mode = "ancova", family_phase = "ancova"
  )
  if (!is.null(res_sens3_uo)) {
    res_sens3_uo <- add_e_values_oz(res_sens3_uo)
    write.csv(res_sens3_uo,
              file.path(TF_SENS3, "results_sens3_ancova_usersonly.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[sens3] users-only A+B (n_users=%d): sig BH q<0.10 intensity=%d creative=%d daily=%d social=%d\n",
                sum(ds$treat == 1L),
                sum(res_sens3_uo$exposure == "intensity" & res_sens3_uo$bh_sig %in% TRUE),
                sum(res_sens3_uo$exposure == "creative"  & res_sens3_uo$bh_sig %in% TRUE),
                sum(res_sens3_uo$exposure == "daily"     & res_sens3_uo$bh_sig %in% TRUE),
                sum(res_sens3_uo$exposure == "social"    & res_sens3_uo$bh_sig %in% TRUE)))
    forest_users_only_ab(res_sens3_uo, n_users = sum(ds$treat == 1L),
                         title_tag = "Sens 3",
                         x_lab = "ANCOVA beta per +1-pt exposure (Y_2025 ~ exposure + Y_2022 + cov)",
                         tf_dir = TF_SENS3,
                         file_prefix = "forest_sens3_ancova_usersonly",
                         oz_xlab = "Standardized ANCOVA estimate (β / SD of outcome at baseline)")
    cat("[sens3] users-only Pattern A+B: saved 1 raw + 1 oz forest plot\n")
  }

  cat(sprintf("[sens3] ANCOVA rows=%d (sig BH q<0.10: use=%d intensity=%d creative=%d daily=%d social=%d)\n",
              nrow(res_sens3),
              sum(res_sens3$exposure == "use"       & res_sens3$bh_sig %in% TRUE),
              sum(res_sens3$exposure == "intensity" & res_sens3$bh_sig %in% TRUE),
              sum(res_sens3$exposure == "creative"  & res_sens3$bh_sig %in% TRUE),
              sum(res_sens3$exposure == "daily"     & res_sens3$bh_sig %in% TRUE),
              sum(res_sens3$exposure == "social"    & res_sens3$bh_sig %in% TRUE)))

  # Forest plots — same set as Sens 2 (raw + oz)
  ns <- sum(ds$treat == 1L); nc <- sum(ds$treat == 0L)
  panels5 <- c("use", "intensity", "creative", "daily", "social")
  ancova_xlab3 <- "ANCOVA beta (Y_2025 ~ exposure + Y_2022 + covariates)"
  oz_xa3       <- "Standardized ANCOVA estimate (β / SD of outcome at baseline)"

  p_sens3 <- forest_panels(prep_5col(res_sens3),
    panel_levels = panels5,
    title = "Sens 3 - ANCOVA Y_2025 ~ exposure + Y_2022 + covariates (full sample)",
    subtitle = sprintf("n_D1 = %d  |  n_D0 = %d  |  white = NS, black = BH q<0.10", ns, nc)) +
    labs(x = ancova_xlab3)
  p_sens3_oz <- forest_panels(swap_to_oz(prep_5col(res_sens3)),
    panel_levels = panels5,
    title = "Sens 3 - ANCOVA - outcome-standardised",
    subtitle = sprintf("n_D1 = %d  |  n_D0 = %d", ns, nc)) + labs(x = oz_xa3)

  ggsave(file.path(TF_SENS3, "forest_sens3_ancova.png"),    p_sens3,    width = 13, height = 11, dpi = 300)
  ggsave(file.path(TF_SENS3, "forest_sens3_ancova_oz.png"), p_sens3_oz, width = 13, height = 11, dpi = 300)
  cat("[sens3] saved 1 raw + 1 oz forest plot\n")

  cat("[sens3] [done]\n")
  sens3_done <- TRUE
}

# ============================================================================
# 11.  D I A G N O S T I C S   (EFA loadings + fit + VIF — from Main)
#      EFA + α: on Main treated arm only (codes 5,6 in 3-wave panel).
#      VIF:     on Main DiD design matrices (full Main cohort, n=7,674).
# ============================================================================

if (main_did_ok) {
  if (!is.null(main_efa_loadings)) {
    write.csv(main_efa_loadings,
              file.path(DG_DIR, "efa_loadings.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[diagnostics] saved %s\n",
                file.path(DG_DIR, "efa_loadings.csv")))
  } else {
    cat("[diagnostics] EFA loadings unavailable\n")
  }

  if (!is.null(main_efa_fit)) {
    write.csv(main_efa_fit,
              file.path(DG_DIR, "efa_fit.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[diagnostics] saved %s\n",
                file.path(DG_DIR, "efa_fit.csv")))
  }

  if (!is.null(main_vif_df) && nrow(main_vif_df) > 0) {
    write.csv(main_vif_df, file.path(DG_DIR, "vif.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("[diagnostics] VIF rows=%d max VIF=%.2f >5: %d  >10: %d\n",
                nrow(main_vif_df),
                max(main_vif_df$vif, na.rm = TRUE),
                sum(main_vif_df$vif > 5,  na.rm = TRUE),
                sum(main_vif_df$vif > 10, na.rm = TRUE)))
  } else {
    cat("[diagnostics] VIF unavailable\n")
  }
  cat("[diagnostics] [done]\n")
} else {
  cat("[diagnostics] skipped: Main DiD did not run\n")
}

cat("[done]\n")

# ============================================================================
# 12.  P U B L I C A T I O N   T A B L E S
#      Tables 1, 2, 3 + Sup 1A, 1B, 2, 3.
#      All saved at MS_DIR top level.
# ============================================================================

# Self-contained re-derivation of confounders + outcomes from panel CSVs
# (mirrors the in-script helpers but with .tc_ prefix to avoid collisions).

.tc_input_file <- PANEL_2W
.tc_top_dir    <- OUTPUT_DIR
.tc_out_dir    <- MS_DIR   # publishable tables go to manuscript/ top level

.tc_as_num <- function(x) suppressWarnings(as.numeric(x))

.tc_row_mean <- function(d, cols, na_rm = TRUE) {
  sub <- as.data.frame(d[, cols, drop = FALSE])
  m <- matrix(.tc_as_num(unlist(sub, use.names = FALSE)),
              nrow = nrow(sub), ncol = length(cols))
  rowMeans(m, na.rm = na_rm)
}
.tc_row_sum_tx <- function(d, cols, fn = identity, na_rm = TRUE) {
  sub <- as.data.frame(d[, cols, drop = FALSE])
  m <- matrix(.tc_as_num(unlist(sub, use.names = FALSE)),
              nrow = nrow(sub), ncol = length(cols))
  m2 <- fn(m)
  if (!is.matrix(m2)) m2 <- matrix(as.numeric(m2), nrow = nrow(m), ncol = ncol(m))
  rowSums(m2, na.rm = na_rm)
}
.tc_k6_recode   <- function(M) pmax(0, pmin(4, 5 - M))
.tc_ucla_recode <- function(M) pmax(0, pmin(3, 4 - M))

.tc_outcome_specs <- list(
  list(id="d1_pm_avg",   label="Physical & mental health (0-10)", domain="1. Mental & Physical Health",
       kind="mean", c22=c("Q78.3_2022","Q78.4_2022"), c23=c("Q86.3_2023","Q86.4_2023"),
       c24=c("Q76.3_2024","Q76.4_2024"), c25=c("Q74.3_2025","Q74.4_2025"), ndp=2),
  list(id="d1_k6",       label="K6 distress (sum 0-24)", domain="1. Mental & Physical Health",
       kind="sum_recode", recode=.tc_k6_recode,
       c22=paste0("Q68.", 1:6, "_2022"), c23=paste0("Q75.", 1:6, "_2023"),
       c24=paste0("Q65.", 1:6, "_2024"), c25=paste0("Q65.", 1:6, "_2025"), ndp=2),
  list(id="d2_hls_avg",  label="Happiness & life satisfaction (0-10)", domain="2. Happiness & Life Satisfaction",
       kind="mean", c22=c("Q78.1_2022","Q78.2_2022"), c23=c("Q86.1_2023","Q86.2_2023"),
       c24=c("Q76.1_2024","Q76.2_2024"), c25=c("Q74.1_2025","Q74.2_2025"), ndp=2),
  list(id="d3_mp_avg",   label="Meaning & purpose (0-10)", domain="3. Meaning & Purpose",
       kind="mean", c22=c("Q78.5_2022","Q78.14_2022"), c23=c("Q86.5_2023","Q86.14_2023"),
       c24=c("Q76.5_2024","Q76.14_2024"), c25=c("Q74.5_2025","Q74.14_2025"), ndp=2),
  list(id="d3_ikigai",   label="Ikigai (0-10)", domain="3. Meaning & Purpose",
       kind="single", c22="Q78.6_2022", c23="Q86.6_2023", c24="Q76.6_2024", c25="Q74.6_2025", ndp=2),
  list(id="d4_cv_avg",   label="Character & virtue (0-10)", domain="4. Character & Virtue",
       kind="mean", c22=c("Q78.8_2022","Q78.15_2022"), c23=c("Q86.8_2023","Q86.15_2023"),
       c24=c("Q76.8_2024","Q76.15_2024"), c25=c("Q74.8_2025","Q74.15_2025"), ndp=2),
  list(id="d5_rel_avg",  label="Relations (0-10)", domain="5. Social Wellbeing",
       kind="mean", c22=c("Q78.9_2022","Q78.10_2022"), c23=c("Q86.9_2023","Q86.10_2023"),
       c24=c("Q76.9_2024","Q76.10_2024"), c25=c("Q74.9_2025","Q74.10_2025"), ndp=2),
  list(id="d5_lsns_total", label="LSNS-6 total (sum 0-30)", domain="5. Social Wellbeing",
       kind="sum",
       c22=paste0("Q18.", 1:6, "_2022"), c23=paste0("Q18.", 1:6, "_2023"),
       c24=paste0("Q17.", 1:6, "_2024"), c25=paste0("Q20.", 1:6, "_2025"), ndp=2),
  list(id="d5_ucla",     label="UCLA-3 loneliness (sum 0-9)", domain="5. Social Wellbeing",
       kind="sum_recode", recode=.tc_ucla_recode,
       c22=paste0("Q68S1.", 1:3, "_2022"), c23=paste0("Q76.", 1:3, "_2023"),
       c24=paste0("Q66.", 1:3, "_2024"), c25=paste0("Q66.", 1:3, "_2025"), ndp=2),
  list(id="d6_worry_avg", label="Financial & safety worry (0-10)", domain="6. Financial & Material Stability",
       kind="mean", c22=c("Q78.11_2022","Q78.12_2022"), c23=c("Q86.11_2023","Q86.12_2023"),
       c24=c("Q76.11_2024","Q76.12_2024"), c25=c("Q74.11_2025","Q74.12_2025"), ndp=2),
  list(id="d7_sleep",    label="Sleep hours/day", domain="7. Health Behavior",
       kind="single", c22="Q28.9_2022", c23="Q31.9_2023", c24="Q28.9_2024", c25="Q32.9_2025", ndp=2,
       is_timeuse=TRUE),
  list(id="d7_sitting",  label="Sitting hours/day", domain="7. Health Behavior",
       kind="single", c22="Q28.5_2022", c23="Q31.5_2023", c24="Q28.5_2024", c25="Q32.5_2025", ndp=2,
       is_timeuse=TRUE),
  list(id="d7_walking",  label="Walking hours/day", domain="7. Health Behavior",
       kind="single", c22="Q28.6_2022", c23="Q31.6_2023", c24="Q28.6_2024", c25="Q32.6_2025", ndp=2,
       is_timeuse=TRUE)
)

.tc_build_y <- function(df, spec, wave_key) {
  cols <- spec[[wave_key]]
  if (is.null(cols)) return(rep(NA_real_, nrow(df)))
  if (!all(cols %in% names(df))) return(rep(NA_real_, nrow(df)))
  v <- if (spec$kind == "single") {
    .tc_as_num(df[[cols]])
  } else if (spec$kind == "mean") {
    .tc_row_mean(df, cols)
  } else if (spec$kind == "sum") {
    .tc_row_sum_tx(df, cols)
  } else if (spec$kind == "sum_recode") {
    .tc_row_sum_tx(df, cols, fn = spec$recode)
  } else stop("unknown kind: ", spec$kind)
  # Time-use band → hours conversion. Re-bin 2025 first if
  # needed (16 → 12 levels), then map band to midpoint hours; わからない → NA.
  if (isTRUE(spec$is_timeuse)) {
    if (wave_key == "c25") v <- recode_timeuse_2025(v)
    v <- band_to_hours(v)
  }
  v
}

.tc_derive_confounders <- function(df, suffix) {
  sex_col <- paste0("SEX_",   suffix)
  age_col <- paste0("AGE_",   suffix)
  edu_col <- paste0("Q21.1_", suffix)
  # Bug #1/#2/#3 fix: per-wave switch for employment + income.
  emp_col <- switch(suffix,
                    "2022" = "Q5.1_2022",
                    "2023" = "Q6.1_2023",
                    "2024" = "Q5.1_2024",
                    "2025" = "Q5.1_2025",
                    paste0("Q5.1_", suffix))
  inc_col <- switch(suffix,
                    "2022" = "Q87.1_2022",
                    "2023" = "Q90.1_2023",
                    "2024" = "Q80.1_2024",
                    "2025" = "Q77.1_2025",
                    paste0("Q87.1_", suffix))
  mar_col <- paste0("Q2_",    suffix)
  liv_col <- paste0("Q1.1_",  suffix)
  sp_col  <- if (suffix == "2022") "Q28.13_2022" else
             if (suffix == "2023") "Q31.13_2023" else
             if (suffix == "2024") "Q28.13_2024" else NA_character_
  pc_col  <- if (suffix == "2022") "Q28.14_2022" else
             if (suffix == "2023") "Q31.14_2023" else
             if (suffix == "2024") "Q28.14_2024" else NA_character_

  out <- list()
  out$is_female <- as.integer(.tc_as_num(df[[sex_col]]) == 2L)
  out$age_raw   <- .tc_as_num(df[[age_col]])
  age_v <- out$age_raw
  age_bands <- list(c(18,24), c(25,29), c(30,34), c(35,39), c(40,44),
                    c(45,49), c(50,54), c(55,59), c(60,64))
  for (b in age_bands)
    out[[sprintf("age_%d_%d", b[1], b[2])]] <-
      as.integer(age_v >= b[1] & age_v <= b[2])
  out$age_65_plus <- as.integer(age_v >= 65)

  # Education: drop Unknown dummy; combine Unknown + NA into
  # the "less than university" reference. Two non-reference dummies remain.
  edu <- .tc_as_num(df[[edu_col]])
  out$edu_univ <- as.integer(edu %in% 6:8)
  out$edu_grad <- as.integer(edu == 9)
  out$edu_hs_or_less <- as.integer(!(out$edu_univ | out$edu_grad))

  emp <- .tc_as_num(df[[emp_col]])
  out$emp_exec    <- as.integer(emp == 1)
  out$emp_self    <- as.integer(emp %in% 2:4)
  out$emp_nonreg  <- as.integer(emp %in% 7:11)
  out$emp_student <- as.integer(emp %in% 12:13)
  out$emp_notwork <- as.integer(emp %in% 14:16 | is.na(emp))
  out$emp_regular <- as.integer(!(out$emp_exec | out$emp_self | out$emp_nonreg |
                                  out$emp_student | out$emp_notwork))

  # Income: fixed yen-band cuts replace quantile-on-discrete-
  # ordinal cut. <2M = codes 1-4; 2-6M = 5-8; 6-10M = 9-12; 10M+ = 13-18;
  # Unknown = 19, 20, NA. Reference = <2M.
  inc <- .tc_as_num(df[[inc_col]])
  inc_band <- ifelse(inc %in% 1:4,   1L,
              ifelse(inc %in% 5:8,   2L,
              ifelse(inc %in% 9:12,  3L,
              ifelse(inc %in% 13:18, 4L, NA_integer_))))
  out$inc_lt2m     <- as.integer(!is.na(inc_band) & inc_band == 1L)
  out$inc_2_6m     <- as.integer(!is.na(inc_band) & inc_band == 2L)
  out$inc_6_10m    <- as.integer(!is.na(inc_band) & inc_band == 3L)
  out$inc_10m_plus <- as.integer(!is.na(inc_band) & inc_band == 4L)
  out$inc_unknown  <- as.integer(is.na(inc_band))

  out$married      <- as.integer(.tc_as_num(df[[mar_col]]) %in% 1:3)
  out$living_alone <- as.integer(!is.na(.tc_as_num(df[[liv_col]])) &
                                 .tc_as_num(df[[liv_col]]) == 1L)
  # Smartphone / PC time: convert band → hours.
  # わからない (band 12) → NA. NAs are KEPT here (unlike the analysis-side
  # derive_confounders which fills to 0 for fill_zero_pre): for Table 1 display,
  # apply_cont/apply_cat already filter is.finite, so わからない responders are
  # automatically excluded from the displayed mean (#4).
  sp_raw <- if (!is.na(sp_col)) .tc_as_num(df[[sp_col]]) else rep(NA_real_, nrow(df))
  out$smartphone <- band_to_hours(sp_raw)
  pc_raw <- if (!is.na(pc_col)) .tc_as_num(df[[pc_col]]) else rep(NA_real_, nrow(df))
  out$pc_tablet  <- band_to_hours(pc_raw)
  as.data.frame(out, stringsAsFactors = FALSE)
}

.tc_fmt_meansd <- function(v, ndp = 2) {
  v <- v[is.finite(v)]
  if (length(v) < 1L) return("--")
  sprintf(paste0("%.", ndp, "f (%.", ndp, "f)"), mean(v), sd(v))
}
.tc_fmt_npct <- function(v) {
  v <- as.integer(v)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return("0 (0.0%)")
  n <- sum(v == 1L); pct <- 100 * mean(v == 1L)
  sprintf("%d (%.1f%%)", n, pct)
}
.tc_fmt_p <- function(p) {
  if (!is.finite(p)) return("--")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

.tc_pval_cont <- function(v, grp) {
  v <- as.numeric(v); grp <- as.integer(grp)
  ok <- is.finite(v) & !is.na(grp)
  if (sum(ok) < 4L || length(unique(grp[ok])) < 2L) return(NA_real_)
  tryCatch(t.test(v[ok] ~ grp[ok])$p.value, error = function(e) NA_real_)
}
.tc_pval_cat <- function(v, grp) {
  v <- as.integer(v); grp <- as.integer(grp)
  ok <- !is.na(v) & !is.na(grp)
  if (sum(ok) < 4L) return(NA_real_)
  tab <- table(v[ok], grp[ok])
  if (any(dim(tab) < 2L)) return(NA_real_)
  tryCatch(suppressWarnings(chisq.test(tab)$p.value), error = function(e) NA_real_)
}

.tc_build_one_table <- function(df, baseline_suffix, stratum_idx_list,
                                pval_grp = NULL,
                                exposure_strata_indices = integer(0)) {
  conf_df <- .tc_derive_confounders(df, baseline_suffix)
  wave_key <- paste0("c", substr(baseline_suffix, 3, 4))

  show_pval <- !is.null(pval_grp)

  rows <- list()
  add_row <- function(label, vals, pval = "") {
    out <- c(list(Variable = label), as.list(vals))
    if (show_pval) out[["p-value"]] <- if (length(pval) == 0L) "" else pval
    rows[[length(rows) + 1L]] <<- as.data.frame(out,
                                                stringsAsFactors = FALSE,
                                                check.names = FALSE)
  }
  add_subheader <- function(label) {
    blanks <- setNames(rep("", length(stratum_idx_list)),
                       names(stratum_idx_list))
    out <- c(list(Variable = label), as.list(blanks))
    if (show_pval) out[["p-value"]] <- ""
    rows[[length(rows) + 1L]] <<- as.data.frame(out,
                                                stringsAsFactors = FALSE,
                                                check.names = FALSE)
  }

  apply_cont <- function(v, ndp = 2) {
    vapply(stratum_idx_list, function(ix) .tc_fmt_meansd(v[ix], ndp = ndp),
           character(1))
  }
  apply_cat <- function(v) {
    vapply(stratum_idx_list, function(ix) .tc_fmt_npct(v[ix]),
           character(1))
  }
  pval_cont <- function(v) {
    if (!show_pval) return("")
    .tc_fmt_p(.tc_pval_cont(v, pval_grp))
  }
  pval_cat <- function(v) {
    if (!show_pval) return("")
    .tc_fmt_p(.tc_pval_cat(v, pval_grp))
  }

  n_vals <- vapply(stratum_idx_list, function(ix) format(length(ix), big.mark=","),
                   character(1))
  add_row("n", n_vals, pval = "")

  add_subheader("[Sex]")
  add_row("    Female, n (%)", apply_cat(conf_df$is_female),
          pval = pval_cat(conf_df$is_female))
  add_row("    Male, n (%)",   apply_cat(1L - conf_df$is_female), pval = "")

  add_subheader("[Age]")
  add_row("    Age, mean (SD)", apply_cont(conf_df$age_raw, 1),
          pval = pval_cont(conf_df$age_raw))
  for (lab in c("18-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64")) {
    cn <- paste0("age_", sub("-", "_", lab))
    add_row(sprintf("    %s, n (%%)", lab), apply_cat(conf_df[[cn]]), pval = "")
  }
  add_row("    >=65 (reference), n (%)", apply_cat(conf_df$age_65_plus), pval = "")

  add_subheader("[Education]")
  # Per : Unknown is folded into the reference (less than university).
  add_row("    < University, incl. Unknown (reference), n (%)",
          apply_cat(conf_df$edu_hs_or_less),
          pval = pval_cat(conf_df$edu_univ + 2L*conf_df$edu_grad))
  add_row("    University, n (%)", apply_cat(conf_df$edu_univ), pval = "")
  add_row("    Graduate, n (%)",   apply_cat(conf_df$edu_grad), pval = "")

  add_subheader("[Employment]")
  emp_full <- 1L*conf_df$emp_exec + 2L*conf_df$emp_self + 3L*conf_df$emp_nonreg +
              4L*conf_df$emp_student + 5L*conf_df$emp_notwork
  add_row("    Regular employee (reference), n (%)", apply_cat(conf_df$emp_regular),
          pval = pval_cat(emp_full))
  add_row("    Executive, n (%)",     apply_cat(conf_df$emp_exec),    pval = "")
  add_row("    Self-employed, n (%)", apply_cat(conf_df$emp_self),    pval = "")
  add_row("    Non-regular, n (%)",   apply_cat(conf_df$emp_nonreg),  pval = "")
  add_row("    Student, n (%)",       apply_cat(conf_df$emp_student), pval = "")
  add_row("    Not working, n (%)",   apply_cat(conf_df$emp_notwork), pval = "")

  add_subheader("[Household income (annual)]")
  # Per : fixed yen-band cuts (replaces uneven quartile cuts).
  inc_full <- 1L*conf_df$inc_lt2m + 2L*conf_df$inc_2_6m + 3L*conf_df$inc_6_10m +
              4L*conf_df$inc_10m_plus + 5L*conf_df$inc_unknown
  add_row("    < 2,000,000 yen (reference), n (%)", apply_cat(conf_df$inc_lt2m),
          pval = pval_cat(inc_full))
  add_row("    2,000,000-5,999,999 yen, n (%)",  apply_cat(conf_df$inc_2_6m),     pval = "")
  add_row("    6,000,000-9,999,999 yen, n (%)",  apply_cat(conf_df$inc_6_10m),    pval = "")
  add_row("    >= 10,000,000 yen, n (%)",        apply_cat(conf_df$inc_10m_plus), pval = "")
  add_row("    Unknown, n (%)",                  apply_cat(conf_df$inc_unknown),  pval = "")

  add_subheader("[Household]")
  add_row("    Married/cohabiting, n (%)", apply_cat(conf_df$married),
          pval = pval_cat(conf_df$married))
  add_row("    Living alone, n (%)", apply_cat(conf_df$living_alone),
          pval = pval_cat(conf_df$living_alone))

  add_subheader(sprintf("[%s screen-time (hours/day, confounders)]", baseline_suffix))
  # Per : time-use bands converted to hours/day (midpoints; わからない
  # → NA, excluded from displayed mean — see footnote in supplement).
  add_row("    Smartphone time (hours/day), mean (SD)",
          apply_cont(conf_df$smartphone, 2),
          pval = pval_cont(conf_df$smartphone))
  add_row("    PC/tablet time (hours/day), mean (SD)",
          apply_cont(conf_df$pc_tablet, 2),
          pval = pval_cont(conf_df$pc_tablet))

  prev_dom <- NA_character_
  for (sp in .tc_outcome_specs) {
    if (!identical(sp$domain, prev_dom)) {
      add_subheader(sprintf("[Baseline %s] %s", baseline_suffix, sp$domain))
      prev_dom <- sp$domain
    }
    v <- .tc_build_y(df, sp, wave_key)
    add_row(sprintf("    %s, mean (SD)", sp$label),
            apply_cont(v, sp$ndp),
            pval = pval_cont(v))
  }

  out_df <- do.call(rbind, rows)
  out_df
}

.tc_csv_to_md <- function(df_tbl) {
  cn <- colnames(df_tbl)
  esc <- function(x) gsub("\\|", "\\\\|", x)
  hdr <- paste0("| ", paste(esc(cn), collapse = " | "), " |")
  sep <- paste0("|", paste(rep(" --- ", length(cn)), collapse = "|"), "|")
  body <- apply(df_tbl, 1, function(row) {
    paste0("| ", paste(esc(as.character(row)), collapse = " | "), " |")
  })
  c(hdr, sep, body)
}

.tc_save_table <- function(df_tbl, basename, label) {
  csv_path <- file.path(.tc_out_dir, paste0(basename, ".csv"))
  md_path  <- file.path(.tc_out_dir, paste0(basename, ".md"))
  write.csv(df_tbl, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(c(paste0("# ", label), "", .tc_csv_to_md(df_tbl)), md_path,
             useBytes = TRUE)
  cat(sprintf("[%s] saved %s  (%d rows x %d cols)\n",
              label, csv_path, nrow(df_tbl), ncol(df_tbl)))
}

.tc_path_3wave   <- PANEL_3W_MAIN
.tc_path_3wave22 <- PANEL_3W_SENS1
.tc_path_2wave   <- PANEL_2W

# ---------- Table 1: Main (codes 5+6 vs 1; baseline 2024) ---------------
if (file.exists(.tc_path_3wave)) {
  df_main <- readr::read_csv(.tc_path_3wave, show_col_types = FALSE)
  ai <- .tc_as_num(df_main$Q37S1_2025)
  keep <- which(ai %in% c(1L, 5L, 6L))
  df_m <- df_main[keep, , drop = FALSE]
  ai2 <- ai[keep]
  is_user <- as.integer(ai2 %in% c(5L, 6L))

  strata <- list(
    Overall          = seq_len(nrow(df_m)),
    `Never users`    = which(ai2 == 1L),
    Users            = which(ai2 %in% c(5L, 6L))
  )
  tbl1 <- .tc_build_one_table(df_m, baseline_suffix = "2024",
                              stratum_idx_list = strata,
                              pval_grp = is_user)
  .tc_save_table(tbl1, "table1_main", "Table 1")

  # Sup Table 1A: Main + 9 purposes (any-use)
  q37s3_cols <- paste0("Q37S3.", 1:9, "_2025")
  purpose_labels_pub <- c(
    "Q37S3.1 Writing/composition",
    "Q37S3.2 Translation",
    "Q37S3.3 Lookup/reference",
    "Q37S3.4 Image/video gen",
    "Q37S3.5 Learning/study",
    "Q37S3.6 Daily-life planning",
    "Q37S3.7 Health information",
    "Q37S3.8 Conversation",
    "Q37S3.9 Emotional support"
  )
  user_rows <- ai2 %in% c(5L, 6L)
  any_user_idx <- function(threshold_test) {
    out <- list()
    for (k in seq_len(9)) {
      v <- .tc_as_num(df_m[[q37s3_cols[k]]])
      out[[purpose_labels_pub[k]]] <- which(user_rows & threshold_test(v))
    }
    out
  }
  strata_1a <- c(
    list(`Never users` = which(ai2 == 1L),
         `Users (any)` = which(ai2 %in% c(5L, 6L))),
    any_user_idx(function(v) is.finite(v) & v >= 2)
  )
  tbl1a <- .tc_build_one_table(df_m, baseline_suffix = "2024",
                               stratum_idx_list = strata_1a,
                               pval_grp = NULL)
  .tc_save_table(tbl1a, "sup_table_1a_main_purposes_anyuse", "Sup Table 1A")

  # Sup Table 1B: Main + 9 purposes (daily-use)
  strata_1b <- c(
    list(`Never users` = which(ai2 == 1L),
         `Users (any)` = which(ai2 %in% c(5L, 6L))),
    any_user_idx(function(v) is.finite(v) & v == 5)
  )
  tbl1b <- .tc_build_one_table(df_m, baseline_suffix = "2024",
                               stratum_idx_list = strata_1b,
                               pval_grp = NULL)
  .tc_save_table(tbl1b, "sup_table_1b_main_purposes_daily", "Sup Table 1B")
} else {
  cat(sprintf("[Table 1] skipped: %s not found\n", .tc_path_3wave))
}

# ---------- Sup Table 2: Sens 1 characteristics (codes 4+5+6 vs 1; baseline 2023)
if (file.exists(.tc_path_3wave22)) {
  df_s1 <- readr::read_csv(.tc_path_3wave22, show_col_types = FALSE)
  ai <- .tc_as_num(df_s1$Q37S1_2025)
  keep <- which(ai %in% c(1L, 4L, 5L, 6L))
  df_s <- df_s1[keep, , drop = FALSE]
  ai2 <- ai[keep]
  is_user <- as.integer(ai2 %in% c(4L, 5L, 6L))
  strata <- list(
    Overall       = seq_len(nrow(df_s)),
    `Never users` = which(ai2 == 1L),
    Users         = which(ai2 %in% c(4L, 5L, 6L))
  )
  tbl2 <- .tc_build_one_table(df_s, baseline_suffix = "2023",
                              stratum_idx_list = strata,
                              pval_grp = is_user)
  .tc_save_table(tbl2, "sup_table_2_sens1_characteristics", "Sup Table 2")
} else {
  cat(sprintf("[Sup Table 2] skipped: %s not found\n", .tc_path_3wave22))
}

# ---------- Sup Table 3: Sens 3 characteristics (full sample; baseline 2022)
if (file.exists(.tc_path_2wave)) {
  df_s3 <- readr::read_csv(.tc_path_2wave, show_col_types = FALSE)
  ai <- .tc_as_num(df_s3$Q37S1_2025)
  # Per : code 2 (past user) excluded for consistency with other tiers.
  D <- ifelse(ai %in% c(3L, 4L, 5L, 6L), 1L,
       ifelse(ai == 1L,                  0L, NA_integer_))
  keep <- which(!is.na(D))
  df_3 <- df_s3[keep, , drop = FALSE]
  D2 <- D[keep]
  strata <- list(
    Overall                              = seq_len(nrow(df_3)),
    `D=0 (Never users, code 1)`          = which(D2 == 0L),
    `D=1 (Any AI use, codes {3,4,5,6})`  = which(D2 == 1L)
  )
  tbl3 <- .tc_build_one_table(df_3, baseline_suffix = "2022",
                              stratum_idx_list = strata,
                              pval_grp = D2)
  .tc_save_table(tbl3, "sup_table_3_sens3_characteristics", "Sup Table 3")
} else {
  cat(sprintf("[Sup Table 3] skipped: %s not found\n", .tc_path_2wave))
}

cat("[characteristics tables] done\n")

# ----------------------------------------------------------------------------
# 12.2  Results tables (Table 2 Main headline, Table 3 cross-tier ladder)
# ----------------------------------------------------------------------------

.tr_outcome_order <- c(
  "d1_pm_avg","d1_k6",
  "d2_hls_avg",
  "d3_mp_avg","d3_ikigai",
  "d4_cv_avg",
  "d5_rel_avg","d5_lsns_total","d5_ucla",
  "d6_worry_avg",
  "d7_sleep","d7_sitting","d7_walking"
)
.tr_outcome_labels <- c(
  d1_pm_avg     = "Physical & mental health",
  d1_k6         = "K6 distress (sum 0-24)",
  d2_hls_avg    = "Happiness & life satisfaction",
  d3_mp_avg     = "Meaning & purpose",
  d3_ikigai     = "Ikigai",
  d4_cv_avg     = "Character & virtue",
  d5_rel_avg    = "Relations",
  d5_lsns_total = "LSNS-6 total (0-30)",
  d5_ucla       = "UCLA-3 loneliness",
  d6_worry_avg  = "Financial & safety worry",
  d7_sleep      = "Sleep hours",
  d7_sitting    = "Sitting time",
  d7_walking    = "Walking time"
)
.tr_domain_of <- c(
  d1_pm_avg     = "D1. Mental & Physical Health",
  d1_k6         = "D1. Mental & Physical Health",
  d2_hls_avg    = "D2. Happiness & Life Satisfaction",
  d3_mp_avg     = "D3. Meaning & Purpose",
  d3_ikigai     = "D3. Meaning & Purpose",
  d4_cv_avg     = "D4. Character & Virtue",
  d5_rel_avg    = "D5. Social Wellbeing",
  d5_lsns_total = "D5. Social Wellbeing",
  d5_ucla       = "D5. Social Wellbeing",
  d6_worry_avg  = "D6. Financial & Material Stability",
  d7_sleep      = "D7. Health Behavior",
  d7_sitting    = "D7. Health Behavior",
  d7_walking    = "D7. Health Behavior"
)
.tr_pattern_label <- c(
  use       = "A0 use",
  intensity = "A intensity",
  creative  = "B creative",
  daily     = "B daily",
  social    = "B social"
)
.tr_pattern_order <- c("use", "intensity", "creative", "daily", "social")

.tr_make_md_table <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) {
    r[is.na(r)] <- ""
    paste0("| ", paste(r, collapse = " | "), " |")
  })
  c(hdr, sep, rows)
}

# Read a tier's combined results CSV (one file per tier per phase).
read_tier_results <- function(top_dir, tier) {
  stopifnot(tier %in% c("main", "sens1", "sens2", "sens3"))
  tf_dir <- file.path(top_dir, "manuscript", "tables_figures", tier)

  path <- switch(tier,
    main  = file.path(tf_dir, "results_main_did.csv"),
    sens1 = file.path(tf_dir, "results_sens1_did.csv"),
    sens2 = file.path(tf_dir, "results_sens2_ancova.csv"),
    sens3 = file.path(tf_dir, "results_sens3_ancova.csv")
  )

  if (!file.exists(path)) {
    warning(sprintf("[read_tier_results] missing: %s", path))
    return(NULL)
  }

  raw <- utils::read.csv(path, stringsAsFactors = FALSE)
  raw <- raw[raw$exposure %in% .tr_pattern_order, , drop = FALSE]
  raw$bh_sig <- as.logical(raw$bh_sig)
  raw$tier <- tier

  raw$outcome  <- factor(raw$outcome,  levels = .tr_outcome_order)
  raw$exposure <- factor(raw$exposure, levels = .tr_pattern_order)
  raw <- raw[order(raw$outcome, raw$exposure), , drop = FALSE]
  raw$outcome  <- as.character(raw$outcome)
  raw$exposure <- as.character(raw$exposure)

  raw
}

.tr_fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}
.tr_fmt_q <- function(x) {
  ifelse(is.na(x), "",
    ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3)))
}

build_table2 <- function(top_dir) {
  d <- read_tier_results(top_dir, "main")
  if (is.null(d) || !nrow(d)) {
    warning("[build_table2] no Main-tier rows available")
    return(invisible(NULL))
  }

  e_val    <- if ("e_value"    %in% names(d)) d$e_value    else rep(NA_real_, nrow(d))
  e_val_ci <- if ("e_value_ci" %in% names(d)) d$e_value_ci else rep(NA_real_, nrow(d))
  fmt2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 2))
  out <- data.frame(
    Domain  = unname(.tr_domain_of[d$outcome]),
    Outcome = unname(.tr_outcome_labels[d$outcome]),
    Pattern = unname(.tr_pattern_label[d$exposure]),
    `beta_oz [95% CI]` = sprintf(
      "%s [%s, %s]",
      .tr_fmt_num(d$estimate_oz),
      .tr_fmt_num(d$ci_lo_oz),
      .tr_fmt_num(d$ci_hi_oz)
    ),
    `q (BH)`  = .tr_fmt_q(d$p_bh),
    `BH-sig`  = ifelse(d$bh_sig %in% TRUE, "yes", ""),
    `E-value` = fmt2(e_val),
    `E-value (CI)` = fmt2(e_val_ci),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  out$.do_idx <- match(d$outcome,  .tr_outcome_order)
  out$.pa_idx <- match(d$exposure, .tr_pattern_order)
  out <- out[order(out$.do_idx, out$.pa_idx), , drop = FALSE]
  out$.do_idx <- NULL
  out$.pa_idx <- NULL
  rownames(out) <- NULL

  out_dir <- file.path(top_dir, "manuscript")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  csv_path <- file.path(out_dir, "table2_main_results.csv")
  md_path  <- file.path(out_dir, "table2_main_results.md")
  utils::write.csv(out, csv_path, row.names = FALSE, fileEncoding = "UTF-8")

  md_lines <- c(
    "# Table 2 - Main DiD headline results",
    "",
    "Tier: **Main** (3-wave DiD; Delta Y 2024 -> 2025; codes 5+6 vs 1).",
    "Coefficients are outcome-standardised (beta_oz = beta / SD(Y_2024)).",
    "BH-sig = q < 0.10 within the 13-test family for that exposure.",
    "E-value: minimum strength of unmeasured confounder (on the risk-ratio scale)",
    "needed to fully explain away the observed association",
    "(VanderWeele & Ding 2017; RR approx via exp(0.91 x |d|), Chinn 2000).",
    "E-value (CI) = E-value evaluated at the CI bound closer to the null;",
    "= 1.00 when the 95% CI crosses zero (no confounder needed).",
    "",
    .tr_make_md_table(out)
  )
  writeLines(md_lines, md_path, useBytes = TRUE)

  cat(sprintf("[Table 2] saved %s and .md (%d rows)\n", csv_path, nrow(out)))
  invisible(out)
}

build_table3 <- function(top_dir) {
  tiers <- c("main", "sens1", "sens2", "sens3")
  reads <- lapply(tiers, function(t) read_tier_results(top_dir, t))
  names(reads) <- tiers

  symbol_for <- function(df, oc, ex) {
    if (is.null(df)) return("?")
    row <- df[df$outcome == oc & df$exposure == ex, , drop = FALSE]
    if (!nrow(row)) return("?")
    if (isTRUE(row$bh_sig[1])) {
      if (!is.na(row$estimate_oz[1]) && row$estimate_oz[1] >= 0) "+" else "-"
    } else {
      "o"
    }
  }

  cells <- matrix("", nrow = length(.tr_outcome_order),
                  ncol = length(.tr_pattern_order),
                  dimnames = list(.tr_outcome_order, .tr_pattern_order))
  for (oc in .tr_outcome_order) {
    for (ex in .tr_pattern_order) {
      s <- vapply(tiers, function(t) symbol_for(reads[[t]], oc, ex), character(1))
      cells[oc, ex] <- paste(s, collapse = "")
    }
  }

  out <- data.frame(
    Outcome   = unname(.tr_outcome_labels[.tr_outcome_order]),
    Use       = cells[, "use"],
    Intensity = cells[, "intensity"],
    Creative  = cells[, "creative"],
    Daily     = cells[, "daily"],
    Social    = cells[, "social"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out_dir <- file.path(top_dir, "manuscript")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  csv_path <- file.path(out_dir, "table3_cross_tier_ladder.csv")
  md_path  <- file.path(out_dir, "table3_cross_tier_ladder.md")
  utils::write.csv(out, csv_path, row.names = FALSE, fileEncoding = "UTF-8")

  legend <- c(
    "",
    "**Symbol legend** - each cell is a 4-character string with one symbol per tier",
    "in order **Main / Sens 1 / Sens 2 / Sens 3**:",
    "",
    "- `+` BH-significant (q < 0.10 in that tier) **with positive beta_oz**",
    "- `-` BH-significant **with negative beta_oz**",
    "- `o` not BH-significant in that tier",
    "- `?` row missing in that tier's CSV",
    "",
    "Tiers: **Main** = 3-wave DiD Delta Y 2024->2025 (codes 5+6 vs 1); ",
    "**Sens 1** = 3-wave DiD Delta Y 2023->2025 (codes 4+5+6 vs 1); ",
    "**Sens 2** = ANCOVA Y_2025 ~ exposure + Y_2024 + cov (Main cohort); ",
    "**Sens 3** = ANCOVA Y_2025 ~ exposure + Y_2022 + cov (full sample)."
  )
  md_lines <- c(
    "# Table 3 - Cross-tier replication ladder",
    "",
    .tr_make_md_table(out),
    legend
  )
  writeLines(md_lines, md_path, useBytes = TRUE)

  cat(sprintf("[Table 3] saved %s and .md (%d rows x %d cols)\n",
              csv_path, nrow(out), ncol(out)))
  invisible(out)
}

local({
  top_dir <- OUTPUT_DIR
  cat(sprintf("[tables_results] OUTPUT_DIR=%s\n", top_dir))
  build_table2(top_dir)
  build_table3(top_dir)
})

cat("[publication tables] done\n")

# ============================================================================
# 13.3  Main-cohort raw outcome levels and changes (#1).
#       Skipped if Main DiD did not run.
#
#   Sup Table 4 — raw Y_2024 / Y_2025 / DeltaY by Overall / Non-users / Users
#                 with t-test p (Non-users vs Users), all 13 outcomes.
# ============================================================================

if (exists("main_did_ok") && isTRUE(main_did_ok)) {

  suppressPackageStartupMessages({
    library(sandwich); library(lmtest); library(car); library(ggplot2)
  })

  # --- Resolve the Main-cohort objects stashed at end of §7 -----------------
  .sg_ds          <- main_ds
  .sg_y24         <- main_y_2024_df
  .sg_y25         <- main_y_2025_df
  .sg_pre24       <- main_pre_2024_df
  .sg_conf24      <- main_conf_2024_df
  .sg_sd_pre      <- main_sd_pre_did
  .sg_n_total     <- nrow(.sg_ds)

  cat(sprintf("[subgroup] Main cohort n=%d (treat=%d, control=%d)\n",
              .sg_n_total, sum(.sg_ds$treat == 1L), sum(.sg_ds$treat == 0L)))

  # --- Helpers ---------------------------------------------------------------

  .sg_outcome_labels <- outcome_labels   # alias to top-of-script labels
  .sg_outcome_order  <- outcome_order

  .sg_fmt_meansd <- function(v, ndp = 2) {
    v <- v[is.finite(v)]
    if (length(v) < 1L) return("--")
    sprintf(paste0("%.", ndp, "f (%.", ndp, "f)"), mean(v), sd(v))
  }
  .sg_fmt_p <- function(p) {
    if (!is.finite(p)) return("--")
    if (p < 0.001) return("<0.001")
    sprintf("%.3f", p)
  }
  .sg_pval_t <- function(v, grp) {
    v <- as.numeric(v); grp <- as.integer(grp)
    ok <- is.finite(v) & !is.na(grp)
    if (sum(ok) < 4L || length(unique(grp[ok])) < 2L) return(NA_real_)
    tryCatch(t.test(v[ok] ~ grp[ok])$p.value, error = function(e) NA_real_)
  }

  .sg_csv_to_md <- function(df_tbl) {
    cn <- colnames(df_tbl)
    esc <- function(x) gsub("\\|", "\\\\|", x)
    hdr <- paste0("| ", paste(esc(cn), collapse = " | "), " |")
    sep <- paste0("|", paste(rep(" --- ", length(cn)), collapse = "|"), "|")
    body <- apply(df_tbl, 1, function(row)
      paste0("| ", paste(esc(as.character(row)), collapse = " | "), " |"))
    c(hdr, sep, body)
  }
  .sg_save_table <- function(df_tbl, basename, label) {
    csv_path <- file.path(MS_DIR, paste0(basename, ".csv"))
    md_path  <- file.path(MS_DIR, paste0(basename, ".md"))
    write.csv(df_tbl, csv_path, row.names = FALSE, fileEncoding = "UTF-8")
    writeLines(c(paste0("# ", label), "", .sg_csv_to_md(df_tbl)), md_path,
               useBytes = TRUE)
    cat(sprintf("[%s] saved %s  (%d rows x %d cols)\n",
                label, csv_path, nrow(df_tbl), ncol(df_tbl)))
  }

  # ============================================================================
  # Sup Table 4 — raw outcome levels and changes (Main cohort)
  # ============================================================================

  .sg_grp <- as.integer(.sg_ds$treat)            # 0 = Non-users, 1 = Users
  .sg_idx_overall   <- seq_len(.sg_n_total)
  .sg_idx_nonusers  <- which(.sg_grp == 0L)
  .sg_idx_users     <- which(.sg_grp == 1L)
  .sg_n_nonusers    <- length(.sg_idx_nonusers)
  .sg_n_users       <- length(.sg_idx_users)

  .sg_t4_rows <- list()
  .sg_t4_add <- function(varname, overall_str, nonuser_str, user_str, pstr) {
    .sg_t4_rows[[length(.sg_t4_rows) + 1L]] <<- data.frame(
      Variable = varname,
      Overall  = overall_str,
      `Non-users` = nonuser_str,
      Users    = user_str,
      `p-value` = pstr,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  .sg_t4_add("n",
             format(.sg_n_total,    big.mark = ","),
             format(.sg_n_nonusers, big.mark = ","),
             format(.sg_n_users,    big.mark = ","),
             "")

  for (oid in .sg_outcome_order) {
    lab <- unname(.sg_outcome_labels[oid])
    y_pre  <- as.numeric(.sg_y24[[oid]])
    y_post <- as.numeric(.sg_y25[[oid]])
    y_dlt  <- y_post - y_pre   # NA if either side NA

    .sg_t4_add(sprintf("[%s]", lab), "", "", "", "")
    .sg_t4_add("    Y_2024 (baseline), mean (SD)",
               .sg_fmt_meansd(y_pre[.sg_idx_overall]),
               .sg_fmt_meansd(y_pre[.sg_idx_nonusers]),
               .sg_fmt_meansd(y_pre[.sg_idx_users]),
               .sg_fmt_p(.sg_pval_t(y_pre, .sg_grp)))
    .sg_t4_add("    Y_2025 (post), mean (SD)",
               .sg_fmt_meansd(y_post[.sg_idx_overall]),
               .sg_fmt_meansd(y_post[.sg_idx_nonusers]),
               .sg_fmt_meansd(y_post[.sg_idx_users]),
               .sg_fmt_p(.sg_pval_t(y_post, .sg_grp)))
    .sg_t4_add("    DeltaY = Y_2025 - Y_2024, mean (SD)",
               .sg_fmt_meansd(y_dlt[.sg_idx_overall]),
               .sg_fmt_meansd(y_dlt[.sg_idx_nonusers]),
               .sg_fmt_meansd(y_dlt[.sg_idx_users]),
               .sg_fmt_p(.sg_pval_t(y_dlt, .sg_grp)))
  }
  sup_table_4 <- do.call(rbind, .sg_t4_rows)
  .sg_save_table(sup_table_4, "sup_table_4_main_raw_levels_changes",
                 "Sup Table 4")

  cat("[Sup Table 4] [done]\n")

} else {
  cat("[Sup Table 4] skipped: Main DiD did not run\n")
}
