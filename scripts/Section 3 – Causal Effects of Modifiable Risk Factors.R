############################################################
# Section 3 – Causal Effects of Modifiable Risk Factors
# Dataset: PPD_dataset_v2.csv  (Bangladesh postpartum depression)
############################################################

## 1. Load packages ----
# Run once if needed:
# install.packages(c("tidyverse", "tmle", "SuperLearner", "gt"))

library(tidyverse)
library(tmle)
library(SuperLearner)
library(gt)

set.seed(20251129)

## 2. Set up output folders (Section 3) ----
out_root      <- "output"
out_section3  <- file.path(out_root, "section3_prediction")
out_s3_tables <- file.path(out_section3, "tables")
out_s3_figs   <- file.path(out_section3, "figures")

dir.create(out_root,      showWarnings = FALSE)
dir.create(out_section3,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s3_tables, recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s3_figs,   recursive    = TRUE, showWarnings = FALSE)

## 3. Read data & define outcome + exposures (original names) ----
data_path <- "PPD_dataset_v2.csv"

ppd_raw <- readr::read_csv(
  data_path,
  show_col_types = FALSE
)

ppd <- ppd_raw %>%
  mutate(
    across(where(is.character), as.factor),
    epds_score  = as.numeric(`EPDS Score`),
    ppd_epds13  = dplyr::if_else(epds_score >= 13, 1L, 0L)
  ) %>%
  mutate(
    abuse_any         = dplyr::if_else(Abuse == "Yes", 1L, 0L),
    low_rec_support   = dplyr::if_else(`Recieved Support` == "Low", 1L, 0L),
    high_need_support = dplyr::if_else(`Need for Support` %in% c("High", "Medium"), 1L, 0L),
    preg_disease_any  = dplyr::if_else(`Diseases during pregnancy` != "None", 1L, 0L),
    birth_complication = dplyr::if_else(`Birth compliancy` == "Yes", 1L, 0L),
    newborn_illness    = dplyr::if_else(`Newborn illness` == "Yes", 1L, 0L)
  )

outcome_var_orig <- "ppd_epds13"

## 4. Confounder set (original names) ----
w_vars_orig <- c(
  "Age",
  "Residence",
  "Education Level",
  "Marital status",
  "Occupation before latest pregnancy",
  "Monthly income before latest pregnancy",
  "Current monthly income",
  "Husband's education level",
  "Husband’s monthly income",
  "Total children",
  "Disease before pregnancy",
  "History of pregnancy loss",
  "Family type",
  "Number of household members",
  "Major changes or losses during pregnancy",
  "Pregnancy plan",
  "Regular checkups",
  "Fear of pregnancy",
  "Depression before pregnancy (PHQ2)",
  "Depression during pregnancy (PHQ2)"
)

# Keep only covariates that actually exist
w_vars_use_orig <- intersect(w_vars_orig, names(ppd))

## 5. Modifiable exposures (original names) ----
exposures_info <- tibble::tribble(
  ~exposure_var_orig, ~exposure_label,
  "abuse_any",         "Intimate partner violence (any abuse)",
  "low_rec_support",   "Low received support",
  "high_need_support", "High/medium need for support",
  "preg_disease_any",  "Any disease during pregnancy",
  "birth_complication","Birth complications",
  "newborn_illness",   "Neonatal illness"
)

stopifnot(all(exposures_info$exposure_var_orig %in% names(ppd)))

## 6. Create a safe-named version of the data for modelling ----
name_map <- tibble(
  orig = names(ppd),
  safe = make.names(names(ppd), unique = TRUE)
)

ppd_mod <- ppd
names(ppd_mod) <- name_map$safe

# Map outcome and covariates to safe names
outcome_var_safe <- name_map$safe[match(outcome_var_orig, name_map$orig)]

w_vars_use_safe <- name_map$safe[match(w_vars_use_orig, name_map$orig)]

# Map exposures to safe names
exposures_info <- exposures_info %>%
  mutate(
    exposure_var_safe = name_map$safe[match(exposure_var_orig, name_map$orig)]
  )

stopifnot(all(!is.na(exposures_info$exposure_var_safe)))

## 7. Helper functions: TMLE and logistic regression on safe-named data ----

sl_lib <- c("SL.mean", "SL.glm")

run_tmle_single <- function(data, exposure_safe, label, outcome_safe, covars_safe, sl_library) {
  vars_needed <- c(outcome_safe, exposure_safe, covars_safe)
  
  dat <- data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na()
  
  if (nrow(dat) == 0L) {
    return(tibble(
      exposure_label         = label,
      exposure_safe          = exposure_safe,
      n_tmle                 = 0L,
      prevalence_exposed_pct = NA_real_,
      risk0_tmle             = NA_real_,
      risk1_tmle             = NA_real_,
      rd_tmle                = NA_real_,
      rd_tmle_l              = NA_real_,
      rd_tmle_u              = NA_real_,
      rr_tmle                = NA_real_,
      rr_tmle_l              = NA_real_,
      rr_tmle_u              = NA_real_,
      or_tmle                = NA_real_,
      or_tmle_l              = NA_real_,
      or_tmle_u              = NA_real_,
      paf_tmle               = NA_real_,
      paf_tmle_l             = NA_real_,
      paf_tmle_u             = NA_real_
    ))
  }
  
  Y <- dat[[outcome_safe]]
  A <- dat[[exposure_safe]]
  
  W <- if (length(covars_safe) == 0L) {
    NULL
  } else {
    dat %>%
      dplyr::select(dplyr::all_of(covars_safe)) %>%
      as.data.frame()
  }
  
  if (length(unique(A)) < 2L) {
    p_exp <- mean(A == 1, na.rm = TRUE)
    return(tibble(
      exposure_label         = label,
      exposure_safe          = exposure_safe,
      n_tmle                 = nrow(dat),
      prevalence_exposed_pct = 100 * p_exp,
      risk0_tmle             = NA_real_,
      risk1_tmle             = NA_real_,
      rd_tmle                = NA_real_,
      rd_tmle_l              = NA_real_,
      rd_tmle_u              = NA_real_,
      rr_tmle                = NA_real_,
      rr_tmle_l              = NA_real_,
      rr_tmle_u              = NA_real_,
      or_tmle                = NA_real_,
      or_tmle_l              = NA_real_,
      or_tmle_u              = NA_real_,
      paf_tmle               = NA_real_,
      paf_tmle_l             = NA_real_,
      paf_tmle_u             = NA_real_
    ))
  }
  
  paf_fn <- function(rr_val, p) {
    if (is.na(rr_val) || is.na(p)) return(NA_real_)
    p * (rr_val - 1) / (p * (rr_val - 1) + 1)
  }
  
  tm_out <- tryCatch(
    {
      tmle::tmle(
        Y = Y,
        A = A,
        W = W,
        family       = "binomial",
        Q.SL.library = sl_library,
        g.SL.library = sl_library,
        verbose      = FALSE
      )
    },
    error = function(e) {
      message("TMLE failed for ", label, " (", exposure_safe, "): ", conditionMessage(e))
      return(NULL)
    }
  )
  
  p_exp <- mean(A == 1, na.rm = TRUE)
  
  if (is.null(tm_out)) {
    return(tibble(
      exposure_label         = label,
      exposure_safe          = exposure_safe,
      n_tmle                 = nrow(dat),
      prevalence_exposed_pct = 100 * p_exp,
      risk0_tmle             = NA_real_,
      risk1_tmle             = NA_real_,
      rd_tmle                = NA_real_,
      rd_tmle_l              = NA_real_,
      rd_tmle_u              = NA_real_,
      rr_tmle                = NA_real_,
      rr_tmle_l              = NA_real_,
      rr_tmle_u              = NA_real_,
      or_tmle                = NA_real_,
      or_tmle_l              = NA_real_,
      or_tmle_u              = NA_real_,
      paf_tmle               = NA_real_,
      paf_tmle_l             = NA_real_,
      paf_tmle_u             = NA_real_
    ))
  }
  
  risk0 <- unname(tm_out$estimates$EY0$psi)
  risk1 <- unname(tm_out$estimates$EY1$psi)
  
  rd    <- unname(tm_out$estimates$ATE$psi)
  rd_CI <- unname(tm_out$estimates$ATE$CI)
  
  rr    <- unname(tm_out$estimates$RR$psi)
  rr_CI <- unname(tm_out$estimates$RR$CI)
  
  or_v  <- unname(tm_out$estimates$OR$psi)
  or_CI <- unname(tm_out$estimates$OR$CI)
  
  paf   <- paf_fn(rr,       p_exp)
  paf_L <- paf_fn(rr_CI[1], p_exp)
  paf_U <- paf_fn(rr_CI[2], p_exp)
  
  tibble(
    exposure_label         = label,
    exposure_safe          = exposure_safe,
    n_tmle                 = nrow(dat),
    prevalence_exposed_pct = 100 * p_exp,
    risk0_tmle             = risk0,
    risk1_tmle             = risk1,
    rd_tmle                = rd,
    rd_tmle_l              = rd_CI[1],
    rd_tmle_u              = rd_CI[2],
    rr_tmle                = rr,
    rr_tmle_l              = rr_CI[1],
    rr_tmle_u              = rr_CI[2],
    or_tmle                = or_v,
    or_tmle_l              = or_CI[1],
    or_tmle_u              = or_CI[2],
    paf_tmle               = paf,
    paf_tmle_l             = paf_L,
    paf_tmle_u             = paf_U
  )
}

run_glm_single <- function(data, exposure_safe, label, outcome_safe, covars_safe) {
  vars_needed <- c(outcome_safe, exposure_safe, covars_safe)
  
  dat <- data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na()
  
  if (nrow(dat) == 0L) {
    return(tibble(
      exposure_label = label,
      exposure_safe  = exposure_safe,
      n_glm          = 0L,
      risk0_glm      = NA_real_,
      risk1_glm      = NA_real_,
      rd_glm         = NA_real_,
      rr_glm         = NA_real_,
      or_glm         = NA_real_,
      or_glm_l       = NA_real_,
      or_glm_u       = NA_real_
    ))
  }
  
  A <- dat[[exposure_safe]]
  
  if (length(unique(A)) < 2L) {
    return(tibble(
      exposure_label = label,
      exposure_safe  = exposure_safe,
      n_glm          = nrow(dat),
      risk0_glm      = NA_real_,
      risk1_glm      = NA_real_,
      rd_glm         = NA_real_,
      rr_glm         = NA_real_,
      or_glm         = NA_real_,
      or_glm_l       = NA_real_,
      or_glm_u       = NA_real_
    ))
  }
  
  form <- stats::reformulate(
    termlabels = c(exposure_safe, covars_safe),
    response   = outcome_safe
  )
  
  glm_fit <- tryCatch(
    glm(formula = form, family = binomial(), data = dat),
    error = function(e) {
      message("GLM failed for ", label, " (", exposure_safe, "): ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(glm_fit)) {
    return(tibble(
      exposure_label = label,
      exposure_safe  = exposure_safe,
      n_glm          = nrow(dat),
      risk0_glm      = NA_real_,
      risk1_glm      = NA_real_,
      rd_glm         = NA_real_,
      rr_glm         = NA_real_,
      or_glm         = NA_real_,
      or_glm_l       = NA_real_,
      or_glm_u       = NA_real_
    ))
  }
  
  sm <- summary(glm_fit)
  
  if (!exposure_safe %in% rownames(sm$coefficients)) {
    beta <- NA_real_; se <- NA_real_
  } else {
    beta <- sm$coefficients[exposure_safe, "Estimate"]
    se   <- sm$coefficients[exposure_safe, "Std. Error"]
  }
  
  if (is.na(beta) || is.na(se)) {
    or_val <- or_l <- or_u <- NA_real_
  } else {
    or_val <- exp(beta)
    or_l   <- exp(beta - 1.96 * se)
    or_u   <- exp(beta + 1.96 * se)
  }
  
  new1 <- dat
  new0 <- dat
  new1[[exposure_safe]] <- 1
  new0[[exposure_safe]] <- 0
  
  pred1 <- predict(glm_fit, newdata = new1, type = "response")
  pred0 <- predict(glm_fit, newdata = new0, type = "response")
  
  risk1_g <- mean(pred1)
  risk0_g <- mean(pred0)
  rd_g    <- risk1_g - risk0_g
  rr_g    <- risk1_g / risk0_g
  
  tibble(
    exposure_label = label,
    exposure_safe  = exposure_safe,
    n_glm          = nrow(dat),
    risk0_glm      = risk0_g,
    risk1_glm      = risk1_g,
    rd_glm         = rd_g,
    rr_glm         = rr_g,
    or_glm         = or_val,
    or_glm_l       = or_l,
    or_glm_u       = or_u
  )
}

## 8. Run TMLE + logistic regression for all exposures ----

tmle_results <- purrr::map_dfr(
  seq_len(nrow(exposures_info)),
  ~ run_tmle_single(
    data        = ppd_mod,
    exposure_safe = exposures_info$exposure_var_safe[.x],
    label       = exposures_info$exposure_label[.x],
    outcome_safe = outcome_var_safe,
    covars_safe = w_vars_use_safe,
    sl_library  = sl_lib
  )
)

glm_results <- purrr::map_dfr(
  seq_len(nrow(exposures_info)),
  ~ run_glm_single(
    data         = ppd_mod,
    exposure_safe = exposures_info$exposure_var_safe[.x],
    label        = exposures_info$exposure_label[.x],
    outcome_safe = outcome_var_safe,
    covars_safe  = w_vars_use_safe
  )
)

results_combined <- tmle_results %>%
  left_join(
    glm_results,
    by = c("exposure_label", "exposure_safe")
  )

readr::write_csv(
  results_combined,
  file.path(out_s3_tables, "tmle_vs_glm_causal_results_section3.csv")
)

## 9. Table 3 – TMLE vs conventional regression (HTML) ----

table3_display <- results_combined %>%
  mutate(
    prevalence_exposed_pct = round(prevalence_exposed_pct, 1),
    risk0_tmle_pct         = round(100 * risk0_tmle, 1),
    risk1_tmle_pct         = round(100 * risk1_tmle, 1),
    rd_tmle_pp             = round(100 * rd_tmle, 1),
    paf_tmle_pct           = round(100 * paf_tmle, 1),
    rr_tmle_round          = round(rr_tmle, 2),
    or_glm_round           = round(or_glm, 2),
    rr_glm_round           = round(rr_glm, 2)
  ) %>%
  select(
    exposure_label,
    prevalence_exposed_pct,
    risk0_tmle_pct,
    risk1_tmle_pct,
    rd_tmle_pp,
    rr_tmle_round,
    paf_tmle_pct,
    or_glm_round,
    rr_glm_round
  )

table3_gt <-
  table3_display %>%
  gt(rowname_col = NULL) %>%
  tab_header(
    title = md("**Table 3. Causal effects of modifiable risk factors on postpartum depression**"),
    subtitle = md("Targeted maximum likelihood estimation (TMLE) versus conventional logistic regression with g-computation.")
  ) %>%
  cols_label(
    exposure_label          = md("**Exposure (modifiable risk factor)**"),
    prevalence_exposed_pct  = md("**% exposed**"),
    risk0_tmle_pct          = md("**Risk unexposed (TMLE), %**"),
    risk1_tmle_pct          = md("**Risk exposed (TMLE), %**"),
    rd_tmle_pp              = md("**Risk difference (TMLE), percentage points**"),
    rr_tmle_round           = md("**Relative risk (TMLE)**"),
    paf_tmle_pct            = md("**Population-attributable fraction (TMLE), %**"),
    or_glm_round            = md("**Odds ratio (logistic)**"),
    rr_glm_round            = md("**Relative risk (g-computation)**")
  ) %>%
  fmt_number(
    columns = c(
      prevalence_exposed_pct,
      risk0_tmle_pct,
      risk1_tmle_pct,
      rd_tmle_pp,
      paf_tmle_pct
    ),
    decimals = 1
  ) %>%
  fmt_number(
    columns = c(rr_tmle_round, or_glm_round, rr_glm_round),
    decimals = 2
  ) %>%
  tab_spanner(
    label = md("**TMLE (targeted learning)**"),
    columns = c(
      risk0_tmle_pct,
      risk1_tmle_pct,
      rd_tmle_pp,
      rr_tmle_round,
      paf_tmle_pct
    )
  ) %>%
  tab_spanner(
    label = md("**Conventional logistic regression**"),
    columns = c(or_glm_round, rr_glm_round)
  ) %>%
  tab_options(
    table.font.size  = px(12),
    data_row.padding = px(2),
    heading.align    = "center"
  )

gtsave(
  table3_gt,
  filename = file.path(out_s3_tables, "table3_tmle_vs_glm_causal_effects.html")
)

## 10. Figure 2 – Forest plot of TMLE relative risks (PNG) ----

forest_data <- tmle_results %>%
  mutate(
    exposure_label = factor(
      exposure_label,
      levels = rev(exposure_label)
    )
  )

p_forest_tmle <-
  ggplot(forest_data,
         aes(y = exposure_label, x = rr_tmle)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_errorbarh(
    aes(xmin = rr_tmle_l, xmax = rr_tmle_u),
    height = 0.2,
    na.rm  = TRUE
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    labels = c("0.5", "0.75", "1", "1.5", "2", "3")
  ) +
  labs(
    x = "Relative risk (TMLE, log scale)",
    y = NULL,
    title = "Figure 2. Causal effects of modifiable risk factors on postpartum depression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y      = element_text(size = 10),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

p_forest_tmle

ggplot2::ggsave(
  filename = file.path(out_s3_figs, "figure2_tmle_forest_relative_risks.png"),
  plot     = p_forest_tmle,
  width    = 7.5,
  height   = 5,
  dpi      = 300
)

############################################################
# Outputs for Section 3:
# - output/section3_prediction/tables/tmle_vs_glm_causal_results_section3.csv
# - output/section3_prediction/tables/table3_tmle_vs_glm_causal_effects.html
# - output/section3_prediction/figures/figure2_tmle_forest_relative_risks.png
############################################################
