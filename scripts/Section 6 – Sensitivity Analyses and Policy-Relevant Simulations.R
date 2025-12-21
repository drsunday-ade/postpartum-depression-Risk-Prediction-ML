############################################################
# Section 6 – Sensitivity Analyses and Policy-Relevant Simulations
# Dataset: PPD_dataset_v2.csv  (Bangladesh postpartum depression)
############################################################

## 1. Packages (auto-install for reproducibility) -------------------------
required_pkgs <- c("tidyverse", "glmnet", "gt")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)   # dplyr, ggplot2, purrr, tibble, forcats, etc.
library(glmnet)
library(gt)

set.seed(20251129)

## 2. Set up output folders (Section 6) -----------------------------------
out_root      <- "output"
out_section6  <- file.path(out_root, "section6_sensitivity_policy")
out_s6_tables <- file.path(out_section6, "tables")
out_s6_figs   <- file.path(out_section6, "figures")

dir.create(out_root,      showWarnings = FALSE)
dir.create(out_section6,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s6_tables, recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s6_figs,   recursive    = TRUE, showWarnings = FALSE)

## 3. Read data & define outcomes / covariates ----------------------------
data_path <- "PPD_dataset_v2.csv"

ppd_raw <- readr::read_csv(
  data_path,
  show_col_types = FALSE
)

ppd <- ppd_raw %>%
  mutate(
    row_id     = dplyr::row_number(),
    across(where(is.character), as.factor),
    epds_score = as.numeric(`EPDS Score`),
    # Primary outcome
    ppd_epds13 = dplyr::if_else(!is.na(epds_score) & epds_score >= 13, 1L, 0L),
    abuse_any  = dplyr::if_else(Abuse == "Yes", 1L, 0L)
  )

# Alternative outcome 1: EPDS >= 11
ppd <- ppd %>%
  mutate(
    ppd_epds11 = dplyr::if_else(!is.na(epds_score) & epds_score >= 11, 1L, 0L)
  )

# Alternative outcome 2: PHQ-9 >= 10 (if available)
phq9_candidates <- c("PHQ-9 Score", "PHQ9 Score", "PHQ9", "PHQ_9", "PHQ.9")
phq9_name <- intersect(phq9_candidates, names(ppd_raw))

if (length(phq9_name) > 0L) {
  phq9_name <- phq9_name[1]
  ppd <- ppd %>%
    mutate(
      phq9_score  = suppressWarnings(as.numeric(.data[[phq9_name]])),
      ppd_phq9_10 = dplyr::if_else(!is.na(phq9_score) & phq9_score >= 10, 1L, 0L)
    )
  has_phq9 <- TRUE
} else {
  message("Section 6: No PHQ-9 score variable found; PHQ-9 sensitivity will be set to NA.")
  ppd <- ppd %>%
    mutate(
      phq9_score  = NA_real_,
      ppd_phq9_10 = NA_integer_
    )
  has_phq9 <- FALSE
}

# Postpartum months for ≤12 months restriction (if present)
pp_cand <- c(
  "Months postpartum", "Months Postpartum",
  "Months_postpartum", "Months_Postpartum",
  "Months since delivery", "Months_since_delivery"
)
postpartum_var <- intersect(pp_cand, names(ppd_raw))

if (length(postpartum_var) > 0L) {
  postpartum_var <- postpartum_var[1]
  ppd <- ppd %>%
    mutate(
      months_postpartum = suppressWarnings(as.numeric(.data[[postpartum_var]])),
      keep_leq12        = dplyr::if_else(!is.na(months_postpartum) &
                                           months_postpartum <= 12,
                                         1L, 0L)
    )
} else {
  message("Section 6: No postpartum months variable found; ≤12 months restriction will be unavailable.")
  ppd <- ppd %>%
    mutate(
      months_postpartum = NA_real_,
      keep_leq12        = NA_integer_
    )
}

## 4. Safe-named version and design matrix for modelling ------------------
outcome_var_orig <- "ppd_epds13"

name_map <- tibble(
  orig = names(ppd),
  safe = make.names(names(ppd), unique = TRUE)
)

ppd_mod <- ppd
names(ppd_mod) <- name_map$safe

outcome_var_safe <- name_map$safe[match(outcome_var_orig, name_map$orig)]
stopifnot(!is.na(outcome_var_safe))

# Exclude outcomes and obvious non-predictors from predictors
exclude_orig <- c(
  "row_id",
  "ppd_epds13", "epds_score", "EPDS Score",
  "ppd_epds11", "phq9_score", "ppd_phq9_10"
)

predictor_orig <- setdiff(names(ppd), exclude_orig)
predictor_safe <- name_map$safe[match(predictor_orig, name_map$orig)]
predictor_safe <- predictor_safe[!is.na(predictor_safe)]

# Drop predictors that are entirely NA so they don’t collapse N after drop_na
if (length(predictor_safe) > 0L) {
  na_full <- vapply(
    ppd_mod[, predictor_safe, drop = FALSE],
    FUN.VALUE = logical(1),
    FUN = function(col) all(is.na(col))
  )
  if (any(na_full)) {
    message(
      "Section 6: Dropping ", sum(na_full),
      " all-NA predictor(s) before building the design matrix."
    )
    predictor_safe <- predictor_safe[!na_full]
  }
}

data_model <- ppd_mod %>%
  dplyr::select(dplyr::all_of(c("row_id", outcome_var_safe, predictor_safe))) %>%
  tidyr::drop_na()

Y_primary <- data_model[[outcome_var_safe]]
Y_primary <- ifelse(as.numeric(Y_primary) > 0.5, 1L, 0L)

predictor_cols <- setdiff(names(data_model), c("row_id", outcome_var_safe))

X_full <- stats::model.matrix(
  ~ . - 1,
  data = data_model[, predictor_cols, drop = FALSE]
)

row_id_model <- data_model$row_id
N_full       <- nrow(X_full)
p_full       <- ncol(X_full)

message("Section 6: N = ", N_full, " observations; p = ", p_full, " predictor features.")

## 5. Metrics helpers ------------------------------------------------------
clip_probs <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

compute_auc <- function(y, p) {
  idx <- which(!is.na(y) & !is.na(p))
  y   <- y[idx]; p <- p[idx]
  if (length(unique(y)) < 2L) return(NA_real_)
  r  <- rank(p, ties.method = "average")
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  as.numeric(auc)
}

compute_brier <- function(y, p) {
  idx <- which(!is.na(y) & !is.na(p))
  y   <- y[idx]; p <- p[idx]
  mean((y - p)^2)
}

compute_calibration <- function(y, p) {
  idx <- which(!is.na(y) & !is.na(p))
  y   <- y[idx]; p <- clip_probs(p[idx])
  if (length(unique(y)) < 2L) {
    return(tibble(intercept = NA_real_, slope = NA_real_))
  }
  df <- data.frame(y = y, lp = stats::qlogis(p))
  fit <- tryCatch(
    suppressWarnings(glm(y ~ lp, data = df, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(intercept = NA_real_, slope = NA_real_))
  }
  coefs <- stats::coef(fit)
  tibble(
    intercept = unname(coefs[1]),
    slope     = ifelse(length(coefs) >= 2, unname(coefs[2]), NA_real_)
  )
}

## 6. Generic elastic-net CV fitter (for all scenarios) -------------------
fit_elastic_cv <- function(X, Y, alpha = 0.5) {
  n <- nrow(X)
  if (is.null(n)) n <- 0L
  events <- if (n > 0) sum(Y == 1) else 0L
  
  if (n == 0L || length(unique(Y)) < 2L) {
    return(list(
      metrics = tibble(
        N             = n,
        events        = events,
        event_rate    = ifelse(n > 0, events / n, NA_real_),
        auroc         = NA_real_,
        brier         = NA_real_,
        cal_intercept = NA_real_,
        cal_slope     = NA_real_
      ),
      p_cv  = rep(NA_real_, n),
      model = NULL
    ))
  }
  
  V <- dplyr::case_when(
    n >= 500 ~ 10L,
    n >= 200 ~ 5L,
    TRUE     ~ 3L
  )
  
  fold_id <- sample(rep(seq_len(V), length.out = n))
  p_cv    <- rep(NA_real_, n)
  
  for (v in seq_len(V)) {
    idx_val   <- which(fold_id == v)
    idx_train <- setdiff(seq_len(n), idx_val)
    
    Y_train <- Y[idx_train]
    X_train <- X[idx_train, , drop = FALSE]
    X_val   <- X[idx_val,   , drop = FALSE]
    
    if (length(unique(Y_train)) < 2L) {
      next
    }
    
    cvfit <- tryCatch(
      glmnet::cv.glmnet(
        x      = X_train,
        y      = Y_train,
        family = "binomial",
        alpha  = alpha,
        nfolds = min(5L, nrow(X_train))
      ),
      error = function(e) {
        message("Section 6: cv.glmnet failed in CV fold ", v,
                " (alpha = ", alpha, "): ", conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(cvfit)) {
      p_hat <- as.numeric(
        stats::predict(
          cvfit,
          newx = X_val,
          s    = "lambda.1se",
          type = "response"
        )
      )
      p_cv[idx_val] <- p_hat
    }
  }
  
  valid_idx <- which(!is.na(p_cv))
  if (length(valid_idx) < 10L || length(unique(Y[valid_idx])) < 2L) {
    metrics <- tibble(
      N             = n,
      events        = events,
      event_rate    = ifelse(n > 0, events / n, NA_real_),
      auroc         = NA_real_,
      brier         = NA_real_,
      cal_intercept = NA_real_,
      cal_slope     = NA_real_
    )
    cvfit_full <- NULL
  } else {
    Y_use <- Y[valid_idx]
    P_use <- p_cv[valid_idx]
    
    auc_val   <- compute_auc(Y_use, P_use)
    brier_val <- compute_brier(Y_use, P_use)
    cal_vals  <- compute_calibration(Y_use, P_use)
    
    metrics <- tibble(
      N             = n,
      events        = events,
      event_rate    = events / n,
      auroc         = auc_val,
      brier         = brier_val,
      cal_intercept = cal_vals$intercept,
      cal_slope     = cal_vals$slope
    )
    
    cvfit_full <- tryCatch(
      glmnet::cv.glmnet(
        x      = X,
        y      = Y,
        family = "binomial",
        alpha  = alpha,
        nfolds = min(10L, n)
      ),
      error = function(e) {
        message("Section 6: Full-sample cv.glmnet failed (alpha = ", alpha,
                "): ", conditionMessage(e))
        NULL
      }
    )
  }
  
  list(
    metrics = metrics,
    p_cv    = p_cv,
    model   = cvfit_full
  )
}

## 7. Build scenario-specific outcome vectors ------------------------------
# Mapping from row_id_model back to ppd
idx_in_ppd <- match(row_id_model, ppd$row_id)

# Primary outcome
Y_primary_scen <- Y_primary

# EPDS >= 11
Y_epds11_full <- ppd$ppd_epds11[idx_in_ppd]
Y_epds11_scen <- Y_epds11_full

# PHQ-9 >= 10 (if available)
Y_phq9_full <- ppd$ppd_phq9_10[idx_in_ppd]
if (!has_phq9 || all(is.na(Y_phq9_full))) {
  Y_phq9_scen <- numeric(0)
  X_phq9      <- X_full[integer(0), , drop = FALSE]
  idx_phq     <- integer(0)
} else {
  idx_phq     <- which(!is.na(Y_phq9_full))
  X_phq9      <- X_full[idx_phq, , drop = FALSE]
  Y_phq9_scen <- Y_phq9_full[idx_phq]
}

# ≤ 12 months postpartum (primary outcome)
keep_leq12_full <- ppd$keep_leq12[idx_in_ppd]
idx_leq12       <- which(keep_leq12_full == 1 & !is.na(keep_leq12_full))

if (length(idx_leq12) == 0L) {
  X_leq12      <- X_full[integer(0), , drop = FALSE]
  Y_leq12_scen <- numeric(0)
} else {
  X_leq12      <- X_full[idx_leq12, , drop = FALSE]
  Y_leq12_scen <- Y_primary_scen[idx_leq12]
}

## 8. Fit models for each sensitivity scenario ----------------------------

# 8.1 Primary model: EPDS >= 13, all postpartum, alpha = 0.5
res_primary <- fit_elastic_cv(
  X     = X_full,
  Y     = Y_primary_scen,
  alpha = 0.5
)

metrics_primary <- res_primary$metrics %>%
  mutate(
    scenario_id    = "primary",
    analysis_label = "Primary model",
    outcome_def    = "EPDS ≥ 13",
    restriction    = "All postpartum",
    model_type     = "Elastic-net (alpha = 0.5)"
  )

# 8.2 EPDS >= 11 sensitivity
idx_epds11 <- which(!is.na(Y_epds11_scen))
X_epds11   <- X_full[idx_epds11, , drop = FALSE]
Y_epds11   <- Y_epds11_scen[idx_epds11]

res_epds11 <- fit_elastic_cv(
  X     = X_epds11,
  Y     = Y_epds11,
  alpha = 0.5
)

metrics_epds11 <- res_epds11$metrics %>%
  mutate(
    scenario_id    = "epds11",
    analysis_label = "EPDS threshold 11",
    outcome_def    = "EPDS ≥ 11",
    restriction    = "All postpartum",
    model_type     = "Elastic-net (alpha = 0.5)"
  )

# 8.3 PHQ-9 >= 10 sensitivity (if available)
if (nrow(X_phq9) > 0L && length(unique(Y_phq9_scen)) >= 2L) {
  res_phq9 <- fit_elastic_cv(
    X     = X_phq9,
    Y     = Y_phq9_scen,
    alpha = 0.5
  )
  
  metrics_phq9 <- res_phq9$metrics %>%
    mutate(
      scenario_id    = "phq9_10",
      analysis_label = "PHQ-9 threshold 10",
      outcome_def    = "PHQ-9 ≥ 10",
      restriction    = "All postpartum (PHQ-9 observed)",
      model_type     = "Elastic-net (alpha = 0.5)"
    )
} else {
  res_phq9    <- NULL
  metrics_phq9 <- tibble(
    N             = 0L,
    events        = 0L,
    event_rate    = NA_real_,
    auroc         = NA_real_,
    brier         = NA_real_,
    cal_intercept = NA_real_,
    cal_slope     = NA_real_,
    scenario_id    = "phq9_10",
    analysis_label = "PHQ-9 threshold 10",
    outcome_def    = "PHQ-9 ≥ 10",
    restriction    = "All postpartum (PHQ-9 unavailable/missing)",
    model_type     = "Elastic-net (alpha = 0.5)"
  )
}

# 8.4 ≤12 months postpartum restriction
res_leq12 <- fit_elastic_cv(
  X     = X_leq12,
  Y     = Y_leq12_scen,
  alpha = 0.5
)

metrics_leq12 <- res_leq12$metrics %>%
  mutate(
    scenario_id    = "leq12m",
    analysis_label = "≤12 months postpartum",
    outcome_def    = "EPDS ≥ 13",
    restriction    = "≤12 months postpartum",
    model_type     = "Elastic-net (alpha = 0.5)"
  )

# 8.5 Parsimonious / alternative learner: Lasso (alpha = 1)
res_lasso <- fit_elastic_cv(
  X     = X_full,
  Y     = Y_primary_scen,
  alpha = 1.0
)

metrics_lasso <- res_lasso$metrics %>%
  mutate(
    scenario_id    = "lasso",
    analysis_label = "Parsimonious model (lasso)",
    outcome_def    = "EPDS ≥ 13",
    restriction    = "All postpartum",
    model_type     = "L1-penalised logistic (alpha = 1.0)"
  )

## 9. Policy-relevant simulations -----------------------------------------

# 9.1 Build policy data frames for each scenario --------------------------
policy_df_primary <- tibble(
  row_id = row_id_model,
  y      = Y_primary_scen,
  p_hat  = res_primary$p_cv
) %>%
  filter(!is.na(p_hat)) %>%
  left_join(ppd %>% dplyr::select(row_id, abuse_any), by = "row_id")

policy_df_epds11 <- tibble(
  row_id = row_id_model[idx_epds11],
  y      = Y_epds11,
  p_hat  = res_epds11$p_cv
) %>%
  filter(!is.na(p_hat)) %>%
  left_join(ppd %>% dplyr::select(row_id, abuse_any), by = "row_id")

if (!is.null(res_phq9) && nrow(X_phq9) > 0L) {
  policy_df_phq9 <- tibble(
    row_id = row_id_model[idx_phq],
    y      = Y_phq9_scen,
    p_hat  = res_phq9$p_cv
  ) %>%
    filter(!is.na(p_hat)) %>%
    left_join(ppd %>% dplyr::select(row_id, abuse_any), by = "row_id")
} else {
  policy_df_phq9 <- tibble()
}

policy_df_leq12 <- tibble(
  row_id = row_id_model[idx_leq12],
  y      = Y_leq12_scen,
  p_hat  = res_leq12$p_cv
) %>%
  filter(!is.na(p_hat)) %>%
  left_join(ppd %>% dplyr::select(row_id, abuse_any), by = "row_id")

policy_df_lasso <- tibble(
  row_id = row_id_model,
  y      = Y_primary_scen,
  p_hat  = res_lasso$p_cv
) %>%
  filter(!is.na(p_hat)) %>%
  left_join(ppd %>% dplyr::select(row_id, abuse_any), by = "row_id")

# 9.2 Estimate IPV effect (TMLE from Section 3 if available, else crude,
# then fall back to literature-based RR if needed) ------------------------
rr_abuse <- NA_real_

tmle_candidates <- c(
  file.path("output", "section3_prediction", "tables", "tmle_exposures_section3.csv"),
  file.path("output", "section3_prediction", "tables", "tmle_results_section3.csv")
)

tmle_files_existing <- tmle_candidates[file.exists(tmle_candidates)]

if (length(tmle_files_existing) > 0L) {
  tmle_file <- tmle_files_existing[1]
  tmle_df <- tryCatch(
    readr::read_csv(tmle_file, show_col_types = FALSE),
    error = function(e) {
      message("Section 6: Failed to read TMLE file '", tmle_file, "': ",
              conditionMessage(e))
      NULL
    }
  )
  if (!is.null(tmle_df) &&
      "exposure_var" %in% names(tmle_df) &&
      "rr_tmle"      %in% names(tmle_df) &&
      any(tmle_df$exposure_var == "abuse_any")) {
    rr_abuse_val <- tmle_df$rr_tmle[tmle_df$exposure_var == "abuse_any"][1]
    if (is.finite(rr_abuse_val) && !is.na(rr_abuse_val)) {
      rr_abuse <- rr_abuse_val
      message("Section 6: Using TMLE RR for IPV from Section 3: RR = ", round(rr_abuse, 3))
    }
  }
}

# If TMLE not available or unusable, use crude 2x2 RR
if (is.na(rr_abuse) || !is.finite(rr_abuse)) {
  message("Section 6: Using crude IPV risk ratio from current data.")
  if (nrow(policy_df_primary) > 0L &&
      any(policy_df_primary$abuse_any == 1, na.rm = TRUE) &&
      any(policy_df_primary$abuse_any == 0, na.rm = TRUE)) {
    risk1 <- mean(policy_df_primary$y[policy_df_primary$abuse_any == 1], na.rm = TRUE)
    risk0 <- mean(policy_df_primary$y[policy_df_primary$abuse_any == 0], na.rm = TRUE)
    if (!is.na(risk1) && !is.na(risk0) && risk0 > 0) {
      rr_abuse <- risk1 / risk0
      message("Section 6: Crude RR(IPV) = ", round(rr_abuse, 3))
    }
  }
}

# Literature-based RR if data-driven estimate is <= 1 or still NA
rr_abuse_lit  <- 1.80   # <<< update with pooled RR from your meta-analysis >>>
if (is.na(rr_abuse) || !is.finite(rr_abuse) || rr_abuse <= 1) {
  rr_abuse <- rr_abuse_lit
  message("Section 6: Using literature-based RR(IPV) = ", rr_abuse,
          " for policy simulation.")
}

# Evidence-based postpartum mental-health package (RR < 1)
effect_rr_pkg <- 0.70   # <<< update with meta-analytic RR for your chosen package >>>

# 9.3 Policy effect library -----------------------------------------------
policy_effects <- tibble(
  policy_id = c("ipvelim", "mh_pkg"),
  label     = c(
    sprintf("IPV elimination (RR = %.2f)", rr_abuse),
    sprintf("Mental health package (RR = %.2f)", effect_rr_pkg)
  ),
  rr        = c(rr_abuse, effect_rr_pkg)
)

# 9.4 Generic simulation function (any scenario, any policy set) ----------
simulate_policy_curves <- function(policy_df,
                                   policy_effects,
                                   target_grid = seq(0.05, 0.50, by = 0.05)) {
  if (nrow(policy_df) == 0L) return(tibble())
  
  total_cases <- sum(policy_df$y == 1, na.rm = TRUE)
  total_n     <- nrow(policy_df)
  
  res_all <- list()
  
  for (t_prop in target_grid) {
    thr        <- stats::quantile(policy_df$p_hat, probs = 1 - t_prop, na.rm = TRUE)
    idx_target <- which(policy_df$p_hat >= thr)
    n_target   <- length(idx_target)
    
    y_t     <- if (n_target > 0L) policy_df$y[idx_target] else numeric(0)
    abuse_t <- if (n_target > 0L && "abuse_any" %in% names(policy_df)) {
      policy_df$abuse_any[idx_target]
    } else {
      rep(NA_integer_, n_target)
    }
    
    cases_t    <- sum(y_t == 1, na.rm = TRUE)
    risk_obs   <- if (n_target > 0L) mean(y_t, na.rm = TRUE) else NA_real_
    frac_cases <- if (total_cases > 0) cases_t / total_cases else 0
    frac_abuse <- if (n_target > 0L) mean(abuse_t, na.rm = TRUE) else 0
    
    for (j in seq_len(nrow(policy_effects))) {
      pol_id    <- policy_effects$policy_id[j]
      pol_label <- policy_effects$label[j]
      rr        <- policy_effects$rr[j]
      
      prevented <- 0
      
      if (!is.na(risk_obs) && n_target > 0L && total_cases > 0 && risk_obs > 0) {
        if (pol_id == "ipvelim") {
          # IPV elimination: adjust mixture risk; if no IPV, effect = 0
          if (!is.na(frac_abuse) && frac_abuse > 0 && !is.na(rr) && rr > 1) {
            denom <- 1 - frac_abuse + frac_abuse * rr
            if (!is.na(denom) && denom > 0) {
              risk0 <- risk_obs / denom
              arr   <- max(0, risk_obs - risk0)
              prevented <- n_target * arr
            }
          }
        } else {
          # Generic mental-health package: multiplicative RR on risk
          if (!is.na(rr) && rr < 1) {
            risk_post <- risk_obs * rr
            arr       <- max(0, risk_obs - risk_post)
            prevented <- n_target * arr
          }
        }
      }
      
      prevented_pct <- if (total_cases > 0) prevented / total_cases else 0
      
      res_all[[length(res_all) + 1]] <- tibble(
        target_prop             = t_prop,
        n_target                = n_target,
        frac_pop_targeted       = if (total_n > 0) n_target / total_n else 0,
        cases_in_target         = cases_t,
        frac_cases_in_target    = frac_cases,
        policy_id               = pol_id,
        policy_label            = pol_label,
        prevented_cases         = prevented,
        prevented_pct_all_cases = prevented_pct
      )
    }
  }
  
  bind_rows(res_all)
}

# 9.5 Run simulations for each scenario -----------------------------------
sim_primary <- simulate_policy_curves(policy_df_primary, policy_effects) %>%
  mutate(scenario_id = "primary")

sim_epds11  <- simulate_policy_curves(policy_df_epds11, policy_effects) %>%
  mutate(scenario_id = "epds11")

sim_phq9    <- simulate_policy_curves(policy_df_phq9, policy_effects) %>%
  mutate(scenario_id = "phq9_10")

sim_leq12   <- simulate_policy_curves(policy_df_leq12, policy_effects) %>%
  mutate(scenario_id = "leq12m")

sim_lasso   <- simulate_policy_curves(policy_df_lasso, policy_effects) %>%
  mutate(scenario_id = "lasso")

sim_results <- bind_rows(
  sim_primary,
  sim_epds11,
  sim_phq9,
  sim_leq12,
  sim_lasso
) %>%
  arrange(scenario_id, policy_id, target_prop)

readr::write_csv(
  sim_results,
  file.path(out_s6_tables, "policy_simulation_curve_points_section6.csv")
)

## 10. Bootstrap CIs for top-20% targeting (mental-health package) --------
bootstrap_policy_20 <- function(policy_df,
                                policy_effects,
                                scenario_id,
                                target_prop_star = 0.20,
                                B = 500L) {
  if (nrow(policy_df) == 0L || sum(policy_df$y == 1, na.rm = TRUE) == 0L) {
    return(tibble(
      scenario_id               = scenario_id,
      prevented_cases_20_lo     = NA_real_,
      prevented_cases_20_hi     = NA_real_,
      prevented_pct_cases_20_lo = NA_real_,
      prevented_pct_cases_20_hi = NA_real_
    ))
  }
  
  mh_label <- policy_effects$label[policy_effects$policy_id == "mh_pkg"]
  
  vals <- replicate(B, {
    idx   <- sample(seq_len(nrow(policy_df)), replace = TRUE)
    df_b  <- policy_df[idx, ]
    sim_b <- simulate_policy_curves(df_b, policy_effects)
    sim_b_mh <- sim_b %>% filter(policy_label == mh_label)
    
    if (nrow(sim_b_mh) == 0L) {
      return(c(0, 0))
    }
    
    row_star <- sim_b_mh %>%
      mutate(diff20 = abs(target_prop - target_prop_star)) %>%
      arrange(diff20) %>%
      slice(1)
    
    pc <- row_star$prevented_cases
    pp <- row_star$prevented_pct_all_cases
    if (is.na(pc) || !is.finite(pc)) pc <- 0
    if (is.na(pp) || !is.finite(pp)) pp <- 0
    c(pc, pp)
  })
  
  vals <- t(vals)
  pc   <- vals[, 1]
  pp   <- vals[, 2]
  
  tibble(
    scenario_id               = scenario_id,
    prevented_cases_20_lo     = stats::quantile(pc, probs = 0.025, na.rm = TRUE),
    prevented_cases_20_hi     = stats::quantile(pc, probs = 0.975, na.rm = TRUE),
    prevented_pct_cases_20_lo = 100 * stats::quantile(pp, probs = 0.025, na.rm = TRUE),
    prevented_pct_cases_20_hi = 100 * stats::quantile(pp, probs = 0.975, na.rm = TRUE)
  )
}

boot_primary <- bootstrap_policy_20(policy_df_primary, policy_effects, "primary")
boot_epds11  <- bootstrap_policy_20(policy_df_epds11,  policy_effects, "epds11")
boot_phq9    <- bootstrap_policy_20(policy_df_phq9,    policy_effects, "phq9_10")
boot_leq12   <- bootstrap_policy_20(policy_df_leq12,   policy_effects, "leq12m")
boot_lasso   <- bootstrap_policy_20(policy_df_lasso,   policy_effects, "lasso")

boot_all <- bind_rows(
  boot_primary, boot_epds11, boot_phq9, boot_leq12, boot_lasso
)

## 11. Combine sensitivity metrics + policy summary into Table 6 ----------
metrics_all <- bind_rows(
  metrics_primary,
  metrics_epds11,
  metrics_phq9,
  metrics_leq12,
  metrics_lasso
) %>%
  # Drop scenarios with N = 0 (e.g., ≤12 months if not recorded)
  dplyr::filter(N > 0) %>%
  dplyr::select(
    scenario_id, analysis_label, outcome_def, restriction, model_type,
    N, events, event_rate, auroc, brier, cal_intercept, cal_slope
  )

# Point estimates for top ≈20% targeting under the mental-health package
mh_label <- policy_effects$label[policy_effects$policy_id == "mh_pkg"]

pol20 <- sim_results %>%
  filter(policy_label == mh_label) %>%
  mutate(diff20 = abs(target_prop - 0.20)) %>%
  group_by(scenario_id) %>%
  slice_min(order_by = diff20, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    scenario_id,
    target_prop_20         = 100 * frac_pop_targeted,
    prevented_cases_20     = prevented_cases,
    prevented_pct_cases_20 = 100 * prevented_pct_all_cases
  )

metrics_all <- metrics_all %>%
  left_join(pol20,   by = "scenario_id") %>%
  left_join(boot_all, by = "scenario_id") %>%
  mutate(
    policy_scenario_20 = dplyr::if_else(
      !is.na(target_prop_20),
      mh_label,
      NA_character_
    )
  )

# Save raw table
readr::write_csv(
  metrics_all,
  file.path(out_s6_tables, "sensitivity_and_policy_summary_section6.csv")
)

# Pretty HTML Table 6 -----------------------------------------------------
table6_gt <-
  metrics_all %>%
  mutate(
    event_rate               = round(event_rate, 3),
    auroc                    = round(auroc, 3),
    brier                    = round(brier, 3),
    cal_intercept            = round(cal_intercept, 3),
    cal_slope                = round(cal_slope, 3),
    target_prop_20           = round(target_prop_20, 1),
    prevented_cases_20       = round(prevented_cases_20, 1),
    prevented_cases_20_lo    = round(prevented_cases_20_lo, 1),
    prevented_cases_20_hi    = round(prevented_cases_20_hi, 1),
    prevented_pct_cases_20   = round(prevented_pct_cases_20, 1),
    prevented_pct_cases_20_lo= round(prevented_pct_cases_20_lo, 1),
    prevented_pct_cases_20_hi= round(prevented_pct_cases_20_hi, 1)
  ) %>%
  arrange(
    factor(
      scenario_id,
      levels = c("primary", "epds11", "phq9_10", "leq12m", "lasso")
    )
  ) %>%
  gt(groupname_col = "analysis_label") %>%
  tab_header(
    title = md("**Table 6. Sensitivity analyses and policy-relevant targeting scenarios**"),
    subtitle = md(
      "Elastic-net models under alternative outcome definitions, postpartum restrictions, and penalty structures; plus an evidence-based policy scenario targeting the top 20% highest predicted risk with a postpartum mental-health package, with bootstrap 95% CIs."
    )
  ) %>%
  cols_label(
    outcome_def              = md("**Outcome definition**"),
    restriction              = md("**Restriction**"),
    model_type               = md("**Model**"),
    N                        = md("**N**"),
    events                   = md("**Events**"),
    event_rate               = md("**Event rate**"),
    auroc                    = md("**AUROC**"),
    brier                    = md("**Brier**"),
    cal_intercept            = md("**Calib. intercept**"),
    cal_slope                = md("**Calib. slope**"),
    target_prop_20           = md("**Population targeted at 20% risk level (%)**"),
    prevented_cases_20       = md("**PPD cases prevented (top 20%, point)**"),
    prevented_cases_20_lo    = md("**PPD cases prevented (95% CI, low)**"),
    prevented_cases_20_hi    = md("**PPD cases prevented (95% CI, high)**"),
    prevented_pct_cases_20   = md("**PPD cases prevented (%, point)**"),
    prevented_pct_cases_20_lo= md("**PPD cases prevented (%, 95% CI low)**"),
    prevented_pct_cases_20_hi= md("**PPD cases prevented (%, 95% CI high)**"),
    policy_scenario_20       = md("**Policy scenario for 20% targeting**")
  ) %>%
  tab_options(
    table.font.size  = px(11),
    data_row.padding = px(2),
    heading.align    = "center"
  )

gtsave(
  table6_gt,
  filename = file.path(out_s6_tables, "table6_sensitivity_and_policy_summary.html")
)

## 12. Figure 6 – Risk-targeting vs % of PPD cases prevented ---------------
# Plot both policies for the primary scenario
fig6_df <- sim_results %>%
  filter(scenario_id == "primary") %>%
  mutate(
    target_pct          = 100 * frac_pop_targeted,
    prevented_pct_cases = 100 * prevented_pct_all_cases
  )

if (nrow(fig6_df) > 0L) {
  p_fig6 <-
    ggplot(
      fig6_df,
      aes(
        x     = target_pct,
        y     = prevented_pct_cases,
        color = policy_label,
        group = policy_label
      )
    ) +
    geom_line(size = 1) +
    geom_point(size = 1.8) +
    geom_vline(xintercept = 20, linetype = "dashed") +
    labs(
      x = "Proportion of population targeted (top % by predicted risk)",
      y = "PPD cases prevented (% of all cases)",
      title = "Figure 6. Policy simulations: targeting high predicted-risk women\nwith IPV reduction and an evidence-based postpartum mental-health package"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.title     = element_blank()
    )
  
  p_fig6
  
  ggplot2::ggsave(
    filename = file.path(out_s6_figs, "figure6_risk_targeting_policy_simulation.png"),
    plot     = p_fig6,
    width    = 7.5,
    height   = 5.5,
    dpi      = 300
  )
} else {
  message("Section 6: Simulation produced no prevention estimates; Figure 6 not drawn.")
}

############################################################
# Outputs for Section 6:
# Tables:
# - output/section6_sensitivity_policy/tables/sensitivity_and_policy_summary_section6.csv
# - output/section6_sensitivity_policy/tables/table6_sensitivity_and_policy_summary.html
# - output/section6_sensitivity_policy/tables/policy_simulation_curve_points_section6.csv
#
# Figures:
# - output/section6_sensitivity_policy/figures/figure6_risk_targeting_policy_simulation.png
############################################################
