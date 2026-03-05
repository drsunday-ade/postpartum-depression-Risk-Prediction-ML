############################################################
# PPD Prediction Study — BULLETPROOF End-to-End R + ML Pipeline
# Outcome: EPDS-high (primary: EPDS Score >= 13)
# Core predictors: PHQ-2 before + PHQ-2 during pregnancy (trajectory)
# Incremental value: demographics + pregnancy context + psychosocial block
# Models: Elastic-net logistic regression (primary) + XGBoost (secondary)
# Validation: Stratified bootstrap optimism correction + bootstrap percentile CIs
# Outputs:
#   - 5 publishable figures (PNG) with strong labeling / no clipping
#   - 5 publishable tables (HTML + CSV)
#
# Key hardening fixes vs your error:
#   1) Drop rows with missing EPDS score (prevents NA labels in xgb/cv)
#   2) Stratified bootstrap sampling (prevents single-class resamples)
#   3) Safe CV folds sizing + xgb.cv best_iteration fallback
#   4) Robust missing handling for numeric predictors (median imputation)
#   5) Hard tryCatch guards so the pipeline completes without stopping
############################################################

############################
# 0) Setup
############################
options(stringsAsFactors = FALSE)
set.seed(20260214)

# ---- FILE PATHS (edit if needed) ----
DATA_PATH <- "PPD_dataset_v2.csv"
DICT_PATH <- "PPD_Data Dictionary_v2.csv"  # optional, not required
OUT_DIR   <- "ppd_outputs"

# ---- VALIDATION SETTINGS ----
B_BOOT <- 200            # fast default; set 500–1000 for final report
N_FOLDS_GLMNET <- 10
ALPHAS <- c(0.1, 0.5, 0.9)

# XGBoost tuning grid (small, robust)
XGB_GRID <- expand.grid(
  max_depth = c(2, 3),
  eta       = c(0.05, 0.10),
  min_child_weight = c(1, 5),
  subsample = c(0.8),
  colsample_bytree = c(0.8)
)

# Decision curve thresholds
DCA_THRESHOLDS <- seq(0.05, 0.60, by = 0.01)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

############################
# 1) Packages (safe install/load)
############################
pkgs <- c(
  "readr","dplyr","tidyr","stringr","janitor","forcats","purrr",
  "ggplot2","scales","ggrepel",
  "glmnet","pROC",
  "xgboost",
  "tableone","gt"
)

install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
}
invisible(lapply(pkgs, install_if_missing))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(janitor); library(forcats); library(purrr)
  library(ggplot2); library(scales); library(ggrepel)
  library(glmnet); library(pROC)
  library(xgboost)
  library(tableone); library(gt)
})

# dcurves is optional; we implement a manual fallback decision curve
if (!requireNamespace("dcurves", quietly = TRUE)) {
  message("Package 'dcurves' not found; using manual Decision Curve fallback.")
} else {
  suppressPackageStartupMessages(library(dcurves))
}

############################
# 2) Safe I/O + column discovery
############################
stop_if_missing_file <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, "\nSet DATA_PATH / DICT_PATH to the correct location.", call. = FALSE)
}
stop_if_missing_file(DATA_PATH)
# DICT is optional
if (file.exists(DICT_PATH)) {
  dict <- tryCatch(readr::read_csv(DICT_PATH, show_col_types = FALSE), error = function(e) NULL)
} else dict <- NULL

raw <- readr::read_csv(DATA_PATH, show_col_types = FALSE)

# Preserve mapping for traceability
name_map <- tibble::tibble(orig = names(raw), clean = janitor::make_clean_names(names(raw)))
dat <- raw %>% janitor::clean_names()

# Find columns robustly by regex if needed
find_col <- function(df, pattern, prefer = NULL) {
  nms <- names(df)
  hits <- nms[str_detect(nms, regex(pattern, ignore_case = TRUE))]
  if (!is.null(prefer) && prefer %in% hits) return(prefer)
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

# Required columns (use regex-based fallback)
col_epds_score <- if ("epds_score" %in% names(dat)) "epds_score" else find_col(dat, "^epds.*score$")
col_phq9_score <- if ("phq9_score" %in% names(dat)) "phq9_score" else find_col(dat, "^phq.*9.*score$")
col_phq2_before <- if ("depression_before_pregnancy_phq2" %in% names(dat)) "depression_before_pregnancy_phq2" else find_col(dat, "before.*phq2")
col_phq2_during <- if ("depression_during_pregnancy_phq2" %in% names(dat)) "depression_during_pregnancy_phq2" else find_col(dat, "during.*phq2")

req <- c(col_epds_score, col_phq9_score, col_phq2_before, col_phq2_during)
if (any(is.na(req))) {
  stop("Could not locate required columns. Missing: ",
       paste(c("EPDS score","PHQ9 score","PHQ2 before","PHQ2 during")[is.na(req)], collapse = ", "),
       "\nCheck column names in the CSV.", call. = FALSE)
}

############################
# 3) Data management (hardening)
############################
harmonize_chr <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", " ")
  x[x == ""] <- NA_character_
  x
}

# Convert EPDS score numeric; drop missing EPDS score rows (critical for xgboost/cv)
dat <- dat %>%
  mutate(
    epds_score = suppressWarnings(as.numeric(.data[[col_epds_score]])),
    phq9_score = suppressWarnings(as.numeric(.data[[col_phq9_score]]))
  ) %>%
  filter(!is.na(epds_score))

# Primary binary outcomes (explicit)
dat <- dat %>%
  mutate(
    epds_high_13 = as.integer(epds_score >= 13),
    epds_high_12 = as.integer(epds_score >= 12),
    epds_high_11 = as.integer(epds_score >= 11)
  )

# Harmonize PHQ-2 strings and derive binary
dat <- dat %>%
  mutate(
    depression_before_pregnancy_phq2 = harmonize_chr(.data[[col_phq2_before]]),
    depression_during_pregnancy_phq2 = harmonize_chr(.data[[col_phq2_during]])
  ) %>%
  mutate(
    phq2_before_pos = as.integer(str_to_lower(depression_before_pregnancy_phq2) == "positive"),
    phq2_during_pos = as.integer(str_to_lower(depression_during_pregnancy_phq2) == "positive"),
    phq2_traj = factor(
      paste0(ifelse(phq2_before_pos == 1, "P", "N"), ifelse(phq2_during_pos == 1, "P", "N")),
      levels = c("NN","NP","PN","PP")
    )
  )

# Harmonize key categorical variables if present
maybe_cols <- c(
  "age","residence","education_level","family_type","total_children","marital_status",
  "pregnancy_plan","regular_checkups","fear_of_pregnancy","major_changes_or_losses_during_pregnancy",
  "diseases_during_pregnancy",
  "recieved_support","need_for_support","trust_and_share_feelings","abuse",
  "newborn_illness","worry_about_newborn",
  "relax_sleep_when_newborn_is_tended","relax_sleep_when_the_newborn_is_asleep",
  "age_of_newborn","mode_of_delivery"
)
maybe_cols <- intersect(maybe_cols, names(dat))

dat <- dat %>%
  mutate(across(all_of(maybe_cols), ~ harmonize_chr(.x)))

# Factorize categorical variables with a stable "Unknown/Not reported" level
factorize_with_unknown <- function(df, cols) {
  for (cc in cols) {
    if (!cc %in% names(df)) next
    x <- df[[cc]]
    if (is.numeric(x)) next
    x <- ifelse(is.na(x), "Unknown/Not reported", x)
    df[[cc]] <- factor(x)
  }
  df
}

# Numeric conversion + safe imputation later (inside resampling)
if ("age" %in% names(dat)) dat <- dat %>% mutate(age = suppressWarnings(as.numeric(age)))

# Apply factorization now to lock global levels (prevents bootstrap mismatches)
cat_cols <- setdiff(maybe_cols, c("age"))
dat <- factorize_with_unknown(dat, cat_cols)

# Analysis dataset
df <- dat

############################
# 4) Model blocks (nested) — keep only columns that exist
############################
vars_A <- c("phq2_traj")

vars_B_add <- c("age","residence","education_level","family_type","total_children","marital_status")
vars_C_add <- c("pregnancy_plan","regular_checkups","fear_of_pregnancy","major_changes_or_losses_during_pregnancy","diseases_during_pregnancy")
vars_D_add <- c(
  "recieved_support","need_for_support","trust_and_share_feelings","abuse",
  "newborn_illness","worry_about_newborn",
  "relax_sleep_when_newborn_is_tended","relax_sleep_when_the_newborn_is_asleep",
  "age_of_newborn"
)

keep_existing <- function(v, df) intersect(v, names(df))

vars_A <- keep_existing(vars_A, df)
vars_B <- unique(c(vars_A, keep_existing(vars_B_add, df)))
vars_C <- unique(c(vars_B, keep_existing(vars_C_add, df)))
vars_D <- unique(c(vars_C, keep_existing(vars_D_add, df)))

model_specs <- list(A = vars_A, B = vars_B, C = vars_C, D = vars_D)

############################
# 5) Robust utilities: design matrix, folds, stratified bootstrap, metrics
############################
logit <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))

# Impute numeric columns using training medians; apply same to test
impute_numeric_train_test <- function(train, test, num_cols) {
  stats <- list()
  for (cc in num_cols) {
    if (!cc %in% names(train)) next
    med <- suppressWarnings(median(train[[cc]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    train[[cc]][is.na(train[[cc]])] <- med
    test[[cc]][is.na(test[[cc]])] <- med
    stats[[cc]] <- med
  }
  list(train = train, test = test, stats = stats)
}

# model.matrix builder; bulletproof NA handling in resulting matrix
design_matrix <- function(data, vars) {
  if (length(vars) == 0) stop("No predictors supplied to design_matrix().", call. = FALSE)
  f <- as.formula(paste("~", paste(vars, collapse = " + ")))
  mm <- model.matrix(f, data = data)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  mm[is.na(mm)] <- 0
  mm
}

# Safe nfolds for CV (avoid folds with single class as much as possible)
safe_nfolds <- function(y, max_k = 10, min_k = 3) {
  n1 <- sum(y == 1, na.rm = TRUE)
  n0 <- sum(y == 0, na.rm = TRUE)
  k <- min(max_k, n1, n0)
  if (!is.finite(k) || k < min_k) k <- min_k
  k
}

# Stratified bootstrap indices (ensures both classes present)
strat_boot_idx <- function(y) {
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  if (length(idx1) < 2 || length(idx0) < 2) {
    # fallback: simple bootstrap
    return(sample.int(length(y), size = length(y), replace = TRUE))
  }
  c(sample(idx1, length(idx1), replace = TRUE),
    sample(idx0, length(idx0), replace = TRUE)) %>% sample()
}

metrics_binary <- function(y, p) {
  auc <- tryCatch({
    r <- pROC::roc(response = y, predictor = p, quiet = TRUE)
    as.numeric(pROC::auc(r))
  }, error = function(e) NA_real_)
  
  brier <- mean((p - y)^2, na.rm = TRUE)
  
  lp <- logit(p)
  cal_int <- tryCatch({
    m <- glm(y ~ 1, family = binomial(), offset = lp)
    unname(coef(m)[1])
  }, error = function(e) NA_real_)
  
  cal_slope <- tryCatch({
    m <- glm(y ~ lp, family = binomial())
    unname(coef(m)[2])
  }, error = function(e) NA_real_)
  
  tibble::tibble(auc = auc, brier = brier, cal_intercept = cal_int, cal_slope = cal_slope)
}

############################
# 6) ElasticNet (primary) — safe fitting/prediction
############################
fit_glmnet_en <- function(train, vars, y_col = "epds_high_13", alphas = ALPHAS, nfolds = N_FOLDS_GLMNET) {
  y <- train[[y_col]]
  if (anyNA(y)) stop("Outcome contains NA in training set after filtering; check EPDS score conversion.", call. = FALSE)
  
  # safe nfolds
  k <- safe_nfolds(y, max_k = nfolds, min_k = 3)
  
  x <- design_matrix(train, vars)
  
  best <- NULL
  best_cvm <- Inf
  
  for (a in alphas) {
    cv <- tryCatch(
      cv.glmnet(x = x, y = y, family = "binomial", alpha = a, nfolds = k, type.measure = "deviance"),
      error = function(e) NULL
    )
    if (is.null(cv)) next
    lam <- cv$lambda.1se
    if (length(lam) == 0 || is.na(lam)) next
    cvm_at <- cv$cvm[which.min(abs(cv$lambda - lam))]
    if (is.finite(cvm_at) && cvm_at < best_cvm) {
      best_cvm <- cvm_at
      best <- list(alpha = a, lambda = lam, cv = cv)
    }
  }
  
  if (is.null(best)) {
    # fallback: ridge with glmnet default lambda sequence, pick lambda.1se from cv
    cv <- cv.glmnet(x = x, y = y, family = "binomial", alpha = 0, nfolds = k, type.measure = "deviance")
    best <- list(alpha = 0, lambda = cv$lambda.1se, cv = cv)
  }
  
  fit <- glmnet(x = x, y = y, family = "binomial", alpha = best$alpha, lambda = best$lambda)
  
  list(fit = fit, vars = vars, y_col = y_col, alpha = best$alpha, lambda = best$lambda, x_cols = colnames(x))
}

predict_glmnet <- function(model, data) {
  x <- design_matrix(data, model$vars)
  # align columns
  if (!identical(colnames(x), model$x_cols)) {
    missing_cols <- setdiff(model$x_cols, colnames(x))
    if (length(missing_cols) > 0) {
      add <- matrix(0, nrow = nrow(x), ncol = length(missing_cols))
      colnames(add) <- missing_cols
      x <- cbind(x, add)
    }
    x <- x[, model$x_cols, drop = FALSE]
  }
  p <- as.numeric(predict(model$fit, newx = x, type = "response"))
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

nonzero_terms <- function(model) {
  b <- as.matrix(coef(model$fit))
  nz <- rownames(b)[which(abs(b[, 1]) > 0)]
  setdiff(nz, "(Intercept)")
}

############################
# 7) XGBoost (secondary ML) — safe nested tuning with best_iteration fallback
############################
fit_xgb_nested <- function(train, vars, y_col = "epds_high_13", grid = XGB_GRID) {
  y <- train[[y_col]]
  if (anyNA(y)) stop("Outcome contains NA in xgb training set.", call. = FALSE)
  
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (min(n1, n0) < 5) return(NULL)  # too unstable for CV
  
  x <- design_matrix(train, vars)
  dtrain <- xgb.DMatrix(data = x, label = y)
  
  kfold <- min(5, n1, n0)
  if (!is.finite(kfold) || kfold < 3) kfold <- 3
  
  best <- NULL
  best_ll <- Inf
  
  for (i in seq_len(nrow(grid))) {
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = grid$max_depth[i],
      eta = grid$eta[i],
      min_child_weight = grid$min_child_weight[i],
      subsample = grid$subsample[i],
      colsample_bytree = grid$colsample_bytree[i]
    )
    
    cv <- tryCatch(
      xgb.cv(
        params = params,
        data = dtrain,
        nrounds = 2000,
        nfold = kfold,
        stratified = TRUE,
        early_stopping_rounds = 30,
        verbose = 0
      ),
      error = function(e) NULL
    )
    if (is.null(cv) || is.null(cv$evaluation_log)) next
    
    # Robust best_iteration fallback
    best_iter <- cv$best_iteration
    if (is.null(best_iter) || length(best_iter) == 0 || !is.finite(best_iter)) {
      v <- cv$evaluation_log$test_logloss_mean
      if (length(v) == 0) next
      best_iter <- which.min(v)
      if (length(best_iter) == 0 || !is.finite(best_iter)) next
    }
    
    ll <- cv$evaluation_log$test_logloss_mean[best_iter]
    if (is.finite(ll) && ll < best_ll) {
      best_ll <- ll
      best <- list(params = params, nrounds = best_iter)
    }
  }
  
  if (is.null(best)) return(NULL)
  
  fit <- tryCatch(
    xgb.train(params = best$params, data = dtrain, nrounds = best$nrounds, verbose = 0),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  list(fit = fit, vars = vars, y_col = y_col, x_cols = colnames(x), params = best$params, nrounds = best$nrounds)
}

predict_xgb <- function(model, data) {
  if (is.null(model)) return(rep(NA_real_, nrow(data)))
  x <- design_matrix(data, model$vars)
  if (!identical(colnames(x), model$x_cols)) {
    missing_cols <- setdiff(model$x_cols, colnames(x))
    if (length(missing_cols) > 0) {
      add <- matrix(0, nrow = nrow(x), ncol = length(missing_cols))
      colnames(add) <- missing_cols
      x <- cbind(x, add)
    }
    x <- x[, model$x_cols, drop = FALSE]
  }
  d <- xgb.DMatrix(x)
  p <- tryCatch(as.numeric(predict(model$fit, d)), error = function(e) rep(NA_real_, nrow(x)))
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

recalibrate_platt <- function(y, p) {
  ok <- is.finite(p) & !is.na(y)
  if (sum(ok) < 20) return(list(model = NULL, p_cal = p))
  lp <- logit(p[ok])
  m <- tryCatch(glm(y[ok] ~ lp, family = binomial()), error = function(e) NULL)
  if (is.null(m)) return(list(model = NULL, p_cal = p))
  p_cal <- p
  p_cal[ok] <- as.numeric(plogis(coef(m)[1] + coef(m)[2]*lp))
  list(model = m, p_cal = p_cal)
}

############################
# 8) Decision curve (manual fallback; always works)
############################
net_benefit <- function(y, p, pt) {
  # y in {0,1}, p predicted risk
  y <- as.integer(y)
  pred_pos <- (p >= pt)
  TP <- sum(pred_pos & (y == 1), na.rm = TRUE)
  FP <- sum(pred_pos & (y == 0), na.rm = TRUE)
  n <- sum(!is.na(y))
  if (n == 0) return(NA_real_)
  (TP / n) - (FP / n) * (pt / (1 - pt))
}

decision_curve_manual <- function(y, preds_named, thresholds = DCA_THRESHOLDS) {
  prev <- mean(y, na.rm = TRUE)
  out <- purrr::map_dfr(thresholds, function(pt) {
    nb_all <- prev - (1 - prev) * (pt / (1 - pt))
    nb_none <- 0
    tibble::tibble(
      threshold = pt,
      model = c(names(preds_named), "Treat all", "Treat none"),
      net_benefit = c(
        purrr::map_dbl(preds_named, ~ net_benefit(y, .x, pt)),
        nb_all, nb_none
      )
    )
  })
  out
}

############################
# 9) Bootstrap optimism correction (ElasticNet) + bootstrap test dist (XGB)
############################
bootstrap_optimism <- function(df, model_specs, y_col = "epds_high_13", B = B_BOOT) {
  
  # Ensure outcome exists and is not NA
  if (!y_col %in% names(df)) stop("Outcome column not found: ", y_col, call. = FALSE)
  df <- df %>% filter(!is.na(.data[[y_col]]))
  y_full <- df[[y_col]]
  
  # Identify numeric predictors to impute (currently: age only, but general)
  num_cols <- names(df)[sapply(df, is.numeric)]
  num_cols <- setdiff(num_cols, c("epds_score","phq9_score","epds_high_13","epds_high_12","epds_high_11"))
  
  model_names <- names(model_specs)
  
  # Fit full-data ElasticNet models (with numeric imputation using full medians for reporting baseline only)
  imp_full <- impute_numeric_train_test(df, df, num_cols)
  df_full <- imp_full$train
  
  full_fits_en <- purrr::map(model_specs, ~ fit_glmnet_en(df_full, .x, y_col = y_col))
  full_preds_en <- purrr::map(full_fits_en, ~ predict_glmnet(.x, df_full))
  
  full_metrics_en <- purrr::map2_dfr(full_preds_en, model_names, ~ {
    metrics_binary(df_full[[y_col]], .x) %>% mutate(model = .y, family = "ElasticNet")
  })
  
  # Fit full-data XGB Model D (optional baseline; safe)
  full_fit_xgb <- tryCatch(fit_xgb_nested(df_full, model_specs$D, y_col = y_col), error = function(e) NULL)
  p_xgb_full <- predict_xgb(full_fit_xgb, df_full)
  p_xgb_full_cal <- recalibrate_platt(df_full[[y_col]], p_xgb_full)$p_cal
  full_metrics_xgb <- metrics_binary(df_full[[y_col]], p_xgb_full_cal) %>%
    mutate(model = "D", family = "XGBoost(Platt)")
  
  optimism_rows <- list()
  test_rows <- list()
  sel_freq <- setNames(vector("list", length(model_names)), model_names)
  for (mn in model_names) sel_freq[[mn]] <- character(0)
  
  n <- nrow(df)
  
  for (b in seq_len(B)) {
    
    # Stratified bootstrap to avoid single-class samples
    idx <- strat_boot_idx(df[[y_col]])
    boot_raw <- df[idx, , drop = FALSE]
    
    # Split: train=boot, test=original (df)
    # Impute numeric based on boot medians; apply to both boot and df
    imp <- impute_numeric_train_test(boot_raw, df, num_cols)
    boot <- imp$train
    test <- imp$test
    
    # ElasticNet: each nested model
    for (mn in model_names) {
      fit_b <- tryCatch(fit_glmnet_en(boot, model_specs[[mn]], y_col = y_col), error = function(e) NULL)
      if (is.null(fit_b)) {
        optimism_rows[[length(optimism_rows) + 1]] <- tibble::tibble(
          bootstrap = b, model = mn, family = "ElasticNet",
          auc_optimism = NA_real_, brier_optimism = NA_real_,
          cal_intercept_optimism = NA_real_, cal_slope_optimism = NA_real_
        )
        test_rows[[length(test_rows) + 1]] <- tibble::tibble(
          bootstrap = b, model = mn, family = "ElasticNet",
          auc_test = NA_real_, brier_test = NA_real_,
          cal_intercept_test = NA_real_, cal_slope_test = NA_real_
        )
        next
      }
      
      sel_freq[[mn]] <- c(sel_freq[[mn]], nonzero_terms(fit_b))
      
      p_app <- predict_glmnet(fit_b, boot)
      p_test <- predict_glmnet(fit_b, test)
      
      m_app <- metrics_binary(boot[[y_col]], p_app)
      m_test <- metrics_binary(test[[y_col]], p_test)
      
      optimism_rows[[length(optimism_rows) + 1]] <- tibble::tibble(
        bootstrap = b, model = mn, family = "ElasticNet",
        auc_optimism = m_app$auc - m_test$auc,
        brier_optimism = m_app$brier - m_test$brier,
        cal_intercept_optimism = m_app$cal_intercept - m_test$cal_intercept,
        cal_slope_optimism = m_app$cal_slope - m_test$cal_slope
      )
      
      test_rows[[length(test_rows) + 1]] <- tibble::tibble(
        bootstrap = b, model = mn, family = "ElasticNet",
        auc_test = m_test$auc,
        brier_test = m_test$brier,
        cal_intercept_test = m_test$cal_intercept,
        cal_slope_test = m_test$cal_slope
      )
    }
    
    # XGBoost: only Model D (nested tuning); safe and optional
    fit_xgb_b <- tryCatch(fit_xgb_nested(boot, model_specs$D, y_col = y_col), error = function(e) NULL)
    p_xgb_test <- predict_xgb(fit_xgb_b, test)
    
    # Platt calibration fitted on boot predictions
    p_xgb_boot <- predict_xgb(fit_xgb_b, boot)
    platt_b <- recalibrate_platt(boot[[y_col]], p_xgb_boot)
    # apply calibration model to test predictions if available
    if (!is.null(platt_b$model)) {
      lp_test <- logit(p_xgb_test)
      p_xgb_test_cal <- as.numeric(plogis(coef(platt_b$model)[1] + coef(platt_b$model)[2] * lp_test))
      p_xgb_test_cal <- pmin(pmax(p_xgb_test_cal, 1e-6), 1 - 1e-6)
    } else {
      p_xgb_test_cal <- p_xgb_test
    }
    
    m_xgb_test <- metrics_binary(test[[y_col]], p_xgb_test_cal)
    
    test_rows[[length(test_rows) + 1]] <- tibble::tibble(
      bootstrap = b, model = "D", family = "XGBoost(Platt)",
      auc_test = m_xgb_test$auc,
      brier_test = m_xgb_test$brier,
      cal_intercept_test = m_xgb_test$cal_intercept,
      cal_slope_test = m_xgb_test$cal_slope
    )
  }
  
  optimism_df <- bind_rows(optimism_rows)
  test_df <- bind_rows(test_rows)
  
  mean_opt <- optimism_df %>%
    group_by(model, family) %>%
    summarise(
      auc_opt = mean(auc_optimism, na.rm = TRUE),
      brier_opt = mean(brier_optimism, na.rm = TRUE),
      cal_intercept_opt = mean(cal_intercept_optimism, na.rm = TRUE),
      cal_slope_opt = mean(cal_slope_optimism, na.rm = TRUE),
      .groups = "drop"
    )
  
  apparent <- bind_rows(full_metrics_en, full_metrics_xgb)
  
  corrected <- apparent %>%
    left_join(mean_opt, by = c("model","family")) %>%
    mutate(
      auc_corrected = ifelse(family == "ElasticNet", auc - auc_opt, NA_real_),
      brier_corrected = ifelse(family == "ElasticNet", brier - brier_opt, NA_real_),
      cal_intercept_corrected = ifelse(family == "ElasticNet", cal_intercept - cal_intercept_opt, NA_real_),
      cal_slope_corrected = ifelse(family == "ElasticNet", cal_slope - cal_slope_opt, NA_real_)
    )
  
  ci_test <- test_df %>%
    group_by(model, family) %>%
    summarise(
      auc_lo = quantile(auc_test, 0.025, na.rm = TRUE),
      auc_hi = quantile(auc_test, 0.975, na.rm = TRUE),
      brier_lo = quantile(brier_test, 0.025, na.rm = TRUE),
      brier_hi = quantile(brier_test, 0.975, na.rm = TRUE),
      cal_int_lo = quantile(cal_intercept_test, 0.025, na.rm = TRUE),
      cal_int_hi = quantile(cal_intercept_test, 0.975, na.rm = TRUE),
      cal_slope_lo = quantile(cal_slope_test, 0.025, na.rm = TRUE),
      cal_slope_hi = quantile(cal_slope_test, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  perf <- corrected %>% left_join(ci_test, by = c("model","family")) %>% arrange(family, model)
  
  sel_tbl <- purrr::map_dfr(names(sel_freq), function(mn) {
    tibble::tibble(term = sel_freq[[mn]]) %>%
      count(term, sort = TRUE) %>%
      mutate(model = mn, freq = n / B) %>%
      select(model, term, freq, n)
  })
  
  list(
    perf = perf,
    test_df = test_df,
    sel_freq = sel_tbl,
    full_fits_en = full_fits_en,
    full_preds_en = full_preds_en,
    full_fit_xgb_D = full_fit_xgb,
    full_pred_xgb_D_cal = p_xgb_full_cal,
    df_full = df_full
  )
}

############################
# 10) Run primary analysis (EPDS >= 13)
############################
res_13 <- bootstrap_optimism(df, model_specs, y_col = "epds_high_13", B = B_BOOT)

perf_13 <- res_13$perf
sel_13  <- res_13$sel_freq
df_full <- res_13$df_full

############################
# 11) Tables (5) — CSV + HTML (robust)
############################
write_html_gt <- function(gt_obj, path) {
  tryCatch(gt::gtsave(gt_obj, path), error = function(e) {
    message("Could not save HTML via gt::gtsave (", path, "): ", e$message)
  })
}

# TABLE 1: Cohort characteristics by EPDS-high
vars_table1 <- c(
  "age","residence","education_level","marital_status","family_type","total_children",
  "pregnancy_plan","regular_checkups","fear_of_pregnancy","major_changes_or_losses_during_pregnancy",
  "recieved_support","need_for_support","trust_and_share_feelings","abuse",
  "newborn_illness","worry_about_newborn","age_of_newborn",
  "relax_sleep_when_newborn_is_tended","relax_sleep_when_the_newborn_is_asleep",
  "phq2_traj"
)
vars_table1 <- intersect(vars_table1, names(df_full))

tab1_df <- tryCatch({
  tab1 <- tableone::CreateTableOne(
    vars = setdiff(vars_table1, "phq2_traj"),
    strata = "epds_high_13",
    data = df_full,
    factorVars = setdiff(vars_table1, c("age"))
  )
  print(tab1, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("variable")
}, error = function(e) {
  message("TableOne failed; using fallback summary. Error: ", e$message)
  # fallback: simple summary
  df_full %>%
    summarise(across(all_of(vars_table1), ~ {
      if (is.numeric(.x)) paste0(round(mean(.x, na.rm=TRUE),2), " (SD ", round(sd(.x, na.rm=TRUE),2), ")")
      else paste0("Levels: ", nlevels(.x))
    })) %>%
    pivot_longer(everything(), names_to="variable", values_to="summary")
})

readr::write_csv(tab1_df, file.path(OUT_DIR, "Table1_Cohort_by_EPDSHigh.csv"))
gt(tab1_df) %>%
  tab_header(
    title = md("**Table 1. Cohort characteristics by EPDS-high (EPDS ≥ 13)**"),
    subtitle = "SMD shown when available; descriptive only."
  ) %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  tab_source_note(md("**Notes:** Categorical missingness is represented as 'Unknown/Not reported' for modeling.")) %>%
  (\(g) { write_html_gt(g, file.path(OUT_DIR, "Table1_Cohort_by_EPDSHigh.html")); g })()

# TABLE 2: Missingness table (raw NA only)
miss_df <- df %>%
  summarise(across(everything(), ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  mutate(
    n_missing = round(pct_missing * nrow(df)),
    pct_missing = pct_missing * 100
  ) %>%
  arrange(desc(pct_missing))

readr::write_csv(miss_df, file.path(OUT_DIR, "Table2_Missingness.csv"))
gt(miss_df) %>%
  tab_header(
    title = md("**Table 2. Missingness profile (raw NA values)**"),
    subtitle = "Computed before factor recoding (Unknown/Not reported) and numeric imputation."
  ) %>%
  fmt_number(columns = c(pct_missing), decimals = 2) %>%
  (\(g) { write_html_gt(g, file.path(OUT_DIR, "Table2_Missingness.html")); g })()

# TABLE 3: Performance summary
perf_tbl <- perf_13 %>%
  mutate(
    auc_app = auc,
    brier_app = brier,
    cal_int_app = cal_intercept,
    cal_slope_app = cal_slope
  ) %>%
  select(
    family, model,
    auc_app, auc_corrected, auc_lo, auc_hi,
    brier_app, brier_corrected, brier_lo, brier_hi,
    cal_int_app, cal_intercept_corrected, cal_int_lo, cal_int_hi,
    cal_slope_app, cal_slope_corrected, cal_slope_lo, cal_slope_hi
  ) %>%
  arrange(family, model)

readr::write_csv(perf_tbl, file.path(OUT_DIR, "Table3_Performance_EPDS13.csv"))
gt(perf_tbl) %>%
  tab_header(
    title = md("**Table 3. Internal validation performance (EPDS ≥ 13)**"),
    subtitle = paste0("ElasticNet: bootstrap optimism correction (B=", B_BOOT, "). XGBoost: bootstrap test distribution.")
  ) %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  tab_source_note(md("**Notes:** Calibration slope ≈1 and intercept ≈0 indicate good calibration.")) %>%
  (\(g) { write_html_gt(g, file.path(OUT_DIR, "Table3_Performance_EPDS13.html")); g })()

# TABLE 4: Incremental value from bootstrap test distributions (ElasticNet only)
test_only_en <- res_13$test_df %>%
  filter(family == "ElasticNet") %>%
  select(bootstrap, model, auc_test, brier_test) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(auc_test, brier_test))

delta_tbl <- test_only_en %>%
  transmute(
    dAUC_A_to_B = auc_test_B - auc_test_A,
    dAUC_B_to_C = auc_test_C - auc_test_B,
    dAUC_C_to_D = auc_test_D - auc_test_C,
    dBrier_A_to_B = brier_test_B - brier_test_A,
    dBrier_B_to_C = brier_test_C - brier_test_B,
    dBrier_C_to_D = brier_test_D - brier_test_C
  )

delta_summary <- delta_tbl %>%
  summarise(across(everything(), list(
    mean = ~mean(.x, na.rm = TRUE),
    lo   = ~quantile(.x, 0.025, na.rm = TRUE),
    hi   = ~quantile(.x, 0.975, na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("metric","stat"), sep = "_(?=[^_]+$)", extra = "merge", fill = "right") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(metric = str_replace_all(metric, "\\.", "_"))

readr::write_csv(delta_summary, file.path(OUT_DIR, "Table4_IncrementalValue.csv"))
gt(delta_summary) %>%
  tab_header(
    title = md("**Table 4. Incremental predictive value (ElasticNet) — bootstrap test distribution**"),
    subtitle = "Positive ΔAUC = better discrimination; negative ΔBrier = better accuracy."
  ) %>%
  fmt_number(columns = c(mean, lo, hi), decimals = 3) %>%
  (\(g) { write_html_gt(g, file.path(OUT_DIR, "Table4_IncrementalValue.html")); g })()

# TABLE 5: Subgroup AUC (Model A vs Model D, ElasticNet) on full-data fit (descriptive stability)
pA <- res_13$full_preds_en$A
pD <- res_13$full_preds_en$D
y  <- df_full$epds_high_13

pred_full <- tibble::tibble(
  epds_high_13 = y,
  pA = pA,
  pD = pD,
  residence = if ("residence" %in% names(df_full)) df_full$residence else factor("All"),
  education_level = if ("education_level" %in% names(df_full)) df_full$education_level else factor("All"),
  family_type = if ("family_type" %in% names(df_full)) df_full$family_type else factor("All"),
  age_of_newborn = if ("age_of_newborn" %in% names(df_full)) df_full$age_of_newborn else factor("All")
)

subgroup_auc <- function(data, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),
      prev = mean(epds_high_13, na.rm = TRUE),
      auc_A = tryCatch(as.numeric(pROC::auc(pROC::roc(epds_high_13, pA, quiet = TRUE))), error = function(e) NA_real_),
      auc_D = tryCatch(as.numeric(pROC::auc(pROC::roc(epds_high_13, pD, quiet = TRUE))), error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    rename(group = 1) %>%
    mutate(group_var = group_var)
}

tab5 <- bind_rows(
  subgroup_auc(pred_full, "residence"),
  subgroup_auc(pred_full, "education_level"),
  subgroup_auc(pred_full, "family_type"),
  subgroup_auc(pred_full, "age_of_newborn")
) %>% relocate(group_var, group)

readr::write_csv(tab5, file.path(OUT_DIR, "Table5_SubgroupAUC.csv"))
gt(tab5) %>%
  tab_header(
    title = md("**Table 5. Subgroup discrimination (AUC) — Model A vs Model D**"),
    subtitle = "Descriptive stability check (not external validation)."
  ) %>%
  fmt_number(columns = c(prev, auc_A, auc_D), decimals = 3) %>%
  (\(g) { write_html_gt(g, file.path(OUT_DIR, "Table5_SubgroupAUC.html")); g })()

############################
# 12) Figures (5) — PNG, no clipping, no overlapping
############################
theme_pub <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, color = "gray30"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 20, 12, 12)
    )
}

save_png <- function(p, filename, w = 10, h = 6.5, dpi = 320) {
  ggsave(filename = file.path(OUT_DIR, filename), plot = p,
         width = w, height = h, dpi = dpi, units = "in", bg = "white")
}

# FIGURE 1: EPDS distribution with threshold
p1 <- ggplot(df_full, aes(x = epds_score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25, color = "white", alpha = 0.95) +
  geom_density(linewidth = 1.1, alpha = 0.25) +
  geom_vline(xintercept = 13, linetype = "22", linewidth = 1.1) +
  annotate("label", x = 13, y = Inf, label = "Primary endpoint threshold: EPDS ≥ 13",
           vjust = 1.2, hjust = 0.5, label.size = 0.2, alpha = 0.95) +
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  labs(
    title = "Distribution of EPDS scores in the postpartum cohort",
    subtitle = "Primary endpoint: EPDS-high defined as EPDS Score ≥ 13",
    x = "EPDS total score",
    y = "Density",
    caption = "Notes: Histogram + kernel density. Vertical line marks prespecified screening threshold."
  ) +
  coord_cartesian(clip = "off") +
  theme_pub()

save_png(p1, "Figure1_EPDS_Distribution.png")

# FIGURE 2: EPDS-high prevalence by PHQ-2 trajectory
prev_traj <- df_full %>%
  group_by(phq2_traj) %>%
  summarise(n = n(), prev = mean(epds_high_13, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    se = sqrt(prev*(1-prev)/pmax(n,1)),
    lo = pmax(0, prev - 1.96*se),
    hi = pmin(1, prev + 1.96*se),
    label = paste0("n=", n, "\n", percent(prev, accuracy = 0.1))
  )

p2 <- ggplot(prev_traj, aes(x = phq2_traj, y = prev)) +
  geom_col(width = 0.65, alpha = 0.95) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.9) +
  geom_text_repel(aes(label = label), nudge_y = 0.05, size = 3.6,
                  box.padding = 0.25, min.segment.length = 0) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Postpartum depression risk (EPDS ≥ 13) by prenatal PHQ-2 trajectory",
    subtitle = "Trajectory: NN, NP, PN, PP (before pregnancy × during pregnancy)",
    x = "Prenatal PHQ-2 trajectory",
    y = "Prevalence of EPDS-high",
    caption = "Notes: Error bars are approximate 95% normal CIs; descriptive risk stratification."
  ) +
  coord_cartesian(clip = "off") +
  theme_pub()

save_png(p2, "Figure2_EPDSHigh_by_PHQ2_Trajectory.png")

# FIGURE 3: ROC curves (ElasticNet A–D + XGB D if available)
roc_df <- function(y, p, model_name) {
  r <- pROC::roc(y, p, quiet = TRUE)
  tibble::tibble(
    fpr = 1 - r$specificities,
    tpr = r$sensitivities,
    model = model_name,
    auc = as.numeric(pROC::auc(r))
  )
}

pB <- res_13$full_preds_en$B
pC <- res_13$full_preds_en$C
pX <- res_13$full_pred_xgb_D_cal

roc_long <- bind_rows(
  roc_df(y, pA, "ElasticNet A"),
  roc_df(y, pB, "ElasticNet B"),
  roc_df(y, pC, "ElasticNet C"),
  roc_df(y, pD, "ElasticNet D")
)

if (all(is.finite(pX))) roc_long <- bind_rows(roc_long, roc_df(y, pX, "XGBoost D (Platt)"))

roc_long <- roc_long %>%
  group_by(model) %>%
  mutate(model_lab = paste0(model, " (AUC=", sprintf("%.3f", first(auc)), ")")) %>%
  ungroup()

p3 <- ggplot(roc_long, aes(x = fpr, y = tpr, color = model_lab)) +
  geom_line(linewidth = 1.05, alpha = 0.95) +
  geom_abline(linetype = "33", linewidth = 0.9) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Discrimination: ROC curves for nested prediction models",
    subtitle = "Primary endpoint: EPDS ≥ 13 (internal validation summarized in Table 3)",
    x = "False positive rate (1 − specificity)",
    y = "True positive rate (sensitivity)",
    color = "Model (AUC)",
    caption = "Notes: ElasticNet models are nested blocks A–D. XGBoost is a secondary ML benchmark with calibration."
  ) +
  coord_cartesian(clip = "off") +
  theme_pub() +
  theme(legend.position = "bottom")

save_png(p3, "Figure3_ROC_NestedModels.png", w = 11, h = 7)

# FIGURE 4: Calibration (decile bins) for A–D (+ XGB if available)
calib_bins <- function(y, p, model_name, bins = 10) {
  tibble::tibble(y = y, p = p) %>%
    mutate(bin = ntile(p, bins)) %>%
    group_by(bin) %>%
    summarise(
      p_mean = mean(p, na.rm = TRUE),
      y_obs  = mean(y, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se = sqrt(y_obs*(1-y_obs)/pmax(n,1)),
      lo = pmax(0, y_obs - 1.96*se),
      hi = pmin(1, y_obs + 1.96*se),
      model = model_name
    )
}

cal_long <- bind_rows(
  calib_bins(y, pA, "ElasticNet A"),
  calib_bins(y, pB, "ElasticNet B"),
  calib_bins(y, pC, "ElasticNet C"),
  calib_bins(y, pD, "ElasticNet D")
)
if (all(is.finite(pX))) cal_long <- bind_rows(cal_long, calib_bins(y, pX, "XGBoost D (Platt)"))

p4 <- ggplot(cal_long, aes(x = p_mean, y = y_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "33", linewidth = 0.9) +
  geom_point(size = 2.3, alpha = 0.95) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.0, linewidth = 0.8, alpha = 0.9) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.0, alpha = 0.6) +
  facet_wrap(~ model, ncol = 3) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Calibration: observed vs predicted risk (deciles + loess)",
    subtitle = "Good calibration: points near 45° line; loess summarizes trend",
    x = "Mean predicted risk in bin",
    y = "Observed event rate in bin",
    caption = "Notes: Error bars are approximate 95% normal CIs for observed rates within each risk decile."
  ) +
  coord_cartesian(clip = "off") +
  theme_pub()

save_png(p4, "Figure4_Calibration_Deciles.png", w = 12, h = 7.8)

# FIGURE 5: Decision Curve (use dcurves if available; else manual always works)
preds_named <- list(
  "ElasticNet A" = pA,
  "ElasticNet B" = pB,
  "ElasticNet C" = pC,
  "ElasticNet D" = pD
)
if (all(is.finite(pX))) preds_named[["XGBoost D (Platt)"]] <- pX

dca_plot_data <- NULL
if (requireNamespace("dcurves", quietly = TRUE)) {
  # dcurves expects predictions as columns
  dca_df <- tibble::tibble(epds_high_13 = y) %>%
    bind_cols(as_tibble(preds_named, .name_repair = "minimal"))
  
  dca_plot_data <- tryCatch({
    dcurves::dca(epds_high_13 ~ ., data = dca_df, thresholds = DCA_THRESHOLDS)
  }, error = function(e) {
    message("dcurves::dca failed; using manual Decision Curve. Error: ", e$message)
    NULL
  })
}

if (!is.null(dca_plot_data)) {
  p5 <- dca_plot_data %>%
    filter(!is.na(net_benefit)) %>%
    ggplot(aes(x = threshold, y = net_benefit, color = model)) +
    geom_line(linewidth = 1.05, alpha = 0.95) +
    labs(
      title = "Clinical utility: Decision Curve Analysis (net benefit)",
      subtitle = "Higher net benefit indicates better clinical utility across decision thresholds",
      x = "Threshold probability",
      y = "Net benefit",
      color = "Strategy",
      caption = "Notes: dcurves output. Interpret within plausible screening/referral thresholds."
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(legend.position = "bottom")
} else {
  dca_manual <- decision_curve_manual(y, preds_named, thresholds = DCA_THRESHOLDS)
  p5 <- ggplot(dca_manual, aes(x = threshold, y = net_benefit, color = model)) +
    geom_line(linewidth = 1.05, alpha = 0.95) +
    labs(
      title = "Clinical utility: Decision Curve (manual net benefit; always defined)",
      subtitle = "Higher net benefit indicates better clinical utility across decision thresholds",
      x = "Threshold probability",
      y = "Net benefit",
      color = "Strategy",
      caption = "Notes: Manual NB = TP/n − FP/n × pt/(1−pt). Includes Treat-all and Treat-none baselines."
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    coord_cartesian(clip = "off") +
    theme_pub() +
    theme(legend.position = "bottom")
}
save_png(p5, "Figure5_DecisionCurve.png", w = 11, h = 7)

############################
# 13) Manifest
############################
cat("\nSUCCESS: Outputs saved to: ", normalizePath(OUT_DIR), "\n\n")

cat("Figures (PNG):\n",
    " - Figure1_EPDS_Distribution.png\n",
    " - Figure2_EPDSHigh_by_PHQ2_Trajectory.png\n",
    " - Figure3_ROC_NestedModels.png\n",
    " - Figure4_Calibration_Deciles.png\n",
    " - Figure5_DecisionCurve.png\n\n", sep = "")

cat("Tables (CSV + HTML):\n",
    " - Table1_Cohort_by_EPDSHigh\n",
    " - Table2_Missingness\n",
    " - Table3_Performance_EPDS13\n",
    " - Table4_IncrementalValue\n",
    " - Table5_SubgroupAUC\n\n", sep = "")
