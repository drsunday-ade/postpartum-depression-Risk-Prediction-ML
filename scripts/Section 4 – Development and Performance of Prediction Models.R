############################################################
# Section 4 – Development and Performance of Prediction Models
# Dataset: PPD_dataset_v2.csv  (Bangladesh postpartum depression)
############################################################

## 1. Load packages ----
# Run once if needed:
# install.packages(c("tidyverse", "glmnet", "SuperLearner", "gt"))

library(tidyverse)
library(glmnet)
library(SuperLearner)
library(gt)

set.seed(20251129)

## 2. Set up output folders (Section 4) ----
out_root      <- "output"
out_section4  <- file.path(out_root, "section4_prediction")
out_s4_tables <- file.path(out_section4, "tables")
out_s4_figs   <- file.path(out_section4, "figures")

dir.create(out_root,      showWarnings = FALSE)
dir.create(out_section4,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s4_tables, recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s4_figs,   recursive    = TRUE, showWarnings = FALSE)

## 3. Read data & define outcome (original names) ----
data_path <- "PPD_dataset_v2.csv"

ppd_raw <- readr::read_csv(
  data_path,
  show_col_types = FALSE
)

ppd <- ppd_raw %>%
  mutate(
    across(where(is.character), as.factor),
    epds_score = as.numeric(`EPDS Score`),
    ppd_epds13 = dplyr::if_else(epds_score >= 13, 1L, 0L)
  )

outcome_var_orig <- "ppd_epds13"

## 4. Create safe-named version for modelling ----
name_map <- tibble(
  orig = names(ppd),
  safe = make.names(names(ppd), unique = TRUE)
)

ppd_mod <- ppd
names(ppd_mod) <- name_map$safe

outcome_var_safe <- name_map$safe[match(outcome_var_orig, name_map$orig)]
stopifnot(!is.na(outcome_var_safe))

## 5. Define predictors (broad, excluding outcome and EPDS total) ----
exclude_orig <- c("ppd_epds13", "epds_score", "EPDS Score")
predictor_orig <- setdiff(names(ppd), exclude_orig)

predictor_safe <- name_map$safe[match(predictor_orig, name_map$orig)]
predictor_safe <- predictor_safe[!is.na(predictor_safe)]

## 6. Build modelling data (complete cases, design matrix) ----
data_model <- ppd_mod %>%
  dplyr::select(dplyr::all_of(c(outcome_var_safe, predictor_safe))) %>%
  tidyr::drop_na()

Y_raw <- data_model[[outcome_var_safe]]
# Ensure numeric 0/1 outcome
Y <- ifelse(as.numeric(Y_raw) > 0.5, 1L, 0L)

predictor_cols <- setdiff(names(data_model), outcome_var_safe)

# Model matrix: dummy coding for all predictors, no intercept
X <- stats::model.matrix(
  ~ . - 1,
  data = data_model[, predictor_cols, drop = FALSE]
)

N <- nrow(X)
p <- ncol(X)

message("Section 4: N = ", N, " observations; p = ", p, " predictor features.")

## 7. Cross-validation folds (outer CV) ----
# Use 10-fold if N is large enough, otherwise fewer folds
V <- dplyr::case_when(
  N >= 500 ~ 10L,
  N >= 200 ~ 5L,
  TRUE     ~ 3L
)

fold_id <- sample(rep(seq_len(V), length.out = N))

# Storage for cross-validated predictions
pred_elnet <- rep(NA_real_, N)
pred_sl    <- rep(NA_real_, N)

## 8. Helper functions for metrics & curves ----

clip_probs <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

compute_auc <- function(y, p) {
  idx <- which(!is.na(y) & !is.na(p))
  y <- y[idx]; p <- p[idx]
  if (length(unique(y)) < 2L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  as.numeric(auc)
}

compute_brier <- function(y, p) {
  idx <- which(!is.na(y) & !is.na(p))
  y <- y[idx]; p <- p[idx]
  mean((y - p)^2)
}

compute_calibration <- function(y, p) {
  idx <- which(!is.na(y) & !is.na(p))
  y <- y[idx]; p <- clip_probs(p[idx])
  if (length(unique(y)) < 2L) {
    return(tibble(intercept = NA_real_, slope = NA_real_))
  }
  df <- data.frame(y = y, lp = stats::qlogis(p))
  fit <- tryCatch(
    glm(y ~ lp, data = df, family = binomial()),
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

calibration_curve <- function(y, p, n_groups = 10, min_per_group = 20) {
  idx <- which(!is.na(y) & !is.na(p))
  y <- y[idx]; p <- p[idx]
  if (length(unique(y)) < 2L) {
    return(tibble(pred_mean = NA_real_, obs_rate = NA_real_, group = NA_integer_))
  }
  
  # Quantile-based groups
  probs <- seq(0, 1, length.out = n_groups + 1)
  brks  <- stats::quantile(p, probs = probs, na.rm = TRUE, type = 8)
  brks  <- sort(unique(brks))
  
  if (length(brks) <= 2L) {
    return(tibble(pred_mean = NA_real_, obs_rate = NA_real_, group = NA_integer_))
  }
  
  g <- cut(p, breaks = brks, include.lowest = TRUE, labels = FALSE)
  df <- tibble(y = y, p = p, group = g)
  
  df %>%
    group_by(group) %>%
    summarise(
      n         = dplyr::n(),
      pred_mean = mean(p),
      obs_rate  = mean(y),
      .groups   = "drop"
    ) %>%
    filter(n >= min_per_group)
}

roc_points <- function(y, p, model_name) {
  idx <- which(!is.na(y) & !is.na(p))
  y <- y[idx]; p <- p[idx]
  if (length(unique(y)) < 2L) {
    return(tibble(model = model_name, fpr = numeric(0), tpr = numeric(0)))
  }
  ord <- order(p, decreasing = TRUE)
  y <- y[ord]; p <- p[ord]
  
  thresholds <- unique(p)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) {
    return(tibble(model = model_name, fpr = numeric(0), tpr = numeric(0)))
  }
  
  tpr <- numeric(length(thresholds) + 1L)
  fpr <- numeric(length(thresholds) + 1L)
  tpr[1] <- 0; fpr[1] <- 0
  
  for (i in seq_along(thresholds)) {
    t <- thresholds[i]
    pred <- as.integer(p >= t)
    tpr[i + 1L] <- sum(pred == 1 & y == 1) / n_pos
    fpr[i + 1L] <- sum(pred == 1 & y == 0) / n_neg
  }
  
  tibble(
    model = model_name,
    fpr   = fpr,
    tpr   = tpr
  )
}

## 9. Outer cross-validation: Elastic-net vs Super Learner ----

for (v in seq_len(V)) {
  idx_val   <- which(fold_id == v)
  idx_train <- setdiff(seq_len(N), idx_val)
  
  Y_train <- Y[idx_train]
  Y_val   <- Y[idx_val]
  
  X_train <- X[idx_train, , drop = FALSE]
  X_val   <- X[idx_val, , drop = FALSE]
  
  if (length(unique(Y_train)) < 2L || length(unique(Y_val)) < 2L) {
    message("Section 4: Fold ", v, " has only one outcome level; skipping this fold.")
    next
  }
  
  ## Elastic-net logistic (alpha = 0.5) ----
  cvfit <- tryCatch(
    glmnet::cv.glmnet(
      x       = X_train,
      y       = Y_train,
      family  = "binomial",
      alpha   = 0.5,
      nfolds  = 5
    ),
    error = function(e) {
      message("Elastic-net failed in fold ", v, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(cvfit)) {
    p_hat_elnet <- as.numeric(
      stats::predict(
        cvfit,
        newx = X_val,
        s    = "lambda.1se",
        type = "response"
      )
    )
    pred_elnet[idx_val] <- p_hat_elnet
  }
  
  ## Super Learner ----
  sl_fit <- tryCatch(
    SuperLearner::SuperLearner(
      Y          = Y_train,
      X          = as.data.frame(X_train),
      family     = binomial(),
      SL.library = c("SL.mean", "SL.glm", "SL.glmnet"),
      cvControl  = list(V = 5)
    ),
    error = function(e) {
      message(
        "SuperLearner (full library) failed in fold ", v, ": ",
        conditionMessage(e),
        " -- retrying with SL.mean + SL.glm."
      )
      NULL
    }
  )
  
  if (is.null(sl_fit)) {
    sl_fit <- tryCatch(
      SuperLearner::SuperLearner(
        Y          = Y_train,
        X          = as.data.frame(X_train),
        family     = binomial(),
        SL.library = c("SL.mean", "SL.glm"),
        cvControl  = list(V = 5)
      ),
      error = function(e) {
        message(
          "SuperLearner (reduced library) failed in fold ", v, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }
  
  if (!is.null(sl_fit)) {
    p_hat_sl <- as.numeric(
      predict(sl_fit, newdata = as.data.frame(X_val))$pred
    )
    pred_sl[idx_val] <- p_hat_sl
  }
}

## 10. Performance metrics: AUROC, Brier, calibration ----

metrics_elnet <- tibble(
  model  = "Elastic net",
  n_eval = sum(!is.na(pred_elnet)),
  auc    = compute_auc(Y, pred_elnet),
  brier  = compute_brier(Y, pred_elnet)
) %>%
  bind_cols(compute_calibration(Y, pred_elnet))

metrics_sl <- tibble(
  model  = "Super Learner",
  n_eval = sum(!is.na(pred_sl)),
  auc    = compute_auc(Y, pred_sl),
  brier  = compute_brier(Y, pred_sl)
) %>%
  bind_cols(compute_calibration(Y, pred_sl))

perf_table <- bind_rows(metrics_elnet, metrics_sl)

## 11. Table 4 – Prediction performance (HTML) ----

table4_gt <-
  perf_table %>%
  mutate(
    auc       = round(auc, 3),
    brier     = round(brier, 3),
    intercept = round(intercept, 3),
    slope     = round(slope, 3)
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Table 4. Discrimination and calibration of postpartum depression prediction models**"),
    subtitle = md("Outer cross-validation comparing elastic-net logistic regression and Super Learner ensemble for EPDS ≥13 postpartum depression.")
  ) %>%
  cols_label(
    model     = md("**Model**"),
    n_eval    = md("**N (used for evaluation)**"),
    auc       = md("**AUROC**"),
    brier     = md("**Brier score**"),
    intercept = md("**Calibration intercept**"),
    slope     = md("**Calibration slope**")
  ) %>%
  tab_options(
    table.font.size  = px(12),
    data_row.padding = px(2),
    heading.align    = "center"
  )

gtsave(
  table4_gt,
  filename = file.path(out_s4_tables, "table4_prediction_performance.html")
)

## 12. Figure 3 – ROC and calibration curves (PNG) ----

# ROC points
roc_elnet <- roc_points(Y, pred_elnet, "Elastic net")
roc_sl    <- roc_points(Y, pred_sl,    "Super Learner")

roc_plot_df <- bind_rows(roc_elnet, roc_sl) %>%
  mutate(
    curve_type = "ROC curve",
    x = fpr,
    y = tpr
  )

# Calibration curves (grouped)
calib_elnet <- calibration_curve(Y, pred_elnet, n_groups = 10, min_per_group = 10) %>%
  mutate(model = "Elastic net")

calib_sl <- calibration_curve(Y, pred_sl, n_groups = 10, min_per_group = 10) %>%
  mutate(model = "Super Learner")

calib_plot_df <- bind_rows(calib_elnet, calib_sl) %>%
  filter(!is.na(pred_mean)) %>%
  mutate(
    curve_type = "Calibration curve",
    x = pred_mean,
    y = obs_rate
  )

plot_df <- bind_rows(
  roc_plot_df %>% dplyr::select(model, curve_type, x, y),
  calib_plot_df %>% dplyr::select(model, curve_type, x, y)
)

p_fig3 <-
  ggplot(plot_df, aes(x = x, y = y, color = model)) +
  geom_line(size = 0.9, alpha = 0.9, na.rm = TRUE) +
  geom_point(size = 1.5, alpha = 0.9, na.rm = TRUE) +
  facet_wrap(~ curve_type, scales = "free") +
  geom_abline(
    data = data.frame(curve_type = "Calibration curve"),
    aes(slope = 1, intercept = 0),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Figure 3. ROC and calibration curves for postpartum depression prediction models"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

p_fig3

ggplot2::ggsave(
  filename = file.path(out_s4_figs, "figure3_roc_and_calibration_curves.png"),
  plot     = p_fig3,
  width    = 9,
  height   = 5.5,
  dpi      = 300
)

## 13. Full-sample final models (for key predictors / functional forms) ----

# Elastic-net on full data
cvfit_full <- glmnet::cv.glmnet(
  x      = X,
  y      = Y,
  family = "binomial",
  alpha  = 0.5,
  nfolds = V
)

coef_elnet <- as.matrix(glmnet::coef.glmnet(cvfit_full, s = "lambda.1se"))

coef_elnet_df <- tibble(
  term     = rownames(coef_elnet),
  estimate = as.numeric(coef_elnet[, 1])
) %>%
  arrange(desc(abs(estimate)))

readr::write_csv(
  coef_elnet_df,
  file.path(out_s4_tables, "elastic_net_full_coefficients.csv")
)

# Super Learner on full data
sl_full <- tryCatch(
  SuperLearner::SuperLearner(
    Y          = Y,
    X          = as.data.frame(X),
    family     = binomial(),
    SL.library = c("SL.mean", "SL.glm", "SL.glmnet"),
    cvControl  = list(V = 10)
  ),
  error = function(e) {
    message(
      "SuperLearner full fit failed (full library): ",
      conditionMessage(e),
      " -- retrying with SL.mean + SL.glm."
    )
    NULL
  }
)

if (is.null(sl_full)) {
  sl_full <- tryCatch(
    SuperLearner::SuperLearner(
      Y          = Y,
      X          = as.data.frame(X),
      family     = binomial(),
      SL.library = c("SL.mean", "SL.glm"),
      cvControl  = list(V = 10)
    ),
    error = function(e) {
      message(
        "SuperLearner full fit failed (reduced library): ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

if (!is.null(sl_full)) {
  weights_df <- tibble(
    learner = names(sl_full$coef),
    weight  = as.numeric(sl_full$coef)
  )
  
  readr::write_csv(
    weights_df,
    file.path(out_s4_tables, "superlearner_full_weights.csv")
  )
}

############################################################
# Outputs for Section 4:
# - output/section4_prediction/tables/table4_prediction_performance.html
# - output/section4_prediction/tables/elastic_net_full_coefficients.csv
# - output/section4_prediction/tables/superlearner_full_weights.csv
# - output/section4_prediction/figures/figure3_roc_and_calibration_curves.png
############################################################
