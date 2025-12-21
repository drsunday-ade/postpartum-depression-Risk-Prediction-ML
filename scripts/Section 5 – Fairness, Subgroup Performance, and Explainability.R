############################################################
# Section 5 – Fairness, Subgroup Performance, and Explainability
# Dataset: PPD_dataset_v2.csv  (Bangladesh postpartum depression)
############################################################

## 1. Packages (auto-install for reproducibility) ----
required_pkgs <- c("tidyverse", "glmnet", "gt", "iml")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)  # dplyr, ggplot2, purrr, tibble, forcats, etc.
library(glmnet)
library(gt)
library(iml)

set.seed(20251129)

## 2. Set up output folders (Section 5) ----
out_root       <- "output"
out_section5   <- file.path(out_root, "section5_fairness")
out_s5_tables  <- file.path(out_section5, "tables")
out_s5_figs    <- file.path(out_section5, "figures")

dir.create(out_root,      showWarnings = FALSE)
dir.create(out_section5,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s5_tables, recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s5_figs,   recursive    = TRUE, showWarnings = FALSE)

## 3. Read data & derive key variables (original names) ----
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
    ppd_epds13 = dplyr::if_else(epds_score >= 13, 1L, 0L),
    abuse_any  = dplyr::if_else(Abuse == "Yes", 1L, 0L)
  )

outcome_var_orig <- "ppd_epds13"

## 4. Define subgroup variables (wealth, residence, education, IPV) ----

# Helper: robust wealth tertiles from current monthly income
make_wealth_group <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  if (all(is.na(x_num))) {
    return(factor(
      rep(NA_character_, length(x_num)),
      levels = c("Low", "Middle", "High")
    ))
  }
  q <- stats::quantile(x_num, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (length(unique(q)) < 4L) {
    # fallback to binary high/low using median
    q2 <- stats::quantile(x_num, probs = c(0, 0.5, 1), na.rm = TRUE)
    if (length(unique(q2)) < 3L) {
      med <- stats::median(x_num, na.rm = TRUE)
      grp <- ifelse(x_num <= med, "Lower", "Higher")
      return(factor(grp, levels = c("Lower", "Higher")))
    } else {
      grp <- cut(
        x_num,
        breaks         = unique(q2),
        include.lowest = TRUE,
        labels         = c("Lower", "Higher")
      )
      return(grp)
    }
  } else {
    grp <- cut(
      x_num,
      breaks         = unique(q),
      include.lowest = TRUE,
      labels         = c("Low", "Middle", "High")
    )
    return(grp)
  }
}

ppd <- ppd %>%
  mutate(
    wealth_group    = make_wealth_group(`Current monthly income`),
    residence_group = factor(Residence),
    education_group = factor(`Education Level`),
    ipv_group       = factor(Abuse)  # Yes / None
  )

## 5. Create safe-named version for modelling ----
name_map <- tibble(
  orig = names(ppd),
  safe = make.names(names(ppd), unique = TRUE)
)

ppd_mod <- ppd
names(ppd_mod) <- name_map$safe

outcome_var_safe <- name_map$safe[match(outcome_var_orig, name_map$orig)]
stopifnot(!is.na(outcome_var_safe))

# Exclude outcome and EPDS total score from predictors
exclude_orig    <- c("ppd_epds13", "epds_score", "EPDS Score")
predictor_orig  <- setdiff(names(ppd), exclude_orig)
predictor_safe  <- name_map$safe[match(predictor_orig, name_map$orig)]
predictor_safe  <- predictor_safe[!is.na(predictor_safe)]

## 6. Build modelling data (complete cases, design matrix) ----
data_model <- ppd_mod %>%
  dplyr::select(dplyr::all_of(c("row_id", outcome_var_safe, predictor_safe))) %>%
  tidyr::drop_na()

Y_raw <- data_model[[outcome_var_safe]]
Y     <- ifelse(as.numeric(Y_raw) > 0.5, 1L, 0L)

predictor_cols <- setdiff(names(data_model), c("row_id", outcome_var_safe))

X <- stats::model.matrix(
  ~ . - 1,
  data = data_model[, predictor_cols, drop = FALSE]
)

row_id_model <- data_model$row_id
N <- nrow(X)
p <- ncol(X)

message("Section 5: N = ", N, " observations; p = ", p, " predictor features.")

## 7. Metrics helpers ----
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

## 8. Outer cross-validation for fairness predictions (Elastic-net) ----
V <- dplyr::case_when(
  N >= 500 ~ 10L,
  N >= 200 ~ 5L,
  TRUE     ~ 3L
)

fold_id   <- sample(rep(seq_len(V), length.out = N))
pred_fair <- rep(NA_real_, N)

for (v in seq_len(V)) {
  idx_val   <- which(fold_id == v)
  idx_train <- setdiff(seq_len(N), idx_val)
  
  Y_train <- Y[idx_train]
  Y_val   <- Y[idx_val]
  
  X_train <- X[idx_train, , drop = FALSE]
  X_val   <- X[idx_val, , drop = FALSE]
  
  if (length(unique(Y_train)) < 2L || length(unique(Y_val)) < 2L) {
    message("Section 5: Fold ", v, " has only one outcome level; skipping this fold.")
    next
  }
  
  cvfit <- tryCatch(
    glmnet::cv.glmnet(
      x      = X_train,
      y      = Y_train,
      family = "binomial",
      alpha  = 0.5,
      nfolds = min(5L, length(unique(fold_id)))
    ),
    error = function(e) {
      message("Elastic-net failed in fold ", v, ": ", conditionMessage(e))
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
    pred_fair[idx_val] <- p_hat
  }
}

message("Section 5: ", sum(!is.na(pred_fair)), " observations with CV predictions.")

## 9. Fairness dataset (attach subgroup variables) ----
fairness_df <- tibble(
  row_id = row_id_model,
  y      = Y,
  p_hat  = pred_fair
) %>%
  filter(!is.na(p_hat)) %>%
  left_join(
    ppd %>%
      dplyr::select(
        row_id,
        wealth_group,
        residence_group,
        education_group,
        ipv_group
      ),
    by = "row_id"
  )

## 10. Subgroup metrics and equal-opportunity gaps ----
compute_subgroup_metrics <- function(df, group_var, grouping_label,
                                     threshold = 0.5, min_n = 30L) {
  if (!group_var %in% names(df)) {
    return(list(
      metrics = tibble(),
      gap     = tibble(
        grouping                  = grouping_label,
        equal_opportunity_gap_abs = NA_real_
      )
    ))
  }
  
  df_group <- df %>%
    filter(!is.na(.data[[group_var]]))
  
  if (nrow(df_group) == 0L) {
    return(list(
      metrics = tibble(),
      gap     = tibble(
        grouping                  = grouping_label,
        equal_opportunity_gap_abs = NA_real_
      )
    ))
  }
  
  metrics <- df_group %>%
    group_by(subgroup = .data[[group_var]]) %>%
    summarise(
      n          = n(),
      n_events   = sum(y == 1),
      n_none     = n - n_events,
      prevalence = mean(y),
      auc        = compute_auc(y, p_hat),
      brier      = compute_brier(y, p_hat),
      sensitivity = {
        if (sum(y == 1) == 0) {
          NA_real_
        } else {
          pred <- as.integer(p_hat >= threshold)
          sum(pred == 1 & y == 1) / sum(y == 1)
        }
      },
      cal = list(compute_calibration(y, p_hat)),
      .groups = "drop"
    ) %>%
    mutate(
      intercept = purrr::map_dbl(cal, ~ .x$intercept),
      slope     = purrr::map_dbl(cal, ~ .x$slope)
    ) %>%
    dplyr::select(-cal) %>%
    filter(
      n >= min_n,
      n_events > 0,
      n_none   > 0
    ) %>%
    mutate(grouping = grouping_label) %>%
    dplyr::select(
      grouping, subgroup, n, n_events, prevalence,
      auc, brier, intercept, slope, sensitivity
    )
  
  if (nrow(metrics) >= 2L && any(!is.na(metrics$sensitivity))) {
    gap_val <- diff(range(metrics$sensitivity, na.rm = TRUE))
  } else {
    gap_val <- NA_real_
  }
  
  gap <- tibble(
    grouping                  = grouping_label,
    equal_opportunity_gap_abs = gap_val
  )
  
  list(metrics = metrics, gap = gap)
}

sub_wealth <- compute_subgroup_metrics(
  fairness_df, "wealth_group",
  "Wealth group (current monthly income)"
)
sub_res   <- compute_subgroup_metrics(
  fairness_df, "residence_group",
  "Place of residence"
)
sub_edu   <- compute_subgroup_metrics(
  fairness_df, "education_group",
  "Maternal education level"
)
sub_ipv   <- compute_subgroup_metrics(
  fairness_df, "ipv_group",
  "Intimate partner violence (IPV)"
)

subgroup_metrics_all <- bind_rows(
  sub_wealth$metrics,
  sub_res$metrics,
  sub_edu$metrics,
  sub_ipv$metrics
)

fairness_gaps_all <- bind_rows(
  sub_wealth$gap,
  sub_res$gap,
  sub_edu$gap,
  sub_ipv$gap
)

readr::write_csv(
  subgroup_metrics_all,
  file.path(out_s5_tables, "subgroup_performance_section5.csv")
)

readr::write_csv(
  fairness_gaps_all,
  file.path(out_s5_tables, "fairness_gaps_section5.csv")
)

## 11. Table 5 – Subgroup performances (HTML) ----
table5_gt <-
  subgroup_metrics_all %>%
  mutate(
    prevalence  = round(prevalence, 3),
    auc         = round(auc, 3),
    brier       = round(brier, 3),
    intercept   = round(intercept, 3),
    slope       = round(slope, 3),
    sensitivity = round(sensitivity, 3)
  ) %>%
  arrange(grouping, desc(n)) %>%
  gt(groupname_col = "grouping") %>%
  tab_header(
    title = md("**Table 5. Subgroup discrimination, calibration, and sensitivity of the postpartum depression model**"),
    subtitle = md("Elastic-net logistic model; performance reported by wealth, residence, education, and intimate partner violence (IPV) subgroups.")
  ) %>%
  cols_label(
    subgroup   = md("**Subgroup**"),
    n          = md("**N**"),
    n_events   = md("**Events**"),
    prevalence = md("**Prevalence**"),
    auc        = md("**AUROC**"),
    brier      = md("**Brier**"),
    intercept  = md("**Calib. intercept**"),
    slope      = md("**Calib. slope**"),
    sensitivity = md("**Sensitivity (threshold = 0.5)**")
  ) %>%
  tab_options(
    table.font.size  = px(11),
    data_row.padding = px(2),
    heading.align    = "center"
  )

gtsave(
  table5_gt,
  filename = file.path(out_s5_tables, "table5_subgroup_performance.html")
)

## 12. Figure 5 – Subgroup fairness (sensitivity by subgroup) ----
overall_df <- fairness_df %>%
  filter(!is.na(p_hat))

overall_sens <- {
  if (sum(overall_df$y == 1) == 0) {
    NA_real_
  } else {
    pred <- as.integer(overall_df$p_hat >= 0.5)
    sum(pred == 1 & overall_df$y == 1) / sum(overall_df$y == 1)
  }
}

plot_fairness_df <- subgroup_metrics_all %>%
  filter(!is.na(sensitivity))

p_fig5 <-
  ggplot(plot_fairness_df,
         aes(x = subgroup, y = sensitivity, fill = grouping)) +
  geom_col(alpha = 0.85, width = 0.7) +
  facet_wrap(~ grouping, scales = "free_x") +
  { if (!is.na(overall_sens))
    geom_hline(yintercept = overall_sens, linetype = "dashed") } +
  labs(
    x = NULL,
    y = "Sensitivity (TPR) at threshold = 0.5",
    title = "Figure 5. Subgroup sensitivity of the postpartum depression prediction model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  )

p_fig5

ggplot2::ggsave(
  filename = file.path(out_s5_figs, "figure5_subgroup_fairness_sensitivity.png"),
  plot     = p_fig5,
  width    = 9,
  height   = 5.5,
  dpi      = 300
)

## 13. Full-sample logistic model (for SHAP explanations) ----
# Use the numeric design matrix X already created above to avoid factor issues

X_glm <- as.data.frame(X)
Y_glm <- Y

# Drop any predictor columns with zero variance (all values identical)
nzv <- vapply(
  X_glm,
  FUN.VALUE = logical(1),
  FUN = function(col) length(unique(col)) > 1L
)

if (!all(nzv)) {
  message("Section 5: Dropping ", sum(!nzv),
          " zero-variance feature(s) before glm SHAP model.")
}

X_glm <- X_glm[, nzv, drop = FALSE]

# *** CRITICAL STEP: make all predictor names syntactically safe (no spaces etc.) ***
colnames(X_glm) <- make.names(colnames(X_glm), unique = TRUE)

# Build glm data and fit y ~ . WITHOUT reformulate() to avoid name parsing issues
glm_data <- data.frame(
  y = Y_glm,
  X_glm
)

glm_fit_full <- tryCatch(
  glm(
    y ~ .,
    data   = glm_data,
    family = binomial()
  ),
  error = function(e) {
    message("Section 5: Full-sample glm failed even after cleaning: ",
            conditionMessage(e))
    stop(e)
  }
)

## 14. SHAP global explanations (iml + Shapley) ----
message("Section 5: Computing SHAP values via 'iml' for the cleaned logistic model.")

predictor_iml <- iml::Predictor$new(
  model            = glm_fit_full,
  data             = X_glm,
  y                = Y_glm,
  predict.function = function(m, newdata) {
    as.numeric(
      stats::predict(m, newdata = newdata, type = "response")
    )
  }
)

# Subsample for SHAP to keep computation reasonable but stable
idx_sample <- sample(seq_len(nrow(X_glm)), size = min(100L, nrow(X_glm)))

shap_list <- lapply(idx_sample, function(i) {
  x_i <- X_glm[i, , drop = FALSE]
  sh  <- tryCatch(
    iml::Shapley$new(
      predictor_iml,
      x.interest  = x_i,
      sample.size = 100
    ),
    error = function(e) {
      message("Section 5: Shapley failed for obs ", i, ": ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(sh)) return(NULL)
  res <- sh$results
  res$obs_id <- i
  res
})

shap_list <- purrr::compact(shap_list)

if (length(shap_list) > 0L) {
  shap_all <- dplyr::bind_rows(shap_list)
  
  if (all(c("feature", "phi") %in% names(shap_all))) {
    global_importance <- shap_all %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(
        mean_abs_phi = mean(abs(phi), na.rm = TRUE),
        .groups      = "drop"
      )
  } else {
    message("Section 5: SHAP results missing 'feature' or 'phi'; ",
            "falling back to coefficient-based importance.")
    coef_vec <- stats::coef(glm_fit_full)
    coef_df  <- tibble::tibble(
      feature = names(coef_vec)[-1],
      coef    = as.numeric(coef_vec[-1])
    )
    global_importance <- coef_df %>%
      dplyr::transmute(
        feature,
        mean_abs_phi = abs(coef)
      )
  }
} else {
  message("Section 5: No valid SHAP results; using coefficient-based importance instead.")
  coef_vec <- stats::coef(glm_fit_full)
  coef_df  <- tibble::tibble(
    feature = names(coef_vec)[-1],
    coef    = as.numeric(coef_vec[-1])
  )
  global_importance <- coef_df %>%
    dplyr::transmute(
      feature,
      mean_abs_phi = abs(coef)
    )
}

readr::write_csv(
  global_importance,
  file.path(out_s5_tables, "global_importance_shap_section5.csv")
)

## 15. Figure 4 – SHAP global summary plot ----
global_importance_top <- global_importance %>%
  dplyr::arrange(dplyr::desc(mean_abs_phi)) %>%
  dplyr::slice_head(n = 20L) %>%
  dplyr::mutate(
    feature = forcats::fct_reorder(feature, mean_abs_phi)
  )

p_fig4 <-
  ggplot(global_importance_top,
         aes(x = mean_abs_phi, y = feature)) +
  geom_col(alpha = 0.9) +
  labs(
    x = "Mean |SHAP| contribution",
    y = NULL,
    title = "Figure 4. Global SHAP feature importance for postpartum depression risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

p_fig4

ggplot2::ggsave(
  filename = file.path(out_s5_figs, "figure4_shap_global_importance.png"),
  plot     = p_fig4,
  width    = 7.5,
  height   = 6,
  dpi      = 300
)

############################################################
# Final outputs for Section 5:
# Tables:
# - output/section5_fairness/tables/table5_subgroup_performance.html
# - output/section5_fairness/tables/subgroup_performance_section5.csv
# - output/section5_fairness/tables/fairness_gaps_section5.csv
# - output/section5_fairness/tables/global_importance_shap_section5.csv
#
# Figures:
# - output/section5_fairness/figures/figure4_shap_global_importance.png
# - output/section5_fairness/figures/figure5_subgroup_fairness_sensitivity.png
############################################################
