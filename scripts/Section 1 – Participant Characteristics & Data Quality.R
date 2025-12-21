############################################################
# Section 1 – Participant Characteristics & Data Quality
# Dataset: PPD_dataset_v2.csv  (Bangladesh postpartum depression)
############################################################

## 1. Load packages ----
# Run this install line once if needed:
# install.packages(c("tidyverse", "gtsummary", "janitor", "naniar", "gt"))

library(tidyverse)
library(gtsummary)
library(janitor)
library(naniar)
library(gt)

theme_gtsummary_journal(journal = "lancet")
theme_gtsummary_compact()

## 2. Set up output folders (for all sections) ----
# Root output directory
out_root <- "output"

# Section-specific subdirectories (we'll use section1 now; others for later)
out_section1    <- file.path(out_root, "section1_descriptives")
out_section2    <- file.path(out_root, "section2_causal")
out_section3    <- file.path(out_root, "section3_prediction")
out_section4    <- file.path(out_root, "section4_fairness")
out_section5    <- file.path(out_root, "section5_sensitivity")
out_section6    <- file.path(out_root, "section6_policy_sim")

# Within Section 1: tables and figures
out_s1_tables <- file.path(out_section1, "tables")
out_s1_figs   <- file.path(out_section1, "figures")

# Create directories (no warning if they already exist)
dir.create(out_root,      showWarnings = FALSE)
dir.create(out_section1,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_section2,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_section3,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_section4,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_section5,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_section6,  recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s1_tables, recursive    = TRUE, showWarnings = FALSE)
dir.create(out_s1_figs,   recursive    = TRUE, showWarnings = FALSE)

## 3. Read data ----
# Adjust path as needed
data_path <- "PPD_dataset_v2.csv"

ppd_raw <- readr::read_csv(
  data_path,
  show_col_types = FALSE
)

## 4. Create analysis variables (EPDS-based PPD) ----
ppd <- ppd_raw %>%
  mutate(
    # EPDS total score (numeric)
    epds_score = as.numeric(`EPDS Score`),
    
    # Binary postpartum depression indicator: EPDS ≥ 13
    ppd_epds13 = dplyr::case_when(
      epds_score >= 13 ~ "EPDS ≥13 (high PPD risk)",
      epds_score < 13  ~ "EPDS <13"
    ),
    ppd_epds13 = factor(
      ppd_epds13,
      levels = c("EPDS <13", "EPDS ≥13 (high PPD risk)")
    )
  )

## Quick sanity check on cohort & outcome prevalence
ppd %>%
  count(ppd_epds13) %>%
  mutate(pct = round(100 * n / sum(n), 1))

## 5. Data quality: missingness profile ----
missing_summary <- ppd %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "variable",
    values_to = "n_missing"
  ) %>%
  mutate(
    pct_missing = round(100 * n_missing / nrow(ppd), 2)
  ) %>%
  arrange(desc(n_missing))

# Print only variables with any missing data (for text description)
missing_summary %>%
  filter(n_missing > 0) %>%
  print(n = Inf)

# (Optional) Save missingness summary for records
readr::write_csv(
  missing_summary,
  file.path(out_s1_tables, "missingness_summary_section1.csv")
)

## 6. Define variables for Table 1 (baseline characteristics) ----
baseline_vars <- c(
  "Age",
  "Residence",
  "Education Level",
  "Marital status",
  "Occupation before latest pregnancy",
  "Monthly income before latest pregnancy",
  "Occupation After Your Latest Childbirth",
  "Current monthly income",
  "Husband's education level",
  "Husband’s monthly income",
  "Total children",
  "Disease before pregnancy",
  "History of pregnancy loss",
  "Family type",
  "Number of household members",
  "Abuse",
  "Pregnancy plan",
  "Diseases during pregnancy",
  "Mode of delivery",
  "Newborn illness",
  "Breastfeed"
)

## 7. Create Table 1 – Baseline characteristics overall & by PPD ----
# Use chi-square tests for categorical variables to avoid Fisher FEXACT workspace issues
table1 <-
  ppd %>%
  select(ppd_epds13, all_of(baseline_vars)) %>%
  tbl_summary(
    by = ppd_epds13,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_n() %>%
  add_p(
    test = list(
      all_categorical() ~ "chisq.test",
      all_continuous()  ~ "t.test"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(
    c("stat_1", "stat_2") ~
      "**Postpartum depression status (EPDS ≥13)**"
  )

# Convert to gt object and add caption
table1_gt <-
  table1 %>%
  as_gt() %>%
  tab_caption(
    "Table 1. Baseline characteristics of postpartum women in Bangladesh, overall and by high postpartum depression risk (EPDS ≥13)."
  )

# Save Table 1 as HTML
gtsave(
  table1_gt,
  filename = file.path(out_s1_tables, "table1_baseline_characteristics.html")
)

## 8. Key distributions for text / figures ----

# 8a. Distribution of EPDS scores overall and by PPD status
p_epds_hist <-
  ppd %>%
  ggplot(aes(x = epds_score, fill = ppd_epds13)) +
  geom_histogram(alpha = 0.55, position = "identity", bins = 20) +
  labs(
    x     = "EPDS total score",
    y     = "Count",
    fill  = "PPD status",
    title = "Distribution of EPDS scores by postpartum depression status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title      = element_text(hjust = 0.5, face = "bold")
  )

# Print to viewer
p_epds_hist

# Save as PNG
ggplot2::ggsave(
  filename = file.path(out_s1_figs, "figure1_epds_distribution_by_ppd.png"),
  plot     = p_epds_hist,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# 8b. Age distribution by PPD status
p_age_density <-
  ppd %>%
  ggplot(aes(x = Age, fill = ppd_epds13)) +
  geom_density(alpha = 0.45) +
  labs(
    x     = "Maternal age (years)",
    y     = "Density",
    fill  = "PPD status",
    title = "Age distribution by postpartum depression status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title      = element_text(hjust = 0.5, face = "bold")
  )

# Print to viewer
p_age_density

# Save as PNG
ggplot2::ggsave(
  filename = file.path(out_s1_figs, "figure2_age_density_by_ppd.png"),
  plot     = p_age_density,
  width    = 7,
  height   = 5,
  dpi      = 300
)

############################################################
# End of Section 1 script
# Outputs:
# - output/section1_descriptives/tables/table1_baseline_characteristics.html
# - output/section1_descriptives/tables/missingness_summary_section1.csv
# - output/section1_descriptives/figures/figure1_epds_distribution_by_ppd.png
# - output/section1_descriptives/figures/figure2_age_density_by_ppd.png
############################################################
