############################################################
# Section 2 – Burden & Social Patterning of Postpartum Depression
# Dataset: PPD_dataset_v2.csv  (Bangladesh postpartum depression)
############################################################

## 1. Load packages ----
# Run once if needed:
# install.packages(c("tidyverse", "gtsummary", "gt"))

library(tidyverse)
library(gtsummary)
library(gt)

theme_gtsummary_journal(journal = "lancet")
theme_gtsummary_compact()

## 2. Set up output folders (Section 2) ----
out_root     <- "output"
out_section2 <- file.path(out_root, "section2_causal")
out_s2_tables <- file.path(out_section2, "tables")
out_s2_figs   <- file.path(out_section2, "figures")

dir.create(out_root,      showWarnings = FALSE)
dir.create(out_section2,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_s2_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_s2_figs,   recursive = TRUE, showWarnings = FALSE)

## 3. Read data & derive outcomes ----
data_path <- "PPD_dataset_v2.csv"

ppd_raw <- readr::read_csv(
  data_path,
  show_col_types = FALSE
)

ppd <- ppd_raw %>%
  mutate(
    epds_score  = as.numeric(`EPDS Score`),
    phq9_score  = as.numeric(`PHQ9 Score`),
    
    # Primary outcome: EPDS ≥13 (high PPD risk)
    ppd_epds13 = case_when(
      epds_score >= 13 ~ 1L,
      epds_score < 13  ~ 0L,
      TRUE             ~ NA_integer_
    ),
    
    # Secondary EPDS threshold: ≥11
    ppd_epds11 = case_when(
      epds_score >= 11 ~ 1L,
      epds_score < 11  ~ 0L,
      TRUE             ~ NA_integer_
    ),
    
    # PHQ-9 moderate+ depression: ≥10
    ppd_phq9_10 = case_when(
      phq9_score >= 10 ~ 1L,
      phq9_score < 10  ~ 0L,
      TRUE             ~ NA_integer_
    )
  )

N_total <- nrow(ppd)

## 4. Overall burden of postpartum depression ----
overall_prev <- ppd %>%
  summarise(
    N                   = N_total,
    prev_epds13_pct     = round(100 * mean(ppd_epds13 == 1, na.rm = TRUE), 1),
    prev_epds11_pct     = round(100 * mean(ppd_epds11 == 1, na.rm = TRUE), 1),
    prev_phq9_10_pct    = round(100 * mean(ppd_phq9_10 == 1, na.rm = TRUE), 1)
  )

overall_prev

readr::write_csv(
  overall_prev,
  file.path(out_s2_tables, "overall_ppd_prevalence_section2.csv")
)

## 5. Social & clinical patterning – helper for stratified prevalence ----

# Key domains for gradients:
# - Residence
# - Education Level
# - Current monthly income
# - Husband’s monthly income
# - Recieved Support
# - Need for Support
# - Abuse (violence)
# - Diseases during pregnancy (antenatal complications)
# - Birth compliancy (birth complications)
# - Newborn illness (neonatal complications)

group_vars <- c(
  "Residence",
  "Education Level",
  "Current monthly income",
  "Husband’s monthly income",
  "Recieved Support",
  "Need for Support",
  "Abuse",
  "Diseases during pregnancy",
  "Birth compliancy",
  "Newborn illness"
)

summarise_ppd_prevalence <- function(data, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),
      pct_of_sample   = 100 * n / N_total,
      prev_epds13_pct = 100 * mean(ppd_epds13 == 1, na.rm = TRUE),
      prev_epds11_pct = 100 * mean(ppd_epds11 == 1, na.rm = TRUE),
      prev_phq9_10_pct = 100 * mean(ppd_phq9_10 == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      domain = group_var,
      level  = as.character(.data[[group_var]])
    ) %>%
    select(
      domain, level, n,
      pct_of_sample,
      prev_epds13_pct,
      prev_epds11_pct,
      prev_phq9_10_pct
    )
}

prevalence_by_group <- purrr::map_dfr(
  group_vars,
  ~ summarise_ppd_prevalence(ppd, .x)
)

# Round nicely for display
prevalence_by_group <- prevalence_by_group %>%
  mutate(
    pct_of_sample    = round(pct_of_sample, 1),
    prev_epds13_pct  = round(prev_epds13_pct, 1),
    prev_epds11_pct  = round(prev_epds11_pct, 1),
    prev_phq9_10_pct = round(prev_phq9_10_pct, 1)
  )

# Save raw stratified prevalence for reproducibility
readr::write_csv(
  prevalence_by_group,
  file.path(out_s2_tables, "ppd_prevalence_by_group_section2.csv")
)

## 6. Table 2 – Stratified prevalence (HTML) ----

# Order domains in a clinically intuitive way
domain_order <- c(
  "Residence",
  "Education Level",
  "Current monthly income",
  "Husband’s monthly income",
  "Recieved Support",
  "Need for Support",
  "Abuse",
  "Diseases during pregnancy",
  "Birth compliancy",
  "Newborn illness"
)

table2_data <- prevalence_by_group %>%
  mutate(
    domain = factor(domain, levels = domain_order)
  ) %>%
  arrange(domain, desc(prev_epds13_pct))

table2_gt <-
  table2_data %>%
  gt(rowname_col = "level", groupname_col = "domain") %>%
  tab_header(
    title = md("**Table 2. Prevalence of postpartum depression by social and clinical characteristics**"),
    subtitle = md("Percentages are shown for EPDS ≥13, EPDS ≥11, and PHQ-9 ≥10.")
  ) %>%
  cols_label(
    n                = md("**N**"),
    pct_of_sample    = md("**% of sample**"),
    prev_epds13_pct  = md("**PPD (EPDS ≥13), %**"),
    prev_epds11_pct  = md("**PPD (EPDS ≥11), %**"),
    prev_phq9_10_pct = md("**Depression (PHQ-9 ≥10), %**")
  ) %>%
  fmt_number(
    columns = c(pct_of_sample, prev_epds13_pct, prev_epds11_pct, prev_phq9_10_pct),
    decimals = 1
  ) %>%
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(2)
  )

gtsave(
  table2_gt,
  filename = file.path(out_s2_tables, "table2_ppd_prevalence_by_social_patterning.html")
)

## 7. Figure – PPD prevalence gradients across key domains (PNG) ----

# Choose a subset of domains for the main gradient figure
domains_for_figure <- c(
  "Residence",
  "Education Level",
  "Recieved Support",
  "Abuse",
  "Birth compliancy",
  "Newborn illness"
)

figure_data <- prevalence_by_group %>%
  filter(domain %in% domains_for_figure) %>%
  mutate(
    domain = factor(domain, levels = domains_for_figure)
  ) %>%
  group_by(domain) %>%
  # Order levels within each domain by PPD prevalence (EPDS ≥13)
  mutate(
    level = forcats::fct_reorder(level, prev_epds13_pct, .desc = TRUE)
  ) %>%
  ungroup()

p_ppd_gradients <-
  ggplot(figure_data,
         aes(x = level, y = prev_epds13_pct)) +
  geom_col(alpha = 0.85) +
  facet_wrap(~ domain, scales = "free_y") +
  coord_flip() +
  labs(
    x = NULL,
    y = "PPD prevalence (EPDS ≥13), %",
    title = "Social and clinical gradients in postpartum depression (EPDS ≥13)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.text.y      = element_text(size = 9),
    panel.spacing.x  = unit(12, "pt"),
    panel.spacing.y  = unit(8, "pt")
  )

# Print to viewer
p_ppd_gradients

# Save figure as PNG
ggplot2::ggsave(
  filename = file.path(out_s2_figs, "figure_ppd_prevalence_gradients.png"),
  plot     = p_ppd_gradients,
  width    = 9,
  height   = 6,
  dpi      = 300
)

############################################################
# Outputs for Section 2:
# - output/section2_causal/tables/overall_ppd_prevalence_section2.csv
# - output/section2_causal/tables/ppd_prevalence_by_group_section2.csv
# - output/section2_causal/tables/table2_ppd_prevalence_by_social_patterning.html
# - output/section2_causal/figures/figure_ppd_prevalence_gradients.png
############################################################
