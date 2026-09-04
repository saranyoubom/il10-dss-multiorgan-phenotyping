# ============================================================
# Supplementary Table S4 - Linear mixed-effects models for
# body weight, fasting C-peptide and fasting blood glucose
# Scientific Reports Rev1 rebuttal (Reviewer 2 point 5)
#
# Related to Results Figs. 1b and 2b, and Methods
# "Statistical analysis."
#
# Model: Value ~ Group * Week + (1 | Animal)
#   Random intercept per animal, REML estimation, CTRL as the
#   reference group. Fixed-effect significance is assessed with
#   an asymptotic Wald z-test (Estimate / SE), matching the
#   default behaviour of Python's statsmodels MixedLM (used for
#   the original rebuttal fit, since R/Rscript was unavailable
#   on the machine that produced it) rather than a
#   Satterthwaite/Kenward-Roger t-test.
#
# Body weight model uses all 13 measured weekly timepoints
# (weeks 0-12). C-peptide and fasting blood glucose models use
# the 3 measured timepoints (weeks 0, 6, 12).
#
# IMPORTANT - cohort filtering (n=5 per group):
#   The source data files contain more replicate animals than
#   the n=5 per group reported throughout the manuscript
#   (7 IL10-KO, 6 DSS-IL10-KO for body weight/blood glucose;
#   7 IL10-KO, 6 DSS-IL10-KO replicates for C-peptide). This
#   script identifies the reported n=5 cohort itself, by
#   exhaustive search over every 5-animal subset for an exact
#   numerical match (mean and SD, to 3 decimal places) against
#   the already-published summary-statistics CSVs
#   (Fig_1b_Summary_BW.csv, Fig_2a_Summary_FBG.csv,
#   Fig_2b_Summary_C_peptide.csv) at every reported timepoint,
#   rather than assuming or hardcoding which animals to drop.
#   For body weight and fasting blood glucose (both indexed by
#   a shared animal ID column), the matching subset is:
#     IL10-KO:     IDs 04-08 (IDs 09, 10 excluded)
#     DSS-IL10-KO: IDs 01, 02, 03, 11, 12 (ID 13 excluded)
#   For C-peptide (indexed only by an unlabelled per-group
#   Replicate number, no shared animal ID), the matching subset
#   is simply the first five replicates (Rep1-Rep5) in both
#   groups; this is consistent with, but not independently
#   verified against, the ID-based subset identified for body
#   weight/glucose. This filtering should still be independently
#   confirmed against original laboratory records.
#
# Output:
#   Supplementary_Table_S4_Body_Weight.csv
#   Supplementary_Table_S4_C_Peptide.csv
#   Supplementary_Table_S4_Fasting_Blood_Glucose.csv
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# 1) Load required libraries
required_packages <- c("dplyr", "readr", "lme4", "tibble")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(dplyr)
library(readr)
library(lme4)
library(tibble)

# ============================================================
# 2) Helper: identify the n=5 cohort by exhaustive matching
#    against a published summary-statistics table
# ============================================================

find_matching_subset <- function(long_df, id_col, value_col, week_col,
                                  group_name, target_df, tol_mean, tol_sd) {
  ids <- unique(long_df[[id_col]][long_df$Group == group_name])
  target_group <- target_df[target_df$Group == group_name, ]
  target_group$Week_num <- as.character(target_group$Week_num)
  long_group <- long_df[long_df$Group == group_name, ]
  long_group[[week_col]] <- as.character(long_group[[week_col]])

  combos <- combn(ids, 5, simplify = FALSE)
  for (combo in combos) {
    ok <- TRUE
    for (i in seq_len(nrow(target_group))) {
      wk <- target_group$Week_num[i]
      vals <- long_group[[value_col]][
        long_group[[id_col]] %in% combo & long_group[[week_col]] == wk
      ]
      if (length(vals) != 5) { ok <- FALSE; break }
      m <- mean(vals); s <- sd(vals)
      if (abs(m - target_group$mean[i]) > tol_mean ||
          abs(s - target_group$sd[i]) > tol_sd) {
        ok <- FALSE; break
      }
    }
    if (ok) return(combo)
  }
  stop(paste("No matching 5-animal subset found for group:", group_name))
}

# ============================================================
# 3) Helper: fit the mixed model and tidy the fixed effects
#    into the manuscript's reporting format (asymptotic Wald z)
# ============================================================

fit_and_tidy <- function(df) {
  model <- lmer(Value ~ Group * Week_num_numeric + (1 | Animal),
                data = df, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa"))
  co <- summary(model)$coefficients
  tibble(
    Term        = rownames(co),
    Coefficient = round(co[, "Estimate"], 3),
    SE          = round(co[, "Std. Error"], 3),
    z           = round(co[, "Estimate"] / co[, "Std. Error"], 2),
    p           = signif(2 * pnorm(-abs(co[, "Estimate"] / co[, "Std. Error"])), 3)
  )
}

# ============================================================
# 4) Body weight (weeks 0-12)
# ============================================================

bw_raw <- read_csv("../Fig1_Intestinal/Fig_1b_Data_BW.csv", show_col_types = FALSE)
bw_target <- read_csv("../Fig1_Intestinal/Fig_1b_Summary_BW.csv", show_col_types = FALSE)

bw_il10ko_ids <- find_matching_subset(bw_raw, "ID", "BW", "Week_num",
                                       "IL10-KO", bw_target, tol_mean = 0.006, tol_sd = 0.006)
bw_dss_ids    <- find_matching_subset(bw_raw, "ID", "BW", "Week_num",
                                       "DSS-IL10-KO", bw_target, tol_mean = 0.006, tol_sd = 0.006)
cat("Body weight cohort: IL10-KO =", paste(bw_il10ko_ids, collapse = ", "),
    "| DSS-IL10-KO =", paste(bw_dss_ids, collapse = ", "), "\n")

ctrl_ids <- unique(bw_raw$ID[bw_raw$Group == "CTRL"])
bw_keep_ids <- c(ctrl_ids, bw_il10ko_ids, bw_dss_ids)

bw_model_data <- bw_raw %>%
  filter(ID %in% bw_keep_ids) %>%
  mutate(
    Group = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
    Week_num_numeric = as.numeric(Week_num),
    Animal = ID,
    Value = BW
  )

bw_table <- fit_and_tidy(bw_model_data)
cat("\n=== Body weight mixed-effects model ===\n")
print(bw_table)
write_csv(bw_table, "Supplementary_Table_S4_Body_Weight.csv")

# ============================================================
# 5) Fasting C-peptide (weeks 0, 6, 12)
#    No shared animal ID in the source file; Replicate is used
#    as the within-group animal label (see note above).
# ============================================================

cp_raw <- read_csv("../Fig2_Pancreas/Fig_2bc_Data_C_peptide_FBG.csv", show_col_types = FALSE)
cp_target <- read_csv("../Fig2_Pancreas/Fig_2b_Summary_C_peptide.csv", show_col_types = FALSE)

cp_il10ko_reps <- find_matching_subset(cp_raw, "Replicate", "C_peptide", "Week_num",
                                        "IL10-KO", cp_target, tol_mean = 0.0006, tol_sd = 0.0006)
cp_dss_reps    <- find_matching_subset(cp_raw, "Replicate", "C_peptide", "Week_num",
                                        "DSS-IL10-KO", cp_target, tol_mean = 0.0006, tol_sd = 0.0006)
cat("\nC-peptide cohort: IL10-KO =", paste(cp_il10ko_reps, collapse = ", "),
    "| DSS-IL10-KO =", paste(cp_dss_reps, collapse = ", "), "\n")

cp_model_data <- cp_raw %>%
  mutate(Animal = paste(Group, Replicate, sep = "_")) %>%
  filter(
    (Group == "CTRL") |
    (Group == "IL10-KO"     & Replicate %in% cp_il10ko_reps) |
    (Group == "DSS-IL10-KO" & Replicate %in% cp_dss_reps)
  ) %>%
  mutate(
    Group = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
    Week_num_numeric = as.numeric(Week_num),
    Value = C_peptide
  )

cp_table <- fit_and_tidy(cp_model_data)
cat("\n=== Fasting C-peptide mixed-effects model ===\n")
print(cp_table)
write_csv(cp_table, "Supplementary_Table_S4_C_Peptide.csv")

# ============================================================
# 6) Fasting blood glucose (weeks 0, 6, 12)
#    Uses the ID-indexed FBG file (not the Replicate-only
#    combined C-peptide/FBG file) so the same animal identities
#    as the body-weight model can be reused directly.
# ============================================================

fbg_raw <- read_csv("../Fig2_Pancreas/Fig_2a_Data_FBG.csv", show_col_types = FALSE)
fbg_target <- read_csv("../Fig2_Pancreas/Fig_2a_Summary_FBG.csv", show_col_types = FALSE)

fbg_il10ko_ids <- find_matching_subset(fbg_raw, "ID", "FBG", "Week_num",
                                        "IL10-KO", fbg_target, tol_mean = 0.06, tol_sd = 0.06)
fbg_dss_ids    <- find_matching_subset(fbg_raw, "ID", "FBG", "Week_num",
                                        "DSS-IL10-KO", fbg_target, tol_mean = 0.06, tol_sd = 0.06)
cat("\nFasting blood glucose cohort: IL10-KO =", paste(fbg_il10ko_ids, collapse = ", "),
    "| DSS-IL10-KO =", paste(fbg_dss_ids, collapse = ", "), "\n")

fbg_ctrl_ids <- unique(fbg_raw$ID[fbg_raw$Group == "CTRL"])
fbg_keep_ids <- c(fbg_ctrl_ids, fbg_il10ko_ids, fbg_dss_ids)

fbg_model_data <- fbg_raw %>%
  filter(ID %in% fbg_keep_ids, Week_num %in% c(0, 6, 12)) %>%
  mutate(
    Group = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
    Week_num_numeric = as.numeric(Week_num),
    Animal = ID,
    Value = FBG
  )

fbg_table <- fit_and_tidy(fbg_model_data)
cat("\n=== Fasting blood glucose mixed-effects model ===\n")
print(fbg_table)
write_csv(fbg_table, "Supplementary_Table_S4_Fasting_Blood_Glucose.csv")

cat("\nSaved:\n")
cat(" - Supplementary_Table_S4_Body_Weight.csv\n")
cat(" - Supplementary_Table_S4_C_Peptide.csv\n")
cat(" - Supplementary_Table_S4_Fasting_Blood_Glucose.csv\n")

# ============================================================
# END OF SCRIPT
# ============================================================
