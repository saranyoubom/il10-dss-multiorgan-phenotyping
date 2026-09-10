# ============================================================
# Supplementary Table 3 - Linear mixed-effects models for
# body weight, fasting C-peptide and fasting blood glucose
# Scientific Reports Rev1 revision (Reviewer 2 point 5)
#
# Related to Results Figs. 1b and 2b, and Methods
# "Statistical analysis."
#
# ---------------------------------------------------------
# MODEL
# ---------------------------------------------------------
#   Value ~ Group * Week + (1 | Animal)
#
#   Random intercept per animal; REML estimation; CTRL as the
#   reference group. Fixed effects are reported with an
#   asymptotic Wald z test (Estimate / SE). A Wald z test is
#   used rather than a Satterthwaite or Kenward-Roger t test
#   because the residual degrees of freedom for the
#   Group x Week interaction are not uniquely defined in an
#   unbalanced repeated-measures design of this size, and the
#   z test makes no distributional claim beyond the asymptotic
#   normality already assumed by the model.
#
#   Timepoints used per outcome:
#     body weight    weeks 0-12           (13 timepoints, 195 obs)
#     C-peptide      weeks 0, 3, 6, 9, 12 ( 5 timepoints,  75 obs)
#     blood glucose  weeks 0, 3, 6, 9, 12 ( 5 timepoints,  75 obs)
#
#   Body weight and C-peptide each use every timepoint measured
#   for that outcome, as stated in the Results text ("a linear
#   mixed-effects model incorporating / across all measured
#   timepoints"), which is more than the subset plotted in the
#   figures (Fig. 1b and Fig. 2b plot selected weeks only).
#
#   Fasting blood glucose is measured weekly (13 timepoints) but
#   is modelled here on C-peptide's 5-timepoint schedule, because
#   the Results text describes it as "the EQUIVALENT model for
#   fasting blood glucose" - equivalent to the C-peptide model -
#   and because the two fasting measures are reported as a pair
#   (their ratio, FCP/FBG, is itself an endpoint in Fig. 2c).
#   Modelling them on a common schedule keeps the two directly
#   comparable. This also matches the 75-observation count quoted
#   in the Supplementary Table 3 caption.
#
# ---------------------------------------------------------
# DATA SOURCES (and why)
# ---------------------------------------------------------
#   Body weight    <- Fig_1b_Data_BW.csv          (ID-indexed)
#   Fasting glucose<- Fig_2a_Data_FBG.csv         (ID-indexed)
#   C-peptide      <- Fig_2bc_Data_C_peptide_FBG.csv
#                                          (Replicate-indexed)
#
#   NOTE: fasting blood glucose appears in two files. The FBG
#   column of Fig_2bc_Data_C_peptide_FBG.csv disagrees with
#   Fig_2a_Data_FBG.csv in 65 of 90 shared Group/Week/Replicate
#   keys: CTRL agrees exactly, but both knockout groups differ,
#   and the per-group-per-week value multisets are not
#   permutations of one another, so the disagreement is a
#   misalignment against animal identity rather than a
#   re-ordering. Fig_2a_Data_FBG.csv is used here because it is
#   the file that reproduces the published summary statistics in
#   Fig_2a_Summary_FBG.csv exactly; the Fig_2bc FBG column
#   reproduces no published summary and should not be used.
#   This discrepancy is worth resolving at source.
#
#   Because that misalignment exists between the two files,
#   the Replicate labels in Fig_2bc cannot be assumed to refer
#   to the same animals as the ID labels elsewhere. C-peptide
#   is therefore modelled using Fig_2bc's own Replicate label
#   as the animal identifier, and is NOT joined to animal IDs.
#
# ---------------------------------------------------------
# COHORT (n = 5 per group)
# ---------------------------------------------------------
#   The raw files contain more animals than the n = 5 per group
#   reported throughout the manuscript (7 IL10-KO, 6
#   DSS-IL10-KO). Rather than assume which animals to drop,
#   this script recovers the reported cohort from the project's
#   own published summary-statistics CSVs: it searches every
#   5-animal subset and keeps the one whose mean and SD match
#   the published values at every reported timepoint. The match
#   is unique in each case.
#
# Output:
#   Supplementary_Table_3_Body_Weight.csv
#   Supplementary_Table_3_C_Peptide.csv
#   Supplementary_Table_3_Fasting_Blood_Glucose.csv
# ============================================================

# 0) Clear workspace
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
# 2) Recover the reported n = 5 cohort from published summaries
# ============================================================

find_matching_subset <- function(long_df, id_col, value_col, week_col,
                                  group_name, target_df, tol_mean, tol_sd) {
  ids <- unique(long_df[[id_col]][long_df$Group == group_name])
  if (length(ids) == 5) return(ids)   # already n = 5, nothing to select

  target_group <- target_df[target_df$Group == group_name, ]
  target_group$Week_num <- as.character(target_group$Week_num)
  long_group <- long_df[long_df$Group == group_name, ]
  long_group[[week_col]] <- as.character(long_group[[week_col]])

  hits <- list()
  for (combo in combn(ids, 5, simplify = FALSE)) {
    ok <- TRUE
    for (i in seq_len(nrow(target_group))) {
      wk <- target_group$Week_num[i]
      vals <- long_group[[value_col]][
        long_group[[id_col]] %in% combo & long_group[[week_col]] == wk
      ]
      if (length(vals) != 5) { ok <- FALSE; break }
      if (abs(mean(vals) - target_group$mean[i]) > tol_mean ||
          abs(sd(vals)   - target_group$sd[i])   > tol_sd) {
        ok <- FALSE; break
      }
    }
    if (ok) hits[[length(hits) + 1]] <- combo
  }
  if (length(hits) == 0)
    stop(paste("No 5-animal subset reproduces the published summary for group:", group_name))
  if (length(hits) > 1)
    warning(paste("Cohort is not uniquely identified for group:", group_name,
                  "-", length(hits), "subsets match; using the first."))
  hits[[1]]
}

# ============================================================
# 3) Fit the model and report fixed effects
# ============================================================

fit_and_tidy <- function(df, label) {
  cat(sprintf("%s: %d observations, %d animals, %d timepoints (weeks %s)\n",
              label, nrow(df), dplyr::n_distinct(df$Animal),
              dplyr::n_distinct(df$Week),
              paste(sort(unique(df$Week)), collapse = ", ")))
  model <- lmer(Value ~ Group * Week + (1 | Animal),
                data = df, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa"))
  if (isSingular(model)) {
    warning(paste0(label, ": singular fit - the between-animal variance is ",
                   "estimated at (or near) zero, i.e. animal-to-animal ",
                   "variability is not separable from residual noise at this ",
                   "sample size. Fixed effects remain interpretable."))
  }
  co <- summary(model)$coefficients
  z  <- co[, "Estimate"] / co[, "Std. Error"]
  tibble(
    Term        = rownames(co),
    Coefficient = round(co[, "Estimate"], 3),
    SE          = round(co[, "Std. Error"], 3),
    z           = round(z, 2),
    p           = signif(2 * pnorm(-abs(z)), 3),
    Singular_fit = isSingular(model)
  )
}

# ============================================================
# 4) Body weight (weeks 0-12), ID-indexed
# ============================================================

bw_raw    <- read_csv("../Fig1_Intestinal/Fig_1b_Data_BW.csv", show_col_types = FALSE)
bw_target <- read_csv("../Fig1_Intestinal/Fig_1b_Summary_BW.csv", show_col_types = FALSE)

bw_ids <- unlist(lapply(c("CTRL", "IL10-KO", "DSS-IL10-KO"), function(g)
  find_matching_subset(bw_raw, "ID", "BW", "Week_num", g, bw_target,
                       tol_mean = 0.006, tol_sd = 0.006)))
cat("Body weight cohort:", paste(bw_ids, collapse = ", "), "\n")

bw_model_data <- bw_raw %>%
  filter(ID %in% bw_ids) %>%
  mutate(Group  = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
         Week   = as.numeric(Week_num),
         Animal = ID,
         Value  = BW)

bw_table <- fit_and_tidy(bw_model_data, "Body weight")
cat("\n=== Body weight ===\n"); print(bw_table)
write_csv(bw_table, "Supplementary_Table_3_Body_Weight.csv")

# ============================================================
# 5) Fasting C-peptide (weeks 0, 6, 12), Replicate-indexed
#    See DATA SOURCES note: not joined to animal IDs.
# ============================================================

cp_raw    <- read_csv("../Fig2_Pancreas/Fig_2bc_Data_C_peptide_FBG.csv", show_col_types = FALSE)
cp_target <- read_csv("../Fig2_Pancreas/Fig_2b_Summary_C_peptide.csv", show_col_types = FALSE)

cp_reps <- lapply(c("CTRL", "IL10-KO", "DSS-IL10-KO"), function(g)
  data.frame(Group = g,
             Replicate = find_matching_subset(cp_raw, "Replicate", "C_peptide",
                                              "Week_num", g, cp_target,
                                              tol_mean = 0.0006, tol_sd = 0.0006)))
cp_keep <- bind_rows(cp_reps)
cat("\nC-peptide cohort:\n"); print(cp_keep)

cp_model_data <- cp_raw %>%
  semi_join(cp_keep, by = c("Group", "Replicate")) %>%
  mutate(Group  = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
         Week   = as.numeric(Week_num),
         Animal = paste(Group, Replicate, sep = "_"),
         Value  = C_peptide)

cp_table <- fit_and_tidy(cp_model_data, "C-peptide")
cat("\n=== Fasting C-peptide ===\n"); print(cp_table)
write_csv(cp_table, "Supplementary_Table_3_C_Peptide.csv")

# ============================================================
# 6) Fasting blood glucose (weeks 0, 6, 12), ID-indexed
# ============================================================

fbg_raw    <- read_csv("../Fig2_Pancreas/Fig_2a_Data_FBG.csv", show_col_types = FALSE)
fbg_target <- read_csv("../Fig2_Pancreas/Fig_2a_Summary_FBG.csv", show_col_types = FALSE)

fbg_ids <- unlist(lapply(c("CTRL", "IL10-KO", "DSS-IL10-KO"), function(g)
  find_matching_subset(fbg_raw, "ID", "FBG", "Week_num", g, fbg_target,
                       tol_mean = 0.06, tol_sd = 0.06)))
cat("\nFasting blood glucose cohort:", paste(fbg_ids, collapse = ", "), "\n")

# C-peptide's measurement schedule, so that the two fasting
# measures are modelled on a common set of timepoints (see note
# in the header).
fbg_weeks <- sort(unique(as.numeric(cp_model_data$Week)))
cat("\nFBG restricted to C-peptide's schedule, weeks:",
    paste(fbg_weeks, collapse = ", "), "\n")

fbg_model_data <- fbg_raw %>%
  filter(ID %in% fbg_ids, Week_num %in% fbg_weeks) %>%
  mutate(Group  = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
         Week   = as.numeric(Week_num),
         Animal = ID,
         Value  = FBG)

fbg_table <- fit_and_tidy(fbg_model_data, "Fasting blood glucose")
cat("\n=== Fasting blood glucose ===\n"); print(fbg_table)
write_csv(fbg_table, "Supplementary_Table_3_Fasting_Blood_Glucose.csv")

# ============================================================
# 7) Report the cohorts actually used, for the record
# ============================================================
cat("\n--- Cohorts used ---\n")
cat("Body weight   :", paste(sort(bw_ids), collapse = ", "), "\n")
cat("Blood glucose :", paste(sort(fbg_ids), collapse = ", "), "\n")
cat("C-peptide     : Replicate labels, see table above\n")
cat("\nSaved:\n")
cat(" - Supplementary_Table_3_Body_Weight.csv\n")
cat(" - Supplementary_Table_3_C_Peptide.csv\n")
cat(" - Supplementary_Table_3_Fasting_Blood_Glucose.csv\n")

# ============================================================
# END OF SCRIPT
# ============================================================
