# ============================================================
# Supplementary Table S5 - Post-hoc effect sizes and
# statistical power for representative group comparisons
# Scientific Reports Rev1 rebuttal (Reviewer 2 points 1/3;
# Reviewer 3 point 2)
#
# Related to Discussion "Limitations."
#
# Method:
#   Cohen's d is computed from each comparison's reported group
#   means and SDs using the pooled-SD formula:
#     d = (mean1 - mean2) / sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2))
#   Achieved power is computed for a two-sample, two-tailed
#   t-test at alpha = 0.05 using base R's power.t.test(), which
#   implements the standard noncentral-t-distribution power
#   calculation (equivalent to Python's
#   statsmodels.stats.power.TTestIndPower used for the original
#   rebuttal numbers, since R/Rscript was unavailable on the
#   machine that produced them).
#
# The four representative comparisons are read directly from
# this project's own summary-statistics CSVs rather than
# retyped, so the group means/SDs stay tied to their source
# figures.
#
# Output:
#   Supplementary_Table_S5_Effect_Sizes_Power.csv
#   Supplementary_Table_S5_Minimum_Detectable_Effect.csv
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# 1) Load required libraries
required_packages <- c("dplyr", "readr", "tibble")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(dplyr)
library(readr)
library(tibble)

# ============================================================
# 2) Helper functions
# ============================================================

cohens_d <- function(m1, s1, n1, m2, s2, n2) {
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  (m1 - m2) / pooled_sd
}

achieved_power <- function(d, n1, n2, alpha = 0.05) {
  power.t.test(n = mean(c(n1, n2)), delta = abs(d), sd = 1,
               sig.level = alpha, type = "two.sample",
               alternative = "two.sided")$power
}

# ============================================================
# 3) Representative comparisons, read from source summaries
# ============================================================

fcp_fbg <- read_csv("../Fig2_Pancreas/Fig_2c_Summary_FCP_FBG_Ratio.csv", show_col_types = FALSE)
ag      <- read_csv("../Fig1_Intestinal/Fig_1d_Summary_AG_Ratio.csv", show_col_types = FALSE)
spleen  <- read_csv("../Fig3_Immune/Fig_3c_Summary_FoldChange_Spleen_GeneExpr.csv", show_col_types = FALSE)
liver   <- read_csv("../Fig4_Organs/Fig_4b_Summary_FoldChange_Liver_GeneExpr.csv", show_col_types = FALSE)

get_row <- function(df, group, gene = NULL) {
  if (!is.null(gene)) df <- df %>% filter(Gene == gene)
  df %>% filter(Group == group)
}

comparisons <- list(
  list(name = "FCP/FBG ratio, CTRL vs IL10-KO", figure = "Fig. 2c",
       g1 = get_row(fcp_fbg, "CTRL"), g2 = get_row(fcp_fbg, "IL10-KO"),
       mean_col = "mean", sd_col = "sd"),
  list(name = "Albumin/globulin ratio, CTRL vs DSS-IL10-KO", figure = "Fig. 1d",
       g1 = get_row(ag, "CTRL"), g2 = get_row(ag, "DSS-IL10-KO"),
       mean_col = "mean", sd_col = "sd"),
  list(name = "Splenic Tcf19 log2FC, CTRL vs DSS-IL10-KO", figure = "Fig. 3c",
       g1 = get_row(spleen, "CTRL", "Tcf19"), g2 = get_row(spleen, "DSS-IL10-KO", "Tcf19"),
       mean_col = "mean_Log2FC", sd_col = "sd_Log2FC"),
  list(name = "Hepatic Zfand6 log2FC, CTRL vs DSS-IL10-KO", figure = "Fig. 4b",
       g1 = get_row(liver, "CTRL", "Zfand6"), g2 = get_row(liver, "DSS-IL10-KO", "Zfand6"),
       mean_col = "mean_Log2FC", sd_col = "sd_Log2FC")
)

es_table <- do.call(rbind, lapply(comparisons, function(cmp) {
  m1 <- cmp$g1[[cmp$mean_col]]; s1 <- cmp$g1[[cmp$sd_col]]; n1 <- cmp$g1$n
  m2 <- cmp$g2[[cmp$mean_col]]; s2 <- cmp$g2[[cmp$sd_col]]; n2 <- cmp$g2$n
  d <- cohens_d(m1, s1, n1, m2, s2, n2)
  pw <- achieved_power(d, n1, n2)
  tibble(
    Comparison   = cmp$name,
    Figure       = cmp$figure,
    `n (group 1)` = n1,
    `n (group 2)` = n2,
    `Cohen's d`   = round(d, 2),
    `Achieved power` = ifelse(pw > 0.999, ">0.999", round(pw, 3))
  )
}))

cat("=== Effect sizes and achieved power ===\n")
print(es_table)
write_csv(es_table, "Supplementary_Table_S5_Effect_Sizes_Power.csv")

# ============================================================
# 4) Minimum detectable effect size at 80% power
# ============================================================

mdes_table <- do.call(rbind, lapply(c(3, 4, 5), function(n) {
  d_min <- power.t.test(n = n, power = 0.80, sig.level = 0.05,
                         type = "two.sample", alternative = "two.sided")$delta
  tibble(`n per group` = n, `Minimum detectable Cohen's d` = round(d_min, 2))
}))

cat("\n=== Minimum detectable effect size (80% power, alpha = 0.05, two-tailed) ===\n")
print(mdes_table)
write_csv(mdes_table, "Supplementary_Table_S5_Minimum_Detectable_Effect.csv")

cat("\nSaved:\n")
cat(" - Supplementary_Table_S5_Effect_Sizes_Power.csv\n")
cat(" - Supplementary_Table_S5_Minimum_Detectable_Effect.csv\n")

# ============================================================
# END OF SCRIPT
# ============================================================
