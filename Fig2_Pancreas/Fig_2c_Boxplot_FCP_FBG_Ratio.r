# ============================================================
# Fasting C-peptide (FCP)/FBG Ratio Analysis - Enhanced Version
# Revised for Diabetologia Journal Submission
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Between-group comparison (independent samples) at Week 12
# ============================================================

# STUDY DESIGN NOTE:
# This is a BETWEEN-GROUP comparison at endpoint (Week 12)
# - Independent samples from 3 different groups
# - Single measurement per animal at Week 12
# - Comparing groups at one time point
# - FCP/FBG ratio reflects beta-cell function relative to glucose
# - Test selected by normality (consistent between plot and CSV)

# 0) Clear workspace and set working directory
rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# 1) Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggpubr)
library(rstatix)

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Read the FCP and FBG data
cat("\n=== LOADING FASTING C-PEPTIDE (FCP) AND FBG DATA ===\n")
dat <- read_csv("Fig_2bc_Data_C_peptide_FBG.csv")
dat <- dat %>% filter(Replicate %in% c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5"))

# Rename C_peptide to FCP for clarity
dat <- dat %>% rename(FCP = C_peptide)

cat("\nData structure:\n")
print(head(dat, 10))
cat("\nWeeks available:", unique(dat$Week_num), "\n")
cat("Groups:", unique(dat$Group), "\n")

# ============================================================
# 4) Calculate FCP/FBG Ratio at Week 12
# ============================================================
cat("\n=== CALCULATING FCP/FBG RATIO AT WEEK 12 ===\n")

ratio_data <- dat %>%
  filter(Week_num == 12) %>%
  filter(!is.na(FCP) & !is.na(FBG) & FBG != 0) %>%
  mutate(
    FCP_FBG_Ratio = FCP / FBG,
    Sample_ID     = paste(Group, Replicate, sep = "_")
  ) %>%
  select(Group, Replicate, FBG, FCP, FCP_FBG_Ratio, Sample_ID)

cat("\nCalculated FCP/FBG Ratios at Week 12:\n")
print(ratio_data %>% select(Group, Replicate, FBG, FCP, FCP_FBG_Ratio))

ratio_data$Group <- factor(ratio_data$Group,
                           levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))

write.csv(ratio_data, "Fig_2c_Data_FCP_FBG_Ratio.csv", row.names = FALSE)
cat("✓ FCP/FBG ratio data exported to Fig_2c_Data_FCP_FBG_Ratio.csv\n")

cat("\nSample sizes per group at Week 12:\n")
print(table(ratio_data$Group))

# ============================================================
# 5) NORMALITY TESTING
# Drives test selection — same test used in both plot and CSV
# ============================================================
cat("\n=== NORMALITY TESTING (Shapiro-Wilk) ===\n")
cat("Testing each group separately\n")
cat("(Selected test applied consistently in plot and exported CSV)\n\n")

normality_results <- data.frame()

for (grp in levels(ratio_data$Group)) {
  data_subset <- ratio_data %>%
    filter(Group == grp) %>%
    pull(FCP_FBG_Ratio)
  
  if (length(data_subset) >= 3) {
    shapiro_result <- shapiro.test(data_subset)
    is_normal      <- shapiro_result$p.value > 0.05
    
    normality_results <- rbind(normality_results, data.frame(
      Group     = grp,
      n         = length(data_subset),
      Shapiro_W = round(shapiro_result$statistic, 4),
      Shapiro_p = round(shapiro_result$p.value, 4),
      Normal    = is_normal
    ))
    
    cat(sprintf("%s: n=%d, W = %.4f, p = %.4f \u2192 %s\n",
                grp, length(data_subset),
                shapiro_result$statistic,
                shapiro_result$p.value,
                ifelse(is_normal, "Normal", "Non-normal")))
  }
}

cat("\n=== NORMALITY RESULTS SUMMARY ===\n")
print(normality_results)

write.csv(normality_results,
          file      = "Fig_2c_Normality_Results.csv",
          row.names = FALSE)
cat("✓ Normality results exported to Fig_2c_Normality_Results.csv\n")

# ============================================================
# 6) DETERMINE OVERALL TEST STRATEGY
# ============================================================
cat("\n=== TEST SELECTION STRATEGY ===\n")

all_normal <- all(normality_results$Normal)

if (all_normal) {
  overall_test  <- "ANOVA"
  pairwise_test <- "t.test"
  cat("Decision: ALL groups are normally distributed\n")
  cat("\u2192 Use parametric approach:\n")
  cat("  - Overall test: One-way ANOVA\n")
  cat("  - Pairwise tests: Unpaired t-tests\n")
} else {
  overall_test  <- "Kruskal-Wallis"
  pairwise_test <- "wilcox.test"
  cat("Decision: At least one group is NOT normally distributed\n")
  cat("\u2192 Use non-parametric approach:\n")
  cat("  - Overall test: Kruskal-Wallis\n")
  cat("  - Pairwise tests: Mann-Whitney U tests\n")
}

test_decision <- data.frame(
  All_Groups_Normal = all_normal,
  Overall_Test      = overall_test,
  Pairwise_Test     = pairwise_test
)

write.csv(test_decision,
          file      = "Fig_2c_Test_Decision.csv",
          row.names = FALSE)
cat("✓ Test decision exported to Fig_2c_Test_Decision.csv\n")

# ============================================================
# 7) PERFORM OVERALL TEST
# ============================================================
cat("\n=== PERFORMING OVERALL TEST ===\n")

if (all_normal) {
  overall_result  <- aov(FCP_FBG_Ratio ~ Group, data = ratio_data)
  overall_summary <- summary(overall_result)
  cat("\nOne-way ANOVA:\n")
  print(overall_summary)
  overall_p <- overall_summary[[1]][["Pr(>F)"]][1]
  cat(sprintf("\nF-statistic: %.4f, p-value: %.4f\n",
              overall_summary[[1]][["F value"]][1], overall_p))
} else {
  overall_result <- kruskal.test(FCP_FBG_Ratio ~ Group, data = ratio_data)
  cat("\nKruskal-Wallis test:\n")
  cat(sprintf("chi-squared = %.4f, df = %d, p-value = %.4f\n",
              overall_result$statistic,
              overall_result$parameter,
              overall_result$p.value))
  overall_p <- overall_result$p.value
}

if (overall_p < 0.05) {
  cat("\nOverall test is SIGNIFICANT \u2192 proceed with pairwise comparisons\n")
} else {
  cat("\nOverall test is NOT significant \u2192 pairwise tests for completeness only\n")
}

# ============================================================
# 8) PAIRWISE POST-HOC COMPARISONS — RAW P-VALUES, NO CORRECTION
#    Same test method used here as in the plot (pairwise_test)
# ============================================================
cat("\n=== PERFORMING PAIRWISE COMPARISONS (raw p-values, no correction) ===\n")

my_comparisons <- list(
  c("CTRL", "IL10-KO"),
  c("CTRL", "DSS-IL10-KO"),
  c("IL10-KO", "DSS-IL10-KO")
)

posthoc_result <- compare_means(
  formula     = FCP_FBG_Ratio ~ Group,
  data        = ratio_data,
  method      = pairwise_test,   # Adaptive — consistent with plot
  paired      = FALSE,
  comparisons = my_comparisons
) %>%
  rename(Variable = .y.)

# Format results — raw p-values only
posthoc_result <- posthoc_result %>%
  mutate(
    p        = signif(p, 3),
    p.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    Test_Method = pairwise_test
  ) %>%
  select(Variable, group1, group2, Test_Method, p, p.signif)

cat("\n=== Pairwise Comparison Results (raw p-values) ===\n")
print(posthoc_result)

write.csv(posthoc_result,
          file      = "Fig_2c_Statistical_Tests.csv",
          row.names = FALSE)
cat("✓ Statistical test results exported to Fig_2c_Statistical_Tests.csv\n")

# ============================================================
# 9) Summary statistics by Group
# ============================================================
cat("\n=== CALCULATING SUMMARY STATISTICS ===\n")

summary_stats <- ratio_data %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    mean   = round(mean(FCP_FBG_Ratio,   na.rm = TRUE), 4),
    sd     = round(sd(FCP_FBG_Ratio,     na.rm = TRUE), 4),
    sem    = round(sd(FCP_FBG_Ratio,     na.rm = TRUE) / sqrt(n()), 4),
    median = round(median(FCP_FBG_Ratio, na.rm = TRUE), 4),
    min    = round(min(FCP_FBG_Ratio,    na.rm = TRUE), 4),
    max    = round(max(FCP_FBG_Ratio,    na.rm = TRUE), 4),
    .groups = 'drop'
  ) %>%
  mutate(Mean_SEM = paste0(mean, " \u00b1 ", sem))

cat("\n=== Summary Statistics (Mean \u00b1 SEM) ===\n")
print(summary_stats)

write.csv(summary_stats,
          file      = "Fig_2c_Summary_FCP_FBG_Ratio.csv",
          row.names = FALSE)
cat("✓ Summary statistics exported to Fig_2c_Summary_FCP_FBG_Ratio.csv\n")

# ============================================================
# 10) Create boxplot
#     pairwise_test ensures plot matches exported CSV exactly
# ============================================================
cat("\n=== GENERATING BOXPLOT ===\n")

x_labels <- c("CTRL", "IL10-KO", "DSS-\nIL10-KO")
names(x_labels) <- c("CTRL", "IL10-KO", "DSS-IL10-KO")

boxplot_ratio <- ggboxplot(
  ratio_data, x = "Group", y = "FCP_FBG_Ratio",
  color  = "Group", fill = "Group",
  alpha  = 0.5, add = "jitter",
  add.params = list(color = "Group", shape = 21,
                    fill = "white", size = 0.75),
  palette       = diabetologia_colors,
  outlier.shape = NA
) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02), limits = c(0, 0.02)) +
  labs(
    y = "FCP/FBG ratio (AU)",
    x = "Experimental group"
  ) +
  scale_x_discrete(labels = x_labels) +
  theme_classic() +
  theme(
    legend.position  = "none",
    axis.title.y     = element_text(size = 8),
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 8, colour = "black"),
    axis.text.y      = element_text(size = 8, colour = "black",
                                    angle = 90, hjust = 0.5),
    axis.line        = element_line(colour = "black",
                                    linewidth = 0.75 / .pt),
    axis.ticks       = element_line(colour = "black",
                                    linewidth = 0.75 / .pt),
    plot.margin      = margin(3.8, 3.8, 3.8, 3.8, "mm")
  ) +
  stat_compare_means(
    comparisons   = my_comparisons,
    method        = pairwise_test,   # Consistent with CSV export
    paired        = FALSE,
    tip.length    = 0.0,
    step.increase = 0.25,
    label         = "p.signif",
    size          = 2.5,
    bracket.size  = 0.3
  )

print(boxplot_ratio)

ggsave(
  filename  = "Fig_2c_Boxplot_FCP_FBG_Ratio.png",
  plot      = boxplot_ratio,
  scale = 1, width = 5.9, height = 5.7,
  units = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename  = "Fig_2c_Boxplot_FCP_FBG_Ratio.svg",
  plot      = boxplot_ratio,
  scale = 1, width = 5.9, height = 5.7,
  units = "cm", limitsize = FALSE
)

cat("✓ Boxplot saved in PNG and SVG formats\n")

# ============================================================
# 11) Biological interpretation summary
# ============================================================
cat("\n=== BIOLOGICAL INTERPRETATION ===\n")

bio_summary <- ratio_data %>%
  group_by(Group) %>%
  summarise(
    n          = n(),
    Mean_FCP   = round(mean(FCP,           na.rm = TRUE), 3),
    Mean_FBG   = round(mean(FBG,           na.rm = TRUE), 1),
    Mean_Ratio = round(mean(FCP_FBG_Ratio, na.rm = TRUE), 4),
    .groups = 'drop'
  )

cat("\nMean values by group at Week 12:\n")
print(bio_summary)

write.csv(bio_summary,
          file      = "Fig_2c_Biological_Summary.csv",
          row.names = FALSE)

cat("\nFCP/FBG ratio interpretation:\n")
cat("  Higher ratio: Better beta-cell compensation relative to glucose levels\n")
cat("  Lower ratio: Impaired beta-cell function relative to glucose demand\n")
cat("  FCP reflects insulin secretion capacity\n")
cat("  FBG reflects glucose homeostasis\n")
cat("  Ratio normalises insulin secretion to glycaemic state\n")
cat("  Units: (ng/mL) / (mg/dL) = dimensionless\n")

# ============================================================
# 12) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat(" 1. Fig_2c_Normality_Results.csv       (Shapiro-Wilk per group)\n")
cat(" 2. Fig_2c_Test_Decision.csv           (test selection strategy)\n")
cat(" 3. Fig_2c_Statistical_Tests.csv       (pairwise tests, raw p-values)\n")
cat(" 4. Fig_2c_Boxplot_FCP_FBG_Ratio.png  (600 dpi)\n")
cat(" 5. Fig_2c_Boxplot_FCP_FBG_Ratio.svg  (vector format)\n")
cat(" 6. Fig_2c_Data_FCP_FBG_Ratio.csv     (calculated ratios)\n")
cat(" 7. Fig_2c_Summary_FCP_FBG_Ratio.csv  (descriptive statistics)\n")
cat(" 8. Fig_2c_Biological_Summary.csv     (biological interpretation data)\n")
cat("\nFigure dimensions: 5.9 cm (width) \u00d7 5.7 cm (height)\n")
cat("\nStatistical approach:\n")
cat(" - Design: BETWEEN-GROUP comparison at Week 12 (endpoint)\n")
cat(" - Normality tested separately for each group (Shapiro-Wilk)\n")
cat(sprintf(" - Overall test: %s\n", overall_test))
cat(sprintf(" - Pairwise tests: %s (consistent in plot and CSV)\n", pairwise_test))
cat(" - Multiple comparison correction: NONE (raw p-values reported)\n")
cat("\nComparisons:\n")
cat(" - CTRL vs IL10-KO\n")
cat(" - CTRL vs DSS-IL10-KO\n")
cat(" - IL10-KO vs DSS-IL10-KO\n")
cat("   Total: 3 pairwise comparisons\n")
cat("\nSignificance levels (raw p-values):\n")
cat(" *** p<0.001\n")
cat("  ** p<0.01\n")
cat("   * p<0.05\n")
cat("  ns  not significant\n")
cat("\nFCP/FBG ratio interpretation:\n")
cat(" - Formula: FCP (ng/mL) / FBG (mg/dL)\n")
cat(" - Reflects beta-cell function relative to glycaemic state\n")
cat(" - Higher ratio: Better insulin secretion capacity\n")
cat(" - Lower ratio: Impaired beta-cell compensation\n")
cat(" - Typical range: 0.001\u20130.020 (dimensionless)\n")
cat(rep("=", 60), "\n", sep = "")
