# ============================================================
# Albumin/Globulin (A/G) Ratio Analysis
# Revised for Diabetologia Journal Submission
# Unpaired t-test throughout - uniform across panel
# One-way ANOVA as overall omnibus test
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# ============================================================

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

# 3) Define pairwise comparisons
my_comparisons <- list(
  c("CTRL", "IL10-KO"),
  c("CTRL", "DSS-IL10-KO"),
  c("IL10-KO", "DSS-IL10-KO")
)

# ============================================================
# 4) Import and calculate A/G ratio
# ============================================================
cat("\n=== LOADING BLOOD CHEMISTRY DATA ===\n")
bc_data <- read_csv("Fig_1d_Data_BC.csv")
bc_data <- bc_data %>%
  filter(Replicate %in% c("Rep1", "Rep2", "Rep3")) %>%
  mutate(Value = as.numeric(Value)) %>%
  filter(!is.na(Value))

cat("\n=== CALCULATING ALBUMIN/GLOBULIN RATIO ===\n")

ag_ratio <- bc_data %>%
  filter(Parameter %in% c("ALB", "GLOB")) %>%
  select(Group, Replicate, Parameter, Value) %>%
  pivot_wider(names_from = Parameter, values_from = Value) %>%
  filter(!is.na(ALB) & !is.na(GLOB) & GLOB != 0) %>%
  mutate(AG_Ratio = ALB / GLOB)

ag_ratio$Group <- factor(ag_ratio$Group,
                         levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))

cat("\nCalculated A/G Ratios:\n")
print(ag_ratio %>% select(Group, Replicate, ALB, GLOB, AG_Ratio))

write.csv(ag_ratio,
          file      = "Fig_1d_Data_AG_Ratio.csv",
          row.names = FALSE)
cat("✓ A/G ratio data exported\n")

cat("\nSample sizes per group:\n")
print(table(ag_ratio$Group))

# ============================================================
# 5) NORMALITY TESTING (informational only)
# ============================================================
cat("\n=== NORMALITY TESTING (Shapiro-Wilk, informational only) ===\n")

normality_results <- data.frame()

for (grp in levels(ag_ratio$Group)) {
  data_subset <- ag_ratio %>%
    filter(Group == grp) %>%
    pull(AG_Ratio)
  
  if (length(data_subset) >= 3) {
    shapiro_result <- shapiro.test(data_subset)
    
    normality_results <- rbind(normality_results, data.frame(
      Group     = grp,
      n         = length(data_subset),
      Shapiro_W = round(shapiro_result$statistic, 4),
      Shapiro_p = round(shapiro_result$p.value, 4),
      Normal    = shapiro_result$p.value > 0.05,
      Test_Used = "Unpaired t-test (uniform across panel)"
    ))
    
    cat(sprintf("%s: W = %.4f, p = %.4f (%s)\n",
                grp,
                shapiro_result$statistic,
                shapiro_result$p.value,
                ifelse(shapiro_result$p.value > 0.05,
                       "Normal", "Non-normal")))
  }
}

write.csv(normality_results,
          file      = "Fig_1d_Normality_Results.csv",
          row.names = FALSE)
cat("\n✓ Normality results exported (informational only)\n")

# ============================================================
# 6) OVERALL TEST - ONE-WAY ANOVA
# ============================================================
cat("\n=== ONE-WAY ANOVA (overall omnibus test) ===\n")

anova_result  <- aov(AG_Ratio ~ Group, data = ag_ratio)
anova_summary <- summary(anova_result)
print(anova_summary)

anova_p <- anova_summary[[1]][["Pr(>F)"]][1]
cat(sprintf("\nF = %.4f, p = %.4f\n",
            anova_summary[[1]][["F value"]][1], anova_p))

if (anova_p < 0.05) {
  cat("Overall ANOVA significant → proceed with pairwise comparisons\n")
} else {
  cat("Overall ANOVA not significant → pairwise tests reported for completeness\n")
}

# ============================================================
# 7) PAIRWISE COMPARISONS - UNPAIRED T-TEST THROUGHOUT
# ============================================================
cat("\n=== PERFORMING PAIRWISE UNPAIRED T-TESTS ===\n")

posthoc_result <- compare_means(
  formula     = AG_Ratio ~ Group,
  data        = ag_ratio,
  method      = "t.test",      # Uniform unpaired t-test throughout
  paired      = FALSE,
  comparisons = my_comparisons
) %>%
  mutate(
    p        = signif(p, 3),
    p.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    Test_Method = "unpaired t-test"
  ) %>%
  rename(Variable = .y.) %>%
  select(Variable, group1, group2, Test_Method, p, p.signif)

cat("\n=== Pairwise Results (uncorrected, unpaired t-test) ===\n")
print(posthoc_result)

write.csv(posthoc_result,
          file      = "Fig_1d_Statistical_Tests.csv",
          row.names = FALSE)
cat("\n✓ Statistical test results exported to Fig_1d_Statistical_Tests.csv\n")

# ============================================================
# 8) Summary statistics
# ============================================================
summary_stats <- ag_ratio %>%
  group_by(Group) %>%
  summarise(
    n      = n(),
    mean   = round(mean(AG_Ratio, na.rm = TRUE), 2),
    sd     = round(sd(AG_Ratio,   na.rm = TRUE), 2),
    sem    = round(sd(AG_Ratio,   na.rm = TRUE) / sqrt(n()), 2),
    median = round(median(AG_Ratio, na.rm = TRUE), 2),
    min    = round(min(AG_Ratio,  na.rm = TRUE), 2),
    max    = round(max(AG_Ratio,  na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(Mean_SEM = paste0(mean, " ± ", sem))

print(summary_stats)

write.csv(summary_stats,
          file      = "Fig_1d_Summary_AG_Ratio.csv",
          row.names = FALSE)
cat("✓ Summary statistics exported to Fig_1d_Summary_AG_Ratio.csv\n")

# ============================================================
# 9) PLOT - UNPAIRED T-TEST THROUGHOUT (consistent with CSV)
# ============================================================
cat("\n=== GENERATING BOXPLOT ===\n")

x_labels <- c("CTRL", "IL10-KO", "DSS-\nIL10-KO")
names(x_labels) <- c("CTRL", "IL10-KO", "DSS-IL10-KO")

boxplot_ag <- ggboxplot(
  ag_ratio, x = "Group", y = "AG_Ratio",
  color  = "Group", fill = "Group",
  alpha  = 0.5, add = "jitter",
  add.params = list(color = "Group", shape = 21,
                    fill  = "white", size  = 0.75),
  palette       = diabetologia_colors,
  outlier.shape = NA
) +
  lims(y = c(0, 6)) +
  labs(y = "Albumin/globulin ratio (AU)",
       x = "Experimental group") +
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
                                    linewidth = 0.75 / .pt)
  ) +
  stat_compare_means(
    comparisons   = my_comparisons,
    method        = "t.test",    # Consistent with CSV export
    paired        = FALSE,
    tip.length    = 0.0,
    step.increase = 0.3,
    label         = "p.signif", # Uncorrected p values
    size          = 2.5,
    bracket.size  = 0.3
  )

print(boxplot_ag)

ggsave(filename = "Fig_1d_Boxplot_AG_Ratio.png",
       plot = boxplot_ag,
       scale = 1, width = 5.95, height = 6.0,
       units = "cm", dpi = 600, limitsize = FALSE)

ggsave(filename = "Fig_1d_Boxplot_AG_Ratio.svg",
       plot = boxplot_ag,
       scale = 1, width = 5.95, height = 6.0,
       units = "cm", limitsize = FALSE)

cat("✓ Boxplot saved in PNG and SVG formats\n")

# ============================================================
# 10) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat(" 1. Fig_1d_Normality_Results.csv  (informational only)\n")
cat(" 2. Fig_1d_Statistical_Tests.csv  (unpaired t-test, uncorrected p)\n")
cat(" 3. Fig_1d_Boxplot_AG_Ratio.png   (600 dpi)\n")
cat(" 4. Fig_1d_Boxplot_AG_Ratio.svg   (vector format)\n")
cat(" 5. Fig_1d_Data_AG_Ratio.csv      (calculated ratios)\n")
cat(" 6. Fig_1d_Summary_AG_Ratio.csv   (descriptive statistics)\n")
cat("\nStatistical approach:\n")
cat(" - Overall test: one-way ANOVA (omnibus)\n")
cat(" - Pairwise: unpaired t-test throughout (uniform across panel)\n")
cat(" - Uncorrected p values reported\n")
cat(" - Normality tests retained for transparency (informational)\n")
cat(rep("=", 60), "\n", sep = "")
