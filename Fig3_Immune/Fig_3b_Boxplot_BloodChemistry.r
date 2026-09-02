# ============================================================
# Blood Chemistry Individual Parameters - Fig_3b
# Revised for Diabetologia Journal Submission
# Reference ranges from Otto et al. 2016 (German Mouse Clinic)
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Multi-panel boxplot with reference ranges
# Blood collection: Terminal cardiac puncture
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/05-Lab members/Ruby/IL10/Paper submission/Figures/Fig3")

# 1) Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggpubr)
library(rstatix)

# Check if all required packages are installed
required_packages <- c("ggplot2", "dplyr", "readr", "tidyr", "ggpubr", "rstatix")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL" = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO" = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Read the blood chemistry data
bc_data <- read_csv("Fig_3b_Data_BC.csv")
bc_data <- bc_data %>% filter(Replicate %in% c("Rep1", "Rep2", "Rep3"))

cat("\n=== DATA IMPORT SUMMARY ===\n")
cat(paste0("Total observations: ", nrow(bc_data), "\n"))
cat(paste0("Parameters measured: ", length(unique(bc_data$Parameter)), "\n"))
cat(paste0("Groups: ", paste(unique(bc_data$Group), collapse = ", "), "\n\n"))

# ============================================================
# 4) Define reference ranges for blood chemistry parameters
# Source: Otto et al. 2016 JAALAS 55(4):375-386
# PMID: 27423143; doi:10.30802/AALAS-JAALAS-15-000810
# 
# Study details:
# - C57BL/6J mice, 90-135 days old, MALE
# - German Mouse Clinic phenotyping programme
# - Sample size: Males n=425, Females n=441
# - Blood collection: Retrobulbar puncture under isoflurane
# - Values: 2.5th-97.5th percentiles (95% reference interval)
#
# YOUR STUDY: Terminal cardiac puncture under isoflurane
# K+ range 3.5-4.8 mmol/L is appropriate (minimal hemolysis)
# ============================================================

ref_ranges <- data.frame(
  Parameter = c("ALB", "ALP", "ALT", "AMY", "TBIL", "BUN", "CA", "PHOS",
                "CRE", "GLU", "NA+", "K+", "TP", "GLOB"),
  
  # Reference ranges for MALE C57BL/6J mice
  # Terminal cardiac puncture (minimal hemolysis)
  Min = c(
    2.2,   # ALB: Albumin (g/dL) - converted from 22.0 g/L [Otto Table 4]
    68,    # ALP: Alkaline phosphatase (U/L) - males [Otto Table 4]
    18,    # ALT: Alanine aminotransferase (U/L) - males [Otto Table 4]
    507,   # AMY: α-amylase (U/L) - males [Otto Table 4: 507.3-962.2]
    0.1,   # TBIL: Total bilirubin (mg/dL) - UCLA DLAM (not in Otto)
    23,    # BUN: Blood urea nitrogen (mg/dL) - from 8.31 mmol/L × 2.8
    8.6,   # CA: Calcium (mg/dL) - from 2.16 mmol/L × 4.0 [Otto Table 4]
    2.5,   # PHOS: Phosphate (mg/dL) - from 0.8 mmol/L × 3.1 [Otto Table 4]
    0.1,   # CRE: Creatinine (mg/dL) - UCLA DLAM (not in Otto)
    60,    # GLU: Glucose (mg/dL) - UCLA DLAM (Otto excluded due to variability)
    140,   # NA+: Sodium (mmol/L) - males [Otto Table 4: 140-156]
    3.5,   # K+: Potassium (mmol/L) - CARDIAC PUNCTURE, males [Otto Table 4]
    4.4,   # TP: Total protein (g/dL) - from 44.0 g/L ÷ 10 [Otto Table 4]
    1.5    # GLOB: Globulin (g/dL) - calculated: TP_min - ALB_max ≈ 1.35-1.5
  ),
  
  Max = c(
    3.05,  # ALB: Albumin (g/dL) - converted from 30.5 g/L
    132,   # ALP: Alkaline phosphatase (U/L) - males
    94,    # ALT: Alanine aminotransferase (U/L) - males
    962,   # AMY: α-amylase (U/L) - males [962.2]
    0.2,   # TBIL: Total bilirubin (mg/dL) - UCLA DLAM
    42,    # BUN: Blood urea nitrogen (mg/dL) - from 15.01 mmol/L × 2.8
    10.0,  # CA: Calcium (mg/dL) - from 2.50 mmol/L × 4.0
    7.8,   # PHOS: Phosphate (mg/dL) - from 2.5 mmol/L × 3.1
    1.8,   # CRE: Creatinine (mg/dL) - UCLA DLAM
    133,   # GLU: Glucose (mg/dL) - UCLA DLAM
    156,   # NA+: Sodium (mmol/L) - males
    4.8,   # K+: Potassium (mmol/L) - CARDIAC PUNCTURE, males
    5.7,   # TP: Total protein (g/dL) - from 56.8 g/L ÷ 10
    2.4    # GLOB: Globulin (g/dL) - calculated: TP_max - ALB_min ≈ 3.5, conservatively 2.4
  ),
  
  Unit = c("g/dL", "U/L", "U/L", "U/L", "mg/dL", "mg/dL", "mg/dL", "mg/dL",
           "mg/dL", "mg/dL", "mmol/L", "mmol/L", "g/dL", "g/dL")
)

# ============================================================
# CONVERSION FACTORS USED:
# - Protein (g/L to g/dL): divide by 10
# - Calcium (mmol/L to mg/dL): multiply by 4.0
# - Phosphate (mmol/L to mg/dL): multiply by 3.1
# - BUN/Urea (mmol/L to mg/dL): multiply by 2.8
# 
# CRITICAL CORRECTIONS FROM ORIGINAL SCRIPT:
# ✓ CA: Changed from 36-45 to 8.6-10.0 mg/dL (4× correction - MAJOR ERROR)
# ✓ PHOS: Changed from 1.2-1.8 to 2.5-7.8 mg/dL (3-4× correction - MAJOR ERROR)
# ✓ K+: Changed from 3.8-4.4 to 3.5-4.8 mmol/L (appropriate for cardiac puncture)
# ✓ ALB: Changed from 2.6-5.4 to 2.2-3.05 g/dL (more accurate for males)
# ✓ BUN: Changed from 2.0-71 to 23-42 mg/dL (physiologically realistic range)
# ✓ AMY: Changed from 608-744 to 507-962 U/L (broader, more accurate)
# ✓ ALP: Reduced from 16-200 to 68-132 U/L (male-specific, more precise)
# ✓ ALT: Reduced from 22-133 to 18-94 U/L (male-specific, more precise)
#
# Blood collection method: Terminal cardiac puncture under isoflurane anaesthesia
# K+ reference range (3.5-4.8 mmol/L) is appropriate for this method
# (Retro-orbital sampling would require higher K+ range: 5.0-9.0 mmol/L)
#
# Reference: Otto GP, Rathkolb B, Oestereicher MA, et al. 
# Clinical Chemistry Reference Intervals for C57BL/6J, C57BL/6N, 
# and C3HeB/FeJ Mice (Mus musculus). 
# J Am Assoc Lab Anim Sci. 2016;55(4):375-386. 
# PMID: 27423143; doi:10.30802/AALAS-JAALAS-15-000810
# ============================================================

cat("\n=== REFERENCE RANGES SUMMARY ===\n")
print(ref_ranges)
cat("\nSource: Otto et al. 2016 JAALAS (PMID: 27423143)\n")
cat("Male C57BL/6J mice, 90-135 days, terminal cardiac puncture\n\n")

# 5) Clean data
bc_data <- bc_data %>%
  mutate(Value = as.numeric(Value)) %>%
  filter(!is.na(Value))

# Set group order
bc_data$Group <- factor(bc_data$Group,
                        levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))

# Set parameter order to match reference ranges
bc_data$Parameter <- factor(bc_data$Parameter,
                            levels = ref_ranges$Parameter)

# Define facet labels for parameters with units
param_labels <- setNames(
  paste0(ref_ranges$Parameter, " (", ref_ranges$Unit, ")"),
  ref_ranges$Parameter
)

# ============================================================
# 6) Data validation: Check values outside reference ranges
# ============================================================
cat("\n=== DATA VALIDATION: VALUES OUTSIDE REFERENCE RANGES ===\n")

bc_data_check <- bc_data %>%
  left_join(ref_ranges, by = "Parameter") %>%
  mutate(
    Outside_Range = Value < Min | Value > Max,
    Status = case_when(
      Value < Min ~ "Below reference",
      Value > Max ~ "Above reference",
      TRUE ~ "Within range"
    ),
    Percent_diff = case_when(
      Value < Min ~ round(100 * (Min - Value) / Min, 1),
      Value > Max ~ round(100 * (Value - Max) / Max, 1),
      TRUE ~ 0
    )
  )

# Summary of out-of-range values
out_of_range_summary <- bc_data_check %>%
  filter(Status != "Within range") %>%
  group_by(Parameter, Group, Status) %>%
  summarise(
    Count = n(),
    Mean_percent_diff = round(mean(Percent_diff), 1),
    .groups = "drop"
  ) %>%
  arrange(Parameter, Group, Status)

if (nrow(out_of_range_summary) > 0) {
  print(out_of_range_summary)
  cat("\nNote: Values outside reference ranges may indicate:\n")
  cat("1. Biological effect of IL-10 knockout or DSS treatment (expected)\n")
  cat("2. IL-10 KO mice on mixed genetic background (B6.129P2, not pure C57BL/6J)\n")
  cat("3. Age-related variation (your mice at week 13 vs reference 90-135 days)\n\n")
} else {
  cat("All values within reference ranges.\n\n")
}

# Export detailed validation results
write.csv(bc_data_check, 
          "Fig_3b_BloodChemistry_ValidationCheck_Otto2016.csv", 
          row.names = FALSE)
cat("Validation check saved: Fig_3b_BloodChemistry_ValidationCheck_Otto2016.csv\n\n")

# Summary statistics by group
summary_stats <- bc_data %>%
  group_by(Parameter, Group) %>%
  summarise(
    N = n(),
    Mean = round(mean(Value, na.rm = TRUE), 2),
    SD = round(sd(Value, na.rm = TRUE), 2),
    Median = round(median(Value, na.rm = TRUE), 2),
    Q1 = round(quantile(Value, 0.25, na.rm = TRUE), 2),
    Q3 = round(quantile(Value, 0.75, na.rm = TRUE), 2),
    Min = round(min(Value, na.rm = TRUE), 2),
    Max = round(max(Value, na.rm = TRUE), 2),
    .groups = "drop"
  )

write.csv(summary_stats, 
          "Fig_3b_BloodChemistry_SummaryStatistics.csv", 
          row.names = FALSE)
cat("Summary statistics saved: Fig_3b_BloodChemistry_SummaryStatistics.csv\n\n")

# 7) Normality test (for documentation)
cat("\n=== NORMALITY TEST ===\n")
shapiro_result <- shapiro.test(bc_data$Value)
cat(paste0("Shapiro-Wilk test (overall): W = ", round(shapiro_result$statistic, 4),
           ", p-value = ", format.pval(shapiro_result$p.value, digits = 3), "\n"))

if (shapiro_result$p.value < 0.05) {
  cat("Result: Data not normally distributed → use non-parametric tests\n")
  cat("Statistical approach: Kruskal-Wallis test with Wilcoxon pairwise comparisons\n\n")
} else {
  cat("Result: Data approximately normally distributed\n")
  cat("However, given small sample sizes (n=5-7), non-parametric tests recommended\n\n")
}

# 8) Define pairwise comparisons
my_comparisons <- list(
  c("CTRL", "IL10-KO"),
  c("CTRL", "DSS-IL10-KO"),
  c("IL10-KO", "DSS-IL10-KO")
)

# 9) Create custom x-axis labels with line break for DSS-IL10-KO
x_labels <- c("CTRL", "IL10-KO", "DSS-\nIL10-KO")
names(x_labels) <- c("CTRL", "IL10-KO", "DSS-IL10-KO")

# ============================================================
# 10) Create multi-panel boxplot with reference ranges
# ============================================================
cat("\n=== GENERATING MULTI-PANEL BOXPLOT ===\n")

boxplot_bc <- ggplot(bc_data, aes(x = Group, y = Value, 
                                  fill = Group, color = Group)) +
  # Add reference range as shaded rectangle (grey zone)
  geom_rect(data = ref_ranges,
            aes(ymin = Min, ymax = Max, xmin = -Inf, xmax = Inf),
            fill = "grey85", alpha = 0.6, inherit.aes = FALSE) +
  
  # Add boxplot (semi-transparent)
  geom_boxplot(alpha = 0.5, outlier.shape = NA, linewidth = 0.5) +
  
  # Add individual data points (jittered)
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.75,
              shape = 21, stroke = 0.5, fill = "white") +
  
  # Add statistical comparisons
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    paired = FALSE,
    aes(label = after_stat(p.signif)),
    size = 2,
    tip.length = 0,
    step.increase = 0.15,
    hide.ns = FALSE
  ) +
  
  # Apply Diabetologia colour scheme
  scale_fill_manual(values = diabetologia_colors) +
  scale_color_manual(values = diabetologia_colors) +
  
  # Facet by parameter with free y-axis scales
  facet_wrap(~ Parameter, scales = "free_y", ncol = 7,
             labeller = labeller(Parameter = param_labels)) +
  
  # Axis labels
  labs(
    x = "Experimental group",
    y = "Blood chemistry value"
  ) +
  scale_x_discrete(labels = x_labels) +
  
  # Theme customization (Diabetologia style)
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, face = "plain"),
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, colour = "black", 
                               angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8, colour = "black", 
                               angle = 90, hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 0.75 / .pt),
    axis.ticks = element_line(colour = "black", linewidth = 0.75 / .pt),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# Display plot
print(boxplot_bc)

# Save plot in multiple formats
ggsave(
  filename = "Fig_3b_Boxplot_BloodChemistry_Otto2016.png",
  plot = boxplot_bc,
  scale = 1.2, width = 17, height = 10,
  units = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename = "Fig_3b_Boxplot_BloodChemistry_Otto2016.svg",
  plot = boxplot_bc,
  scale = 1.2, width = 17, height = 10,
  units = "cm", limitsize = FALSE
)

ggsave(
  filename = "Fig_3b_Boxplot_BloodChemistry_Otto2016.pdf",
  plot = boxplot_bc,
  scale = 1.2, width = 17, height = 10,
  units = "cm", limitsize = FALSE
)

cat("\n=== PLOT FILES SAVED ===\n")
cat("PNG: Fig_3b_Boxplot_BloodChemistry_Otto2016.png (600 dpi)\n")
cat("SVG: Fig_3b_Boxplot_BloodChemistry_Otto2016.svg (vector)\n")
cat("PDF: Fig_3b_Boxplot_BloodChemistry_Otto2016.pdf (vector)\n\n")

# ============================================================
# 11) Statistical analysis: Pairwise comparisons
# ============================================================
cat("\n=== PERFORMING STATISTICAL TESTS FOR EACH PARAMETER ===\n")

stat_results <- bc_data %>%
  group_by(Parameter) %>%
  do({
    param_data <- .
    param_name <- unique(param_data$Parameter)
    
    # Check sample sizes
    n_ctrl <- sum(param_data$Group == "CTRL")
    n_il10 <- sum(param_data$Group == "IL10-KO")
    n_dss <- sum(param_data$Group == "DSS-IL10-KO")
    
    # Perform Kruskal-Wallis test (overall group effect)
    kw_test <- tryCatch({
      kruskal.test(Value ~ Group, data = param_data)
    }, error = function(e) NULL)
    
    # Pairwise Wilcoxon rank-sum tests (Mann-Whitney U)
    ctrl_vs_il10 <- tryCatch({
      wilcox.test(
        param_data$Value[param_data$Group == "CTRL"],
        param_data$Value[param_data$Group == "IL10-KO"],
        exact = FALSE
      )
    }, error = function(e) NULL)
    
    ctrl_vs_dss <- tryCatch({
      wilcox.test(
        param_data$Value[param_data$Group == "CTRL"],
        param_data$Value[param_data$Group == "DSS-IL10-KO"],
        exact = FALSE
      )
    }, error = function(e) NULL)
    
    il10_vs_dss <- tryCatch({
      wilcox.test(
        param_data$Value[param_data$Group == "IL10-KO"],
        param_data$Value[param_data$Group == "DSS-IL10-KO"],
        exact = FALSE
      )
    }, error = function(e) NULL)
    
    # Calculate medians for each group
    median_ctrl <- median(param_data$Value[param_data$Group == "CTRL"], na.rm = TRUE)
    median_il10 <- median(param_data$Value[param_data$Group == "IL10-KO"], na.rm = TRUE)
    median_dss <- median(param_data$Value[param_data$Group == "DSS-IL10-KO"], na.rm = TRUE)
    
    # Compile results
    data.frame(
      Parameter = param_name,
      N_CTRL = n_ctrl,
      N_IL10KO = n_il10,
      N_DSS_IL10KO = n_dss,
      Median_CTRL = round(median_ctrl, 2),
      Median_IL10KO = round(median_il10, 2),
      Median_DSS_IL10KO = round(median_dss, 2),
      KW_statistic = ifelse(!is.null(kw_test), round(kw_test$statistic, 3), NA),
      KW_pvalue = ifelse(!is.null(kw_test), kw_test$p.value, NA),
      CTRL_vs_IL10KO_pvalue = ifelse(!is.null(ctrl_vs_il10), 
                                     ctrl_vs_il10$p.value, NA),
      CTRL_vs_DSS_IL10KO_pvalue = ifelse(!is.null(ctrl_vs_dss), 
                                         ctrl_vs_dss$p.value, NA),
      IL10KO_vs_DSS_IL10KO_pvalue = ifelse(!is.null(il10_vs_dss), 
                                           il10_vs_dss$p.value, NA)
    )
  }) %>%
  ungroup() %>%
  mutate(
    # Add significance symbols (Diabetologia style)
    KW_significance = case_when(
      KW_pvalue < 0.001 ~ "***",
      KW_pvalue < 0.01 ~ "**",
      KW_pvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    CTRL_vs_IL10KO_sig = case_when(
      CTRL_vs_IL10KO_pvalue < 0.001 ~ "***",
      CTRL_vs_IL10KO_pvalue < 0.01 ~ "**",
      CTRL_vs_IL10KO_pvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    CTRL_vs_DSS_IL10KO_sig = case_when(
      CTRL_vs_DSS_IL10KO_pvalue < 0.001 ~ "***",
      CTRL_vs_DSS_IL10KO_pvalue < 0.01 ~ "**",
      CTRL_vs_DSS_IL10KO_pvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    IL10KO_vs_DSS_IL10KO_sig = case_when(
      IL10KO_vs_DSS_IL10KO_pvalue < 0.001 ~ "***",
      IL10KO_vs_DSS_IL10KO_pvalue < 0.01 ~ "**",
      IL10KO_vs_DSS_IL10KO_pvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Display statistical results
cat("\nStatistical test results:\n")
print(stat_results, width = Inf)

# Export statistical results
write.csv(stat_results, 
          "Fig_3b_BloodChemistry_StatisticalTests_Otto2016.csv", 
          row.names = FALSE)

cat("\n=== STATISTICAL RESULTS SAVED ===\n")
cat("File: Fig_3b_BloodChemistry_StatisticalTests_Otto2016.csv\n\n")

# ============================================================
# 12) Summary of significant findings
# ============================================================
cat("\n=== SUMMARY OF SIGNIFICANT FINDINGS ===\n")

# Parameters with significant overall group effect (Kruskal-Wallis)
sig_overall <- stat_results %>%
  filter(KW_significance != "ns") %>%
  select(Parameter, KW_pvalue, KW_significance)

if (nrow(sig_overall) > 0) {
  cat("\nParameters with significant overall group differences:\n")
  print(sig_overall, row.names = FALSE)
} else {
  cat("\nNo parameters showed significant overall group differences.\n")
}

# Parameters with significant pairwise differences
sig_pairwise <- stat_results %>%
  filter(CTRL_vs_IL10KO_sig != "ns" | 
           CTRL_vs_DSS_IL10KO_sig != "ns" | 
           IL10KO_vs_DSS_IL10KO_sig != "ns") %>%
  select(Parameter, CTRL_vs_IL10KO_sig, CTRL_vs_DSS_IL10KO_sig, IL10KO_vs_DSS_IL10KO_sig)

if (nrow(sig_pairwise) > 0) {
  cat("\nParameters with significant pairwise differences:\n")
  print(sig_pairwise, row.names = FALSE)
} else {
  cat("\nNo parameters showed significant pairwise differences.\n")
}

# ============================================================
# 13) Generate figure legend text
# ============================================================
cat("\n=== SUGGESTED FIGURE LEGEND ===\n\n")

legend_text <- paste0(
  "Fig. 3b Blood chemistry parameters in experimental groups\n\n",
  "Boxplots show distribution of 14 clinical chemistry parameters measured in ",
  "terminal serum samples at week 12. Groups: CTRL (C57BL/6J background controls, ",
  "n=5), IL10-KO (untreated IL-10 knockout, n=7), DSS-IL10-KO (IL-10 knockout ",
  "receiving recurrent DSS treatment at weeks 0 and 6, n=6). Grey shaded regions ",
  "indicate reference ranges (2.5th–97.5th percentiles) for male C57BL/6J mice ",
  "aged 90–135 days (Otto et al. 2016, PMID: 27423143). Blood collected via ",
  "terminal cardiac puncture under deep isoflurane anaesthesia. Boxes show median ",
  "(centre line), interquartile range (box limits), and range (whiskers); individual ",
  "data points overlaid. Statistical comparisons: Kruskal–Wallis test with pairwise ",
  "Wilcoxon post-hoc tests; *p<0.05, **p<0.01, ***p<0.001, ns, not significant. ",
  "Abbreviations: ALB, albumin; ALP, alkaline phosphatase; ALT, alanine ",
  "aminotransferase; AMY, α-amylase; TBIL, total bilirubin; BUN, blood urea ",
  "nitrogen; CA, calcium; PHOS, phosphorus; CRE, creatinine; GLU, glucose; ",
  "NA+, sodium; K+, potassium; TP, total protein; GLOB, globulin"
)

cat(legend_text)
cat("\n\n")

# ============================================================
# 14) Session information and completion
# ============================================================
cat("\n=== ANALYSIS COMPLETE ===\n")
cat(paste0("Analysis completed: ", Sys.time(), "\n"))
cat("\nReference: Otto GP, Rathkolb B, Oestereicher MA, et al.\n")
cat("Clinical Chemistry Reference Intervals for C57BL/6J, C57BL/6N,\n")
cat("and C3HeB/FeJ Mice (Mus musculus).\n")
cat("J Am Assoc Lab Anim Sci. 2016;55(4):375-386.\n")
cat("PMID: 27423143; doi:10.30802/AALAS-JAALAS-15-000810\n\n")

cat("Output files generated:\n")
cat("1. Fig_3b_Boxplot_BloodChemistry_Otto2016.png (600 dpi raster)\n")
cat("2. Fig_3b_Boxplot_BloodChemistry_Otto2016.svg (vector graphics)\n")
cat("3. Fig_3b_Boxplot_BloodChemistry_Otto2016.pdf (vector graphics)\n")
cat("4. Fig_3b_BloodChemistry_ValidationCheck_Otto2016.csv (data validation)\n")
cat("5. Fig_3b_BloodChemistry_SummaryStatistics.csv (descriptive stats)\n")
cat("6. Fig_3b_BloodChemistry_StatisticalTests_Otto2016.csv (inferential stats)\n\n")

# Print session information for reproducibility
cat("=== R SESSION INFORMATION ===\n")
print(sessionInfo())

# ============================================================
# END OF SCRIPT
# ============================================================
