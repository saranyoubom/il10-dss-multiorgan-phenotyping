# ============================================================
# Serum Biochemistry – Supplementary Figure S1 (consolidated)
# Diabetologia Journal Submission
#
# Reference ranges: Otto et al. 2016 JAALAS 55(4):375–386
#   PMID: 27423143 | male C57BL/6J | 90–135 days
#   2.5th–97.5th percentile reference intervals
#
# Groups    : CTRL, IL10-KO, DSS-IL10-KO (n=3; Rep1–Rep3)
# Parameters: ALT, ALB, GLOB, PHOS, BUN, CA, GLU, AMY, NA+, TP
#   Layout  : 2 rows × 5 columns
#   Row 1   : ALT | ALB | GLOB | PHOS | BUN
#   Row 2   : CA  | GLU | AMY  | NA+  | TP
#
# EXCLUDED parameters (systematic iQC failures — not reportable):
#   K+   — haemolysis artefact (iQC 0808/HEM) across all groups
#   CRE  — below detection limit (iQC 0001/0801)
#   ALP  — artefactual zero in all IL10-KO replicates (iQC 0804)
#   TBIL — substrate interference in all IL10-KO replicates (iQC 0801)
#
# Special data handling (not annotated in plot):
#   IL10-KO ALT Rep2 — EXCLUDED (iQC 0810 haemolysis)
#     → IL10-KO ALT n=2 (Rep1=317, Rep3=837)
#   CTRL    ALT Rep3 = 640  — minimum value (≥640 U/L;  iQC 0003)
#   IL10-KO ALT Rep3 = 837  — minimum value (≥837 U/L;  iQC 0003)
#   CTRL    GLOB Rep3 = 1.5 — derived (TP − ALB)
#   IL10-KO GLOB Rep1 = 1.4 — derived (TP − ALB)
#   IL10-KO GLOB Rep2 = 1.8 — derived (TP − ALB)
#   NA+: iQC 0080 (approximate) in CTRL Rep3, IL10-KO Rep1 & Rep3
#   All above retained in analysis; no special shape annotation.
#
# Statistical test: Welch two-sample t-test (var.equal = FALSE)
#   BH correction applied within each parameter (k=3 tests)
#   Significant findings (per-parameter BH):
#     PHOS CTRL vs IL10-KO           p_adj = 0.009  **
#     PHOS IL10-KO vs DSS-IL10-KO    p_adj = 0.009  **
#     GLOB CTRL vs DSS-IL10-KO       p_adj = 0.022  *
#     ALB  CTRL vs DSS-IL10-KO       p_adj = 0.029  *
#
# Bug fixes applied:
#   (1) stat_results and stat_table_plot are SEPARATE objects;
#       write.csv uses stat_results BEFORE bracket positioning
#   (2) add_xy_position() REPLACED by manual y-position:
#       y.position = y_max × (1 + step_size × rank)
#       rank: 1=CTRL vs IL10-KO, 2=CTRL vs DSS-IL10-KO,
#             3=IL10-KO vs DSS-IL10-KO
#   (3) fill = Group and colour = Group removed from global aes();
#       scoped locally to geom_boxplot() and geom_jitter() only —
#       prevents stat_pvalue_manual() inheriting Group aesthetics
#       absent from stat_table_plot → eliminates:
#         Warning: Duplicated aesthetics after name standardisation
#         Error:   object 'Group' not found
#   (4) inherit.aes = FALSE added to stat_pvalue_manual() and
#       geom_rect() as belt-and-braces guard
#   (5) ggsave() uses ifelse() for dpi — NULL not accepted
#   (6) K+, CRE, ALP, TBIL excluded before all analysis
#   (7) IL10-KO Rep2 ALT set to NA before analysis
#   (8) point_type annotation removed — all points as shape 21
#   (9) ALL three pairwise brackets drawn on every panel;
#       label = "p.adj.signif" shows ns / * / ** / ***
#       Increase step_size to 0.25 if bracket text overlaps
#       strip labels
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/05-Lab members/Ruby/IL10/Paper submission/Figures/Supplementary")

# 1) Load required libraries
required_packages <- c("ggplot2", "dplyr", "readr", "tidyr", "ggpubr", "rstatix")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggpubr)
library(rstatix)

# 2) Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Read blood chemistry data
bc_raw <- read_csv("Fig_S1_Data_BloodChemistry.csv")
bc_raw <- bc_raw %>% filter(Replicate %in% c("Rep1", "Rep2", "Rep3"))

cat("\n=== DATA IMPORT SUMMARY ===\n")
cat(paste0("Raw observations (Rep1–Rep3): ", nrow(bc_raw), "\n"))

# ============================================================
# 4) Reference ranges — 10 retained parameters
# Source: Otto et al. 2016 JAALAS 55(4):375–386
#   PMID: 27423143 | male C57BL/6J | 90–135 days
#   Reference interval: 2.5th–97.5th percentile (n=425 males)
# ============================================================

param_order <- c("ALT", "ALB", "GLOB", "PHOS", "BUN",
                 "CA",  "GLU", "AMY",  "NA+",  "TP")

ref_ranges <- data.frame(
  Parameter = param_order,
  Min = c(18,   2.2,  1.5, 2.5, 23,  8.6,  60,  507, 140, 4.4),
  Max = c(94,   3.05, 2.4, 7.8, 42, 10.0, 133,  962, 156, 5.7),
  Unit = c("U/L", "g/dL", "g/dL", "mg/dL", "mg/dL",
           "mg/dL", "mg/dL", "U/L", "mmol/L", "g/dL")
)

cat("\n=== REFERENCE RANGES SUMMARY ===\n")
print(ref_ranges)
cat("\nSource: Otto et al. 2016 JAALAS (PMID: 27423143)\n")
cat("Male C57BL/6J mice, 90–135 days, retrobulbar puncture under isoflurane\n\n")

# ============================================================
# 5) Data preparation
#   - Restrict to 10 retained parameters
#   - Exclude IL10-KO Rep2 ALT (iQC 0810 haemolysis suppression)
#   - All points plotted uniformly (shape 21, no annotation)
# ============================================================

bc_s1 <- bc_raw %>%
  mutate(Value = as.numeric(Value)) %>%
  filter(Parameter %in% param_order) %>%
  mutate(
    Value = if_else(
      Parameter == "ALT" & Group == "IL10-KO" & Replicate == "Rep2",
      NA_real_, Value
    )
  ) %>%
  filter(!is.na(Value))

bc_s1$Group     <- factor(bc_s1$Group,     levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))
bc_s1$Parameter <- factor(bc_s1$Parameter, levels = param_order)

cat("=== DATA PREPARATION SUMMARY ===\n")
cat("Excluded from all analysis: K+, CRE, ALP, TBIL (iQC failures)\n")
cat("Excluded IL10-KO Rep2 ALT (iQC 0810) → IL10-KO ALT n=2\n")
cat("ALT minimum values (CTRL Rep3, IL10-KO Rep3) retained as-is\n")
cat("GLOB derived values (CTRL Rep3, IL10-KO Rep1/2) retained as-is\n")
cat("NA+ approximate values (CTRL Rep3, IL10-KO Rep1/3) retained as-is\n\n")

n_summary <- bc_s1 %>%
  group_by(Parameter, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n)
cat("n per parameter per group:\n")
print(n_summary, n = Inf)
cat("\n")

# Facet labels with units
param_labels <- setNames(
  paste0(ref_ranges$Parameter, " (", ref_ranges$Unit, ")"),
  ref_ranges$Parameter
)

# ============================================================
# 6) Data validation and descriptive statistics
# ============================================================

bc_check <- bc_s1 %>%
  left_join(ref_ranges, by = "Parameter") %>%
  mutate(
    Status = case_when(
      Value < Min ~ "Below reference",
      Value > Max ~ "Above reference",
      TRUE        ~ "Within range"
    ),
    Percent_diff = case_when(
      Value < Min ~ round(100 * (Min - Value) / Min, 1),
      Value > Max ~ round(100 * (Value - Max) / Max, 1),
      TRUE        ~ 0
    )
  )

out_of_range <- bc_check %>%
  filter(Status != "Within range") %>%
  group_by(Parameter, Group, Status) %>%
  summarise(Count = n(),
            Mean_percent_diff = round(mean(Percent_diff), 1),
            .groups = "drop") %>%
  arrange(Parameter, Group, Status)

cat("\n=== DATA VALIDATION: VALUES OUTSIDE REFERENCE RANGES ===\n")
if (nrow(out_of_range) > 0) {
  print(out_of_range, n = Inf)
  cat("\nNote: Values outside reference ranges may reflect:\n")
  cat("  1. Biological effect of IL-10 deficiency or DSS treatment\n")
  cat("  2. Mixed genetic background (B6.129P2 vs pure C57BL/6J)\n")
  cat("  3. ALB: VetScan BCG method may overestimate vs human standard\n")
  cat("  4. PHOS: Otto range may underestimate aged VetScan CDP values\n")
  cat("  5. GLU: terminal non-fasted samples expected above fasting reference\n\n")
} else {
  cat("All values within reference ranges.\n\n")
}

write.csv(bc_check,
          "Fig_S1_BloodChemistry_ValidationCheck_Otto2016.csv",
          row.names = FALSE)
cat("Validation check saved: Fig_S1_BloodChemistry_ValidationCheck_Otto2016.csv\n\n")

summary_stats <- bc_s1 %>%
  group_by(Parameter, Group) %>%
  summarise(
    N      = n(),
    Mean   = round(mean(Value,           na.rm = TRUE), 2),
    SD     = round(sd(Value,             na.rm = TRUE), 2),
    Median = round(median(Value,         na.rm = TRUE), 2),
    Q1     = round(quantile(Value, 0.25, na.rm = TRUE), 2),
    Q3     = round(quantile(Value, 0.75, na.rm = TRUE), 2),
    Min    = round(min(Value,            na.rm = TRUE), 2),
    Max    = round(max(Value,            na.rm = TRUE), 2),
    .groups = "drop"
  )

write.csv(summary_stats,
          "Fig_S1_BloodChemistry_SummaryStatistics.csv",
          row.names = FALSE)
cat("Summary statistics saved: Fig_S1_BloodChemistry_SummaryStatistics.csv\n\n")

# ============================================================
# 7) Statistical analysis: Welch t-test + BH correction
#
# OBJECT 1: stat_results    — pure t-test output; CSV export
# OBJECT 2: stat_table_plot — manual y-positions; plot only
#
# Manual y-position formula (replaces add_xy_position()):
#   y.position = y_max × (1 + step_size × rank)
#   rank 1 = CTRL vs IL10-KO
#   rank 2 = CTRL vs DSS-IL10-KO
#   rank 3 = IL10-KO vs DSS-IL10-KO
#   Increase step_size to 0.25 if bracket text overlaps strip labels.
#
# Bracket label convention:
#   label = "p.adj.signif" → ns / * / ** / *** on every bracket
#   All three pre-specified comparisons drawn on every panel.
# ============================================================

cat("\n=== STATISTICAL ANALYSIS: Welch t-test + BH correction ===\n")
cat("Method: t.test(var.equal = FALSE) via rstatix::t_test()\n")
cat("BH correction: within each parameter (k=3 comparisons)\n\n")

# OBJECT 1 ─────────────────────────────────────────────────────
stat_results <- bc_s1 %>%
  group_by(Parameter) %>%
  t_test(
    Value ~ Group,
    var.equal       = FALSE,
    p.adjust.method = "BH",
    detailed        = FALSE
  ) %>%
  add_significance(
    "p.adj",
    cutpoints = c(0,     0.001, 0.01, 0.05, 1),
    symbols   = c("***", "**",  "*",  "ns")
  )

cat("Welch t-test results (BH-adjusted p within each parameter):\n")
print(
  stat_results %>%
    select(Parameter, group1, group2, n1, n2,
           statistic, df, p, p.adj, p.adj.signif),
  n = Inf, width = Inf
)

write.csv(
  stat_results %>%
    select(Parameter, group1, group2, n1, n2,
           statistic, df, p, p.adj, p.adj.signif),
  "Fig_S1_BloodChemistry_WelchTtest_BH.csv",
  row.names = FALSE
)
cat("\nStatistical results saved: Fig_S1_BloodChemistry_WelchTtest_BH.csv\n\n")

# OBJECT 2 ─────────────────────────────────────────────────────
step_size <- 0.20

y_maxima <- bc_s1 %>%
  group_by(Parameter) %>%
  summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")

# ALL comparisons passed to plot (no significance filter)
# label = "p.adj.signif" renders ns / * / ** / *** on every bracket
stat_table_plot <- stat_results %>%
  left_join(y_maxima, by = "Parameter") %>%
  group_by(Parameter) %>%
  mutate(
    row_num    = row_number(),
    y.position = y_max * (1 + step_size * row_num)
  ) %>%
  ungroup()

cat(paste0("Total bracket annotations (all comparisons): ",
           nrow(stat_table_plot), "\n"))
cat(paste0("  Significant (labelled):      ",
           sum(stat_table_plot$p.adj.signif != "ns"), "\n"))
cat(paste0("  Not significant (ns label):  ",
           sum(stat_table_plot$p.adj.signif == "ns"), "\n\n"))

cat("Significant comparisons:\n")
print(
  stat_table_plot %>%
    filter(p.adj.signif != "ns") %>%
    select(Parameter, group1, group2, p.adj, p.adj.signif, y.position),
  n = Inf, width = Inf
)
cat("\n")

# ============================================================
# 8) Multi-panel boxplot — 10 parameters, 2 rows × 5 columns
#
# KEY FIX — aes() scoping:
#   fill and colour NOT in global aes(); scoped locally to
#   geom_boxplot() and geom_jitter() only.
#   stat_pvalue_manual() has inherit.aes = FALSE.
#   Eliminates "Duplicated aesthetics" warning and
#   "object 'Group' not found" error.
# ============================================================

cat("\n=== GENERATING MULTI-PANEL BOXPLOT ===\n")

boxplot_s1 <- ggplot(bc_s1,
                     aes(x = Group, y = Value)) +   # fill/colour NOT global
  
  # Otto 2016 reference range (grey band)
  geom_rect(
    data        = ref_ranges,
    aes(ymin = Min, ymax = Max, xmin = -Inf, xmax = Inf),
    fill        = "grey85",
    alpha       = 0.6,
    inherit.aes = FALSE
  ) +
  
  # Box-and-whisker — fill and colour scoped locally
  geom_boxplot(
    aes(fill = Group, colour = Group),
    alpha         = 0.5,
    outlier.shape = NA,
    linewidth     = 0.5
  ) +
  
  # Individual data points — colour scoped locally; fill fixed white
  geom_jitter(
    aes(colour = Group),
    width  = 0.2,
    alpha  = 0.75,
    size   = 0.75,
    shape  = 21,
    stroke = 0.5,
    fill   = "white"
  ) +
  
  scale_fill_manual(values = diabetologia_colors, guide = "none") +
  scale_colour_manual(values = diabetologia_colors, guide = "none") +
  
  facet_wrap(
    ~ Parameter,
    scales   = "free_y",
    ncol     = 5,
    labeller = labeller(Parameter = param_labels)
  ) +
  
  scale_x_discrete(
    labels = c("CTRL", "IL10-KO", "DSS-\nIL10-KO")
  ) +
  
  labs(x = NULL, y = "Concentration in serum") +
  
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 8, face = "plain"),
    strip.background = element_blank(),
    panel.spacing    = unit(0.3, "cm"),
    axis.title.y     = element_text(size = 8),
    axis.text.x      = element_text(size = 8, colour = "black",
                                    angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y      = element_text(size = 8, colour = "black", angle= 90, 
                                    vjust = 0.5, hjust = 0.5),
    axis.line        = element_line(colour = "black",
                                    linewidth = 0.75 / .pt),
    axis.ticks       = element_line(colour = "black",
                                    linewidth = 0.75 / .pt),
    plot.margin      = margin(5, 5, 5, 5, "mm")
  ) +
  
  # All three pairwise brackets on every panel
  # label = "p.adj.signif": ns / * / ** / *** as appropriate
  stat_pvalue_manual(
    data        = stat_table_plot,
    label       = "p.adj.signif",
    xmin        = "group1",
    xmax        = "group2",
    y.position  = "y.position",
    facet.by    = "Parameter",
    tip.length  = 0.0,
    size        = 2.5,
    colour      = "black",
    inherit.aes = FALSE
  )

print(boxplot_s1)

for (ext in c("png", "svg", "pdf")) {
  ggsave(
    filename  = paste0("Fig_S1_Boxplot_BloodChemistry.", ext),
    plot      = boxplot_s1,
    scale     = 1.0,
    width     = 15.5,
    height    = 12,
    units     = "cm",
    dpi       = ifelse(ext == "png", 600, 300),
    limitsize = FALSE
  )
}

cat("\n=== PLOT FILES SAVED ===\n")
cat("PNG: Fig_S1_Boxplot_BloodChemistry.png (600 dpi)\n")
cat("SVG: Fig_S1_Boxplot_BloodChemistry.svg (vector)\n")
cat("PDF: Fig_S1_Boxplot_BloodChemistry.pdf (vector)\n\n")

# ============================================================
# 9) Summary of significant findings
# ============================================================

cat("\n=== SUMMARY OF SIGNIFICANT FINDINGS ===\n")

sig_summary <- stat_results %>%
  filter(p.adj.signif != "ns") %>%
  select(Parameter, group1, group2, n1, n2,
         p, p.adj, p.adj.signif)

if (nrow(sig_summary) > 0) {
  cat("\nSignificant pairwise comparisons (BH-adjusted within each parameter):\n")
  print(sig_summary, row.names = FALSE, width = Inf)
} else {
  cat("\nNo pairwise comparisons reached significance after BH correction.\n")
}

# ============================================================
# 10) Session information
# ============================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(paste0("Completed: ", Sys.time(), "\n"))
cat("\nReference: Otto GP et al. J Am Assoc Lab Anim Sci. 2016;55(4):375-386.\n")
cat("PMID: 27423143; doi:10.30802/AALAS-JAALAS-15-000810\n\n")

cat("Output files:\n")
cat("1. Fig_S1_Boxplot_BloodChemistry.png              (600 dpi)\n")
cat("2. Fig_S1_Boxplot_BloodChemistry.svg              (vector)\n")
cat("3. Fig_S1_Boxplot_BloodChemistry.pdf              (vector)\n")
cat("4. Fig_S1_BloodChemistry_ValidationCheck_Otto2016.csv\n")
cat("5. Fig_S1_BloodChemistry_SummaryStatistics.csv\n")
cat("6. Fig_S1_BloodChemistry_WelchTtest_BH.csv\n\n")

print(sessionInfo())

# ============================================================
# END OF SCRIPT
# ============================================================
