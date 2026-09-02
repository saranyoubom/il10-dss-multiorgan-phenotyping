# ============================================================
# Large Intestine Gene Expression Analysis - Fig_1f
# FOLD CHANGE VERSION - Using facet_grid
# Integrated Script for Diabetologia Journal Submission
# Log2(Fold-change) relative to CTRL baseline
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Five genes: Mafa, Glp1r, Insr, Zfand6, Tcf19
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/05-Lab members/Ruby/IL10/Paper submission/Figures/Fig1")

# 1) Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL" = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO" = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define large intestine genes for Fig_1f
genes_fig1f <- c("Mafa", "Glp1r", "Insr", "Zfand6", "Tcf19")

# 4) Define comparisons
my_comparisons <- list(
  c("CTRL", "IL10-KO"),
  c("CTRL", "DSS-IL10-KO"),
  c("IL10-KO", "DSS-IL10-KO")
)

# ============================================================
# 5) Import and prepare data
# ============================================================
cat("\n=== LOADING LARGE INTESTINE GENE EXPRESSION DATA ===\n")
dat_raw <- read.csv("Fig_1f_Data_LargeIntestine_GeneExpr.csv", stringsAsFactors = FALSE)

# Filter for the 5 selected genes
dat_raw <- dat_raw %>%
  filter(Gene %in% genes_fig1f) %>%
  select(Group, Gene, Relative_expression, Tissue)

# Convert to numeric and factor
dat_raw$Relative_expression <- as.numeric(dat_raw$Relative_expression)
dat_raw$Group <- factor(dat_raw$Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))
dat_raw$Gene <- factor(dat_raw$Gene, levels = genes_fig1f)

cat("Total observations (raw):", nrow(dat_raw), "\n")
cat("Genes:", paste(genes_fig1f, collapse = ", "), "\n")

# ============================================================
# 6) Calculate FOLD CHANGE relative to CTRL mean
# ============================================================
cat("\n=== CALCULATING FOLD CHANGE ===\n")

# Calculate CTRL mean for each gene
ctrl_means <- dat_raw %>%
  filter(Group == "CTRL") %>%
  group_by(Gene) %>%
  summarise(
    CTRL_mean = mean(Relative_expression, na.rm = TRUE),
    CTRL_n = sum(!is.na(Relative_expression)),
    .groups = 'drop'
  )

cat("\nCTRL baseline means:\n")
print(ctrl_means)

# Calculate fold change for each sample
dat_fc <- dat_raw %>%
  left_join(ctrl_means, by = "Gene") %>%
  mutate(
    Fold_Change = Relative_expression / CTRL_mean,
    Log2_Fold_Change = log2(Fold_Change)
  ) %>%
  filter(!is.na(Log2_Fold_Change))  # Remove NA values

cat("\nTotal observations after fold change calculation:", nrow(dat_fc), "\n")

# Check sample sizes
cat("\nSample sizes per gene per group:\n")
print(table(dat_fc$Gene, dat_fc$Group))

# ============================================================
# 7) Create faceted boxplot using facet_grid
# ============================================================
cat("\n=== GENERATING FACETED PLOT ===\n")

# Create custom x-axis labels
x_labels <- c("CTRL" = "CTRL", 
              "IL10-KO" = "IL10-KO", 
              "DSS-IL10-KO" = "DSS-\nIL10-KO")

# Create the plot
plot_fig1f <- ggboxplot(
  dat_fc,
  x = "Group",
  y = "Log2_Fold_Change",
  color = "Group",
  fill = "Group",
  alpha = 0.5,
  add = "jitter",
  add.params = list(color = "Group", shape = 21,
                    fill = "white", size = 0.75, alpha = 0.8),
  palette = diabetologia_colors,
  outlier.shape = NA,
  width = 0.6
) +
  facet_grid(~ Gene, scales = "fixed")+
  scale_color_manual(values = diabetologia_colors) +
  scale_fill_manual(values = diabetologia_colors) +
  labs(
    y = expression(paste("Log"[2], "(Fold change)")),
    x = NULL
  ) +
  scale_x_discrete(labels = x_labels) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 8, face = "italic"),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 8, colour = "black", 
                               angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 8, colour = "black", 
                               angle = 90, hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 0.75 / .pt),
    axis.ticks = element_line(colour = "black", linewidth = 0.75 / .pt),
    panel.spacing = unit(0.3, "cm")
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    paired = FALSE,
    tip.length = 0.0,
    step.increase = 0.15,
    size = 2.5,
    bracket.size = 0.3
  ) +
  scale_y_continuous(
    n.breaks = 4,
    expand = expansion(mult = c(0.05, 0.15))
  )

# Display plot
print(plot_fig1f)

# ============================================================
# 8) Save plot in multiple formats
# ============================================================
cat("\n=== SAVING PLOT FILES ===\n")

ggsave(
  filename = "Fig_1f_FoldChange_LargeIntestine_GeneExpr.png",
  plot = plot_fig1f,
  scale = 1, width = 9.5, height = 6,
  units = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename = "Fig_1f_FoldChange_LargeIntestine_GeneExpr.svg",
  plot = plot_fig1f,
  scale = 1, width = 9.5, height = 6,
  units = "cm", limitsize = FALSE
)

cat("✓ Plots saved\n")

# ============================================================
# 9) Statistical analysis with BH correction
# ============================================================
cat("\n=== PERFORMING STATISTICAL ANALYSIS ===\n")

posthoc_result <- dat_fc %>%
  group_by(Gene) %>%
  wilcox_test(Log2_Fold_Change ~ Group, 
              comparisons = my_comparisons) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(
    p = signif(p, 3),
    p.adj = signif(p.adj, 3)
  )

cat("\nStatistical test results with BH-adjusted p-values:\n")
print(posthoc_result)

write.csv(
  posthoc_result,
  file = "Fig_1f_Posthoc_FoldChange_BH_corrected.csv",
  row.names = FALSE
)

cat("✓ Statistical results saved\n")

# ============================================================
# 10) Summary statistics
# ============================================================
cat("\n=== CALCULATING SUMMARY STATISTICS ===\n")

summary_stats <- dat_fc %>%
  group_by(Gene, Group) %>%
  summarise(
    n = n(),
    mean_Log2FC = round(mean(Log2_Fold_Change, na.rm = TRUE), 3),
    sd_Log2FC = round(sd(Log2_Fold_Change, na.rm = TRUE), 3),
    sem_Log2FC = round(sd_Log2FC / sqrt(n), 3),
    mean_FC = round(2^mean_Log2FC, 2),  # Convert back to fold change
    .groups = 'drop'
  )

cat("\nSummary statistics (Log2 and linear fold change):\n")
print(summary_stats)

write.csv(
  summary_stats,
  file = "Fig_1f_Summary_FoldChange_LargeIntestine_GeneExpr.csv",
  row.names = FALSE
)

# Also save the fold change data
write.csv(
  dat_fc %>% select(Group, Gene, Relative_expression, Fold_Change, Log2_Fold_Change),
  file = "Fig_1f_FoldChange_Data.csv",
  row.names = FALSE
)

cat("✓ Summary statistics saved\n")

# ============================================================
# 11) Summary
# ============================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE - LARGE INTESTINE GENES (FOLD CHANGE)\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nGenerated files:\n")
cat("  Figures:\n")
cat("    - Fig_1f_FoldChange_LargeIntestine_GeneExpr.png (600 dpi, 17.27×5.7 cm)\n")
cat("    - Fig_1f_FoldChange_LargeIntestine_GeneExpr.svg (vector format)\n")
cat("  Data tables:\n")
cat("    - Fig_1f_Posthoc_FoldChange_BH_corrected.csv\n")
cat("    - Fig_1f_Summary_FoldChange_LargeIntestine_GeneExpr.csv\n")
cat("    - Fig_1f_FoldChange_Data.csv\n")
cat("\nKey features:\n")
cat("  ✓ Log2(Fold change) relative to CTRL baseline\n")
cat("  ✓ facet_grid layout (5 panels, 1 row)\n")
cat("  ✓ Shared y-axis across all genes\n")
cat("  ✓ DSS-\\nIL10-KO label with line break\n")
cat("  ✓ Diabetologia formatting\n")
cat("  ✓ BH correction for multiple comparisons\n")
cat("  ✓ Statistical significance stars on plot\n")
cat("\nGenes analyzed:\n")
cat("  ", paste(genes_fig1f, collapse = ", "), "\n")
cat(rep("=", 70), "\n", sep = "")
