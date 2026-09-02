# ============================================================
# Pancreatic Stress Response Gene Expression Analysis - Fig 2g
# INDIVIDUAL PLOTS with Fig_1c theme - SEPARATE EXPORTS
# Revised for Diabetologia Journal Submission
# Log2(Fold-change) relative to CTRL baseline
# Genes: Zfand6, Srr
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Y-AXIS: Shown on both individual plots
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# 1) Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr)
library(rstatix)
library(tidyr)
library(patchwork)

# 2) Define Diabetologia colour palette (matching Fig_1c)
diabetologia_colors <- c(
  "CTRL" = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO" = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define stress response genes
stress_genes <- c("Zfand6", "Srr")

# 4) Define comparisons
my_comparisons <- list(
  c("CTRL", "IL10-KO"),
  c("CTRL", "DSS-IL10-KO"),
  c("IL10-KO", "DSS-IL10-KO")
)

# ============================================================
# 5) Import pancreas gene expression data
# ============================================================
cat("\n=== LOADING PANCREAS GENE EXPRESSION DATA ===\n")
dat <- read_csv("Fig_2ef_Pancreas_FoldChange_Data.csv")

# Filter for stress response genes only
dat <- dat %>%
  filter(Gene %in% stress_genes) %>%
  filter(!is.na(Log2_Fold_Change))

# Convert variables to factors
dat$Group <- as.factor(dat$Group)
dat$Gene <- as.factor(dat$Gene)

# Set the order of groups and genes
z_order <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
dat$Group <- factor(dat$Group, levels = z_order)
dat$Gene <- factor(dat$Gene, levels = stress_genes)

# Check data structure
cat("\n=== DATA STRUCTURE FOR STRESS GENES ===\n")
print(head(dat, 15))
cat("\nTotal observations:", nrow(dat), "\n")
cat("\nSample sizes per gene per group:\n")
print(table(dat$Gene, dat$Group))

# ============================================================
# 6) Calculate individual y-axis limits for each gene
# ============================================================
cat("\n=== CALCULATING Y-AXIS LIMITS ===\n")

y_limits_list <- list()

for (gene_name in stress_genes) {
  dat_gene <- dat %>% filter(Gene == gene_name)

  if (nrow(dat_gene) > 0) {
    y_min <- min(dat_gene$Log2_Fold_Change, na.rm = TRUE)
    y_max <- max(dat_gene$Log2_Fold_Change, na.rm = TRUE)
    y_buffer <- (y_max - y_min) * 0.15

    y_limits_list[[gene_name]] <- c(y_min - y_buffer, (y_max + y_buffer)*1.8)

    cat("  ", gene_name, ": [", round(y_limits_list[[gene_name]][1], 2), 
        ", ", round(y_limits_list[[gene_name]][2], 2), "]\n", sep = "")
  }
}

# ============================================================
# 7) Create and export INDIVIDUAL plots for each gene
# ============================================================
cat("\n=== GENERATING AND EXPORTING INDIVIDUAL PLOTS ===\n")

plot_list <- list()

for (gene_name in stress_genes) {
  cat("\nCreating plot for:", gene_name, "\n")

  # Filter data for this gene
  dat_gene <- dat %>% filter(Gene == gene_name)

  if (nrow(dat_gene) == 0) {
    cat("  WARNING: No data for", gene_name, "- skipping\n")
    next
  }

  # Use all comparisons
  comparisons_to_use <- my_comparisons

  # Y-axis title - shown on all individual plots
  y_title <- expression(paste("Log"[2], "(Fold change)"))

  # Get gene-specific y-axis limits
  y_limits <- y_limits_list[[gene_name]]

  # Create individual plot with Fig_1c theme
  p <- ggboxplot(
    dat_gene,
    x = "Group",
    y = "Log2_Fold_Change",
    color = "Group",
    fill = "Group",
    alpha = 0.5,
    add = "jitter",
    add.params = list(color = "Group", shape = 21,
                      fill = "white", size = 0.75),
    palette = diabetologia_colors,
    outlier.shape = NA,
    width = 0.6
  ) +
    labs(
      y = y_title,
      x = NULL,
      title = bquote(italic(.(gene_name)))
    ) +
    # Set gene-specific y-axis limits
    ylim(y_limits[1], y_limits[2]) +
    # Custom x-axis labels with line break for DSS-IL10-KO
    scale_x_discrete(labels = c("CTRL" = "CTRL", 
                                "IL10-KO" = "IL10-KO", 
                                "DSS-IL10-KO" = "DSS-\nIL10-KO")) +
    # EXACT THEME FROM FIG_1C
    theme_classic() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      strip.background = element_blank(),
      axis.title = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8, colour = "black",
                                 angle = 0, hjust = 0.5, vjust = 1),
      axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.75 / .pt),
      axis.ticks = element_line(colour = "black", linewidth = 0.75 / .pt),
      plot.title = element_text(size = 8, face = "italic", hjust = 0.5),
      plot.margin = margin(3.8, 3.8, 3.8, 3.8, "mm")
    )

  # Add statistical comparisons (Fig_1c style)
  if (length(comparisons_to_use) > 0) {
    p <- p + stat_compare_means(
      comparisons = comparisons_to_use,
      method = "wilcox.test",
      paired = FALSE,
      tip.length = 0.0,
      step.increase = 0.2,
      aes(label = after_stat(p.signif)),
      size = 2
    )
  }

  # Store plot in list
  plot_list[[gene_name]] <- p

  # Export individual plot
  filename_png <- paste0("Fig_2g_", gene_name, ".png")
  filename_svg <- paste0("Fig_2g_", gene_name, ".svg")

  ggsave(
    filename = filename_png,
    plot = p,
    scale = 1, width = 4.7, height = 5.7,
    units = "cm", dpi = 600, limitsize = FALSE
  )

  ggsave(
    filename = filename_svg,
    plot = p,
    scale = 1, width = 4.7, height = 5.7,
    units = "cm", limitsize = FALSE
  )

  cat("  ✓ Saved:", filename_png, "\n")
  cat("  ✓ Saved:", filename_svg, "\n")

  # Display plot
  print(p)
}

# ============================================================
# 8) Statistical analysis
# ============================================================
cat("\n=== PERFORMING STATISTICAL ANALYSIS ===\n")

# Use all data
dat_real <- dat %>% filter(!is.na(Log2_Fold_Change))

posthoc_result_list <- list()

for (gene in stress_genes) {
  dat_gene <- dat_real %>% filter(Gene == gene)

  if (nrow(dat_gene) > 0) {
    posthoc_gene <- dat_gene %>%
      wilcox_test(Log2_Fold_Change ~ Group, comparisons = my_comparisons) %>%
      mutate(Gene = gene)

    posthoc_result_list[[gene]] <- posthoc_gene
  }
}

posthoc_result <- bind_rows(posthoc_result_list) %>%
  group_by(Gene) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(
    p = signif(p, 3),
    p.adj = signif(p.adj, 3)
  )

print(posthoc_result)

write.csv(
  posthoc_result,
  file = "Fig_2g_Posthoc_Stress_Genes.csv",
  row.names = FALSE
)

# ============================================================
# 9) Summary statistics
# ============================================================
cat("\n=== CALCULATING SUMMARY STATISTICS ===\n")

summary_stats <- dat_real %>%
  group_by(Gene, Group) %>%
  summarise(
    n = n(),
    mean_Log2FC = round(mean(Log2_Fold_Change, na.rm = TRUE), 3),
    sd_Log2FC = round(sd(Log2_Fold_Change, na.rm = TRUE), 3),
    sem_Log2FC = round(sd_Log2FC / sqrt(n), 3),
    .groups = 'drop'
  )

print(summary_stats)

write.csv(
  summary_stats,
  file = "Fig_2g_Summary_Stress_Genes.csv",
  row.names = FALSE
)

# ============================================================
# 10) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE - STRESS RESPONSE GENES\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat("  Individual Figures:\n")
cat("    - Fig_2g_Zfand6.png (600 dpi, 4.2×5.7 cm)\n")
cat("    - Fig_2g_Zfand6.svg (vector format)\n")
cat("    - Fig_2g_Srr.png (600 dpi, 4.2×5.7 cm)\n")
cat("    - Fig_2g_Srr.svg (vector format)\n")
cat("  Data tables:\n")
cat("    - Fig_2g_Posthoc_Stress_Genes.csv\n")
cat("    - Fig_2g_Summary_Stress_Genes.csv\n")
cat("\nKey features:\n")
cat("  ✓ Individual y-axis limits for each gene\n")
cat("  ✓ Y-axis title and axis shown on both plots\n")
cat("  ✓ Separate PNG and SVG files for each gene\n")
cat("  ✓ DSS-\\nIL10-KO label with line break\n")
cat("  ✓ Fig_1c theme applied consistently\n")
cat(rep("=", 60), "\n", sep = "")
