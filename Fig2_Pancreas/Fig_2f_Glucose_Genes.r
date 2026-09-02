# ============================================================
# Pancreatic Glucose Metabolism Gene Expression Analysis - Fig 2f
# SEPARATED PLOTS with Fig_1c theme - COMBINED EXPORT ONLY
# Revised for Diabetologia Journal Submission
# Log2(Fold-change) relative to CTRL baseline
# Genes: Insr, Tbc1d4, Glut1
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Statistical method: Wilcoxon rank-sum, BH correction within gene
# Plot labels (p.adj.signif) consistent with exported CSV
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/05-Lab members/Ruby/IL10/Paper submission/Figures/Fig2")

# 1) Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr)
library(rstatix)
library(tidyr)
library(patchwork)

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define glucose metabolism genes and all pairwise comparisons
glucose_genes <- c("Insr", "Tbc1d4", "Glut1")

my_comparisons <- list(
  c("CTRL", "IL10-KO"),
  c("CTRL", "DSS-IL10-KO"),
  c("IL10-KO", "DSS-IL10-KO")
)

z_order <- c("CTRL", "IL10-KO", "DSS-IL10-KO")

# ============================================================
# 4) Import pancreas gene expression data
# ============================================================
cat("\n=== LOADING PANCREAS GENE EXPRESSION DATA ===\n")
dat <- read_csv("Fig_2ef_Pancreas_FoldChange_Data.csv")

dat <- dat %>%
  filter(Gene %in% glucose_genes) %>%
  filter(!is.na(Log2_Fold_Change))

dat$Group <- factor(dat$Group, levels = z_order)
dat$Gene  <- factor(dat$Gene,  levels = glucose_genes)

cat("\n=== DATA STRUCTURE FOR GLUCOSE GENES ===\n")
print(head(dat, 15))
cat("\nTotal observations:", nrow(dat), "\n")
cat("\nSample sizes per gene per group:\n")
print(table(dat$Gene, dat$Group))

# Real data for statistics
dat_real <- dat %>% filter(!is.na(Log2_Fold_Change))

# ============================================================
# 5) Calculate shared y-axis limits
# ============================================================
y_min    <- min(dat_real$Log2_Fold_Change, na.rm = TRUE)
y_max    <- max(dat_real$Log2_Fold_Change, na.rm = TRUE)
y_buffer <- (y_max - y_min) * 0.15
y_limits <- c(y_min - y_buffer, (y_max + y_buffer) * 2.0)
y_range  <- diff(y_limits)
cat("\nShared Y-axis limits:", round(y_limits, 2), "\n")

# ============================================================
# 6) STATISTICAL ANALYSIS — run BEFORE plot
#    Wilcoxon rank-sum test; BH correction within each gene
#    droplevels() guards against empty factor level errors
#    Same results used for both plot labels and CSV export
# ============================================================
cat("\n=== PERFORMING STATISTICAL ANALYSIS ===\n")
cat("Method: Wilcoxon rank-sum test (independent samples)\n")
cat("Correction: Benjamini-Hochberg applied within each gene\n\n")

posthoc_result_list <- list()

for (gene in glucose_genes) {
  dat_gene <- dat_real %>%
    filter(Gene == gene) %>%
    mutate(Group = droplevels(Group))   # Guard against empty factor levels
  
  if (nrow(dat_gene) > 0) {
    posthoc_gene <- dat_gene %>%
      wilcox_test(Log2_Fold_Change ~ Group,
                  comparisons = my_comparisons) %>%
      mutate(Gene = gene)
    posthoc_result_list[[gene]] <- posthoc_gene
  }
}

# BH correction within each gene, then add significance labels
posthoc_result <- bind_rows(posthoc_result_list) %>%
  group_by(Gene) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(
    p            = signif(p, 3),
    p.adj        = signif(p.adj, 3),
    Test_Method  = "Wilcoxon rank-sum",      # ← added
    p_adj_method = "Benjamini-Hochberg"      # ← added
  ) %>%
  ungroup()

cat("\n=== Statistical Results (Wilcoxon, BH-corrected within gene) ===\n")
print(posthoc_result %>%
        select(Gene, group1, group2, n1, n2, statistic,
               Test_Method, p, p_adj_method, p.adj, p.adj.signif))

write.csv(
  posthoc_result %>%
    select(Gene, group1, group2, n1, n2, statistic,
           Test_Method, p, p_adj_method, p.adj, p.adj.signif),   # ← updated
  file      = "Fig_2f_Posthoc_Glucose_Genes.csv",
  row.names = FALSE
)
cat("\u2713 Statistical results exported to Fig_2f_Posthoc_Glucose_Genes.csv\n")


# ============================================================
# 7) GENERATE INDIVIDUAL PLOTS
#    Significance brackets use p.adj.signif from posthoc_result
#    via stat_pvalue_manual() — consistent with CSV export
#    Bracket order matches my_comparisons via factor-level sort
# ============================================================
cat("\n=== GENERATING INDIVIDUAL PLOTS ===\n")

plot_list  <- list()
gene_index <- 0

for (gene_name in glucose_genes) {
  gene_index <- gene_index + 1
  cat("Creating plot for:", gene_name, "\n")
  
  dat_gene    <- dat %>% filter(Gene == gene_name)
  show_y_axis <- (gene_index == 1)
  
  if (nrow(dat_gene) == 0) {
    cat("  WARNING: No data for", gene_name, "— skipping\n")
    next
  }
  
  y_title <- if (show_y_axis) expression(paste("Log"[2], "(Fold change)")) else NULL
  
  # Retrieve BH-corrected stats for this gene and assign bracket y positions
  gene_stats <- posthoc_result %>% filter(Gene == gene_name)
  n_brackets <- nrow(gene_stats)
  
  if (n_brackets > 0) {
    gene_data_max  <- max(dat_real$Log2_Fold_Change[dat_real$Gene == gene_name],
                          na.rm = TRUE)
    bracket_start  <- gene_data_max + y_range * 0.08
    bracket_height <- y_range * 0.14
    
    # Sort by z_order factor levels (not alphabetically) to match my_comparisons
    gene_stats <- gene_stats %>%
      mutate(
        group1 = factor(group1, levels = z_order),
        group2 = factor(group2, levels = z_order)
      ) %>%
      arrange(group1, group2) %>%
      mutate(
        group1     = as.character(group1),   # Restore to character for stat_pvalue_manual
        group2     = as.character(group2),
        y.position = bracket_start + (seq_len(n()) - 1) * bracket_height
      )
  }
  
  # Build boxplot
  p <- ggboxplot(
    dat_gene,
    x = "Group", y = "Log2_Fold_Change",
    color = "Group", fill = "Group",
    alpha = 0.5, add = "jitter",
    add.params = list(color = "Group", shape = 21,
                      fill = "white", size = 0.75),
    palette       = diabetologia_colors,
    outlier.shape = NA,
    width         = 0.6
  ) +
    labs(y = y_title, x = NULL,
         title = bquote(italic(.(gene_name)))) +
    ylim(y_limits[1], y_limits[2]) +
    scale_x_discrete(labels = c("CTRL"        = "CTRL",
                                "IL10-KO"     = "IL10-KO",
                                "DSS-IL10-KO" = "DSS-\nIL10-KO")) +
    theme_classic() +
    theme(
      legend.position  = "none",
      strip.background = element_blank(),
      axis.title       = element_text(size = 8),
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 8, colour = "black",
                                      angle = 90, hjust = 1, vjust = 0.5),
      axis.line        = element_line(colour = "black", linewidth = 0.75 / .pt),
      axis.ticks       = element_line(colour = "black", linewidth = 0.75 / .pt),
      plot.title       = element_text(size = 8, face = "italic", hjust = 0.5),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
  
  # Y-axis: show only on first panel (Insr)
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank()
    )
  } else {
    p <- p + theme(
      axis.text.y = element_text(size = 8, colour = "black",
                                 angle = 90, hjust = 0.5)
    )
  }
  
  # Add BH-corrected significance brackets (consistent with CSV)
  if (n_brackets > 0) {
    p <- p + stat_pvalue_manual(
      gene_stats,
      label        = "p.adj.signif",   # BH-corrected — matches CSV
      tip.length   = 0,
      size         = 2,
      bracket.size = 0.3
    )
  }
  
  plot_list[[gene_name]] <- p
}

# ============================================================
# 8) Combine and export plots (combined only)
# ============================================================
cat("\n=== COMBINING AND EXPORTING PLOTS ===\n")

combined_plot <- wrap_plots(plot_list, nrow = 1, ncol = 3)

ggsave(
  filename  = "Fig_2f_Glucose_Genes.png",
  plot      = combined_plot,
  scale = 1, width = 9.7, height = 5.7,
  units = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename  = "Fig_2f_Glucose_Genes.svg",
  plot      = combined_plot,
  scale = 1, width = 9.7, height = 5.7,
  units = "cm", limitsize = FALSE
)

cat("\u2713 Combined plot saved in PNG and SVG formats\n")
print(combined_plot)

# ============================================================
# 9) Summary statistics
# ============================================================
cat("\n=== CALCULATING SUMMARY STATISTICS ===\n")

summary_stats <- dat_real %>%
  group_by(Gene, Group) %>%
  summarise(
    n           = n(),
    mean_Log2FC = round(mean(Log2_Fold_Change, na.rm = TRUE), 3),
    sd_Log2FC   = round(sd(Log2_Fold_Change,   na.rm = TRUE), 3),
    sem_Log2FC  = round(sd_Log2FC / sqrt(n), 3),
    .groups = 'drop'
  )

print(summary_stats)

write.csv(
  summary_stats,
  file      = "Fig_2f_Summary_Glucose_Genes.csv",
  row.names = FALSE
)
cat("\u2713 Summary statistics exported to Fig_2f_Summary_Glucose_Genes.csv\n")

# ============================================================
# 10) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE \u2014 GLUCOSE METABOLISM GENES\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat("  Figures:\n")
cat("    - Fig_2f_Glucose_Genes.png (600 dpi, 12.5\u00d75.7 cm)\n")
cat("    - Fig_2f_Glucose_Genes.svg (vector format)\n")
cat("  Data tables:\n")
cat("    - Fig_2f_Posthoc_Glucose_Genes.csv (Wilcoxon, BH-corrected within gene)\n")
cat("    - Fig_2f_Summary_Glucose_Genes.csv\n")
cat("\nKey features:\n")
cat("  \u2713 Statistics run BEFORE plot (plot labels match CSV exactly)\n")
cat("  \u2713 droplevels() guards against empty factor level errors\n")
cat("  \u2713 Bracket order matches my_comparisons via z_order factor-level sort\n")
cat("  \u2713 stat_pvalue_manual() with p.adj.signif (BH-corrected)\n")
cat("  \u2713 Same y-axis limits across all genes\n")
cat("  \u2713 Y-axis title and axis only on Insr (first panel)\n")
cat("  \u2713 Y-axis removed from Tbc1d4 and Glut1 (panels 2 and 3)\n")
cat("  \u2713 Fig_1c theme applied consistently\n")
cat("\nStatistical approach:\n")
cat("  \u2022 Test: Wilcoxon rank-sum (independent samples, between-group)\n")
cat("  \u2022 Correction: Benjamini-Hochberg within each gene\n")
cat("  \u2022 Plot labels: p.adj.signif (BH-corrected, consistent with CSV)\n")
cat("  \u2022 Significance: *** p<0.001, ** p<0.01, * p<0.05, ns\n")
cat("  \u2022 Bracket order: CTRL/IL10-KO \u2192 CTRL/DSS-IL10-KO \u2192 IL10-KO/DSS-IL10-KO\n")
cat(rep("=", 60), "\n", sep = "")
