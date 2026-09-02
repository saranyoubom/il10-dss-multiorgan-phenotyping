# ============================================================
# Pancreatic β-Cell Gene Expression Analysis - Fig 2e
# SEPARATED PLOTS with Fig_1c theme - COMBINED EXPORT ONLY
# Revised for Diabetologia Journal Submission
# Log2(Fold-change) relative to CTRL baseline
# Genes: Mafa, Ins1, Ins2 (Glp1r excluded - all CTRL ND)
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

# 3) Define β-cell genes and all pairwise comparisons
beta_cell_genes <- c("Mafa", "Ins1", "Ins2")

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
  filter(Gene %in% beta_cell_genes) %>%
  filter(!is.na(Log2_Fold_Change))

dat$Group <- factor(dat$Group, levels = z_order)
dat$Gene  <- factor(dat$Gene,  levels = beta_cell_genes)

cat("\n=== DATA STRUCTURE FOR β-CELL GENES ===\n")
print(head(dat, 15))
cat("\nTotal observations:", nrow(dat), "\n")
cat("\nSample sizes per gene per group:\n")
print(table(dat$Gene, dat$Group))

# Real data for statistics (before dummy rows are added)
dat_real <- dat %>% filter(!is.na(Log2_Fold_Change))

# ============================================================
# 5) Add dummy NA rows for Ins1 IL10-KO if ND
#    (ensures IL10-KO appears on x-axis without plotted points)
# ============================================================
ins1_il10ko_count <- sum(dat_real$Gene == "Ins1" & dat_real$Group == "IL10-KO")

if (ins1_il10ko_count == 0) {
  cat("\n\u26a0 Ins1 IL10-KO is ND \u2014 adding NA placeholder for x-axis display\n")
  dummy_rows <- data.frame(
    Gene             = factor("Ins1", levels = beta_cell_genes),
    Group            = factor("IL10-KO", levels = z_order),
    Log2_Fold_Change = NA_real_,
    Sample           = paste0("IL10-KO_dummy_", 1:4)
  )
  dat <- bind_rows(dat, dummy_rows)
  dat$Gene  <- factor(dat$Gene,  levels = beta_cell_genes)
  dat$Group <- factor(dat$Group, levels = z_order)
}

# ============================================================
# 6) Calculate shared y-axis limits
# ============================================================
y_min    <- min(dat_real$Log2_Fold_Change, na.rm = TRUE)
y_max    <- max(dat_real$Log2_Fold_Change, na.rm = TRUE)
y_buffer <- (y_max - y_min) * 0.15
y_limits <- c(y_min - y_buffer, (y_max + y_buffer) * 1.7)
y_range  <- diff(y_limits)
cat("\nShared Y-axis limits:", round(y_limits, 2), "\n")

# ============================================================
# 7) STATISTICAL ANALYSIS — run BEFORE plot
#    Wilcoxon rank-sum test; BH correction within each gene
#    droplevels() prevents wilcox_test error on empty factor levels
#    Same results used for both plot labels and CSV export
# ============================================================
cat("\n=== PERFORMING STATISTICAL ANALYSIS ===\n")
cat("Method: Wilcoxon rank-sum test (independent samples)\n")
cat("Correction: Benjamini-Hochberg applied within each gene\n\n")

posthoc_result_list <- list()

for (gene in beta_cell_genes) {
  dat_gene <- dat_real %>%
    filter(Gene == gene) %>%
    mutate(Group = droplevels(Group))   # Remove empty factor levels
  
  il10ko_n <- sum(dat_gene$Group == "IL10-KO")
  
  comparisons_gene <- if (gene == "Ins1" && il10ko_n == 0) {
    cat("  Ins1: Using only CTRL vs DSS-IL10-KO (IL10-KO ND)\n")
    list(c("CTRL", "DSS-IL10-KO"))
  } else {
    my_comparisons
  }
  
  if (nrow(dat_gene) > 0) {
    posthoc_gene <- dat_gene %>%
      wilcox_test(Log2_Fold_Change ~ Group,
                  comparisons = comparisons_gene) %>%
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
  file      = "Fig_2e_Posthoc_Beta_Cell_Genes.csv",
  row.names = FALSE
)
cat("\u2713 Statistical results exported to Fig_2e_Posthoc_Beta_Cell_Genes.csv\n")


# ============================================================
# 8) GENERATE INDIVIDUAL PLOTS
#    Significance brackets use p.adj.signif from posthoc_result
#    via stat_pvalue_manual() — consistent with CSV export
#    Bracket order matches my_comparisons via factor-level sort
# ============================================================
cat("\n=== GENERATING INDIVIDUAL PLOTS ===\n")

plot_list  <- list()
gene_index <- 0

for (gene_name in beta_cell_genes) {
  gene_index <- gene_index + 1
  cat("Creating plot for:", gene_name, "\n")
  
  dat_gene      <- dat %>% filter(Gene == gene_name)
  il10ko_count  <- sum(!is.na(dat_gene$Log2_Fold_Change[dat_gene$Group == "IL10-KO"]))
  show_nd_label <- (gene_name == "Ins1" && il10ko_count == 0)
  show_y_axis   <- (gene_index == 1)
  
  y_title <- if (show_y_axis) expression(paste("Log"[2], "(Fold change)")) else NULL
  
  # Retrieve BH-corrected stats for this gene
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
  
  # Y-axis: show only on first panel (Mafa)
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
  
  # ND annotation for Ins1 IL10-KO
  if (show_nd_label) {
    p <- p + annotate(
      "text",
      x        = 2,
      y        = y_limits[1],
      label    = "ND",
      color    = "grey60",
      size     = 2.5,
      fontface = "bold"
    )
  }
  
  plot_list[[gene_name]] <- p
}

# ============================================================
# 9) Combine and export plots (combined only)
# ============================================================
cat("\n=== COMBINING AND EXPORTING PLOTS ===\n")

combined_plot <- wrap_plots(plot_list, nrow = 1, ncol = 3)

ggsave(
  filename  = "Fig_2e_Beta_Cell_Genes.png",
  plot      = combined_plot,
  scale = 1, width = 9.7, height = 5.7,
  units = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename  = "Fig_2e_Beta_Cell_Genes.svg",
  plot      = combined_plot,
  scale = 1, width = 9.7, height = 5.7,
  units = "cm", limitsize = FALSE
)

cat("\u2713 Combined plot saved in PNG and SVG formats\n")
print(combined_plot)

# ============================================================
# 10) Summary statistics
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
  file      = "Fig_2e_Summary_Beta_Cell_Genes.csv",
  row.names = FALSE
)
cat("\u2713 Summary statistics exported to Fig_2e_Summary_Beta_Cell_Genes.csv\n")

# ============================================================
# 11) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE \u2014 \u03b2-CELL GENES\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat("  Figures:\n")
cat("    - Fig_2e_Beta_Cell_Genes.png (600 dpi, 12.5\u00d75.7 cm)\n")
cat("    - Fig_2e_Beta_Cell_Genes.svg (vector format)\n")
cat("  Data tables:\n")
cat("    - Fig_2e_Posthoc_Beta_Cell_Genes.csv (Wilcoxon, BH-corrected within gene)\n")
cat("    - Fig_2e_Summary_Beta_Cell_Genes.csv\n")
cat("\nKey features:\n")
cat("  \u2713 Statistics run BEFORE plot (plot labels match CSV exactly)\n")
cat("  \u2713 droplevels() prevents wilcox_test error on empty Ins1 IL10-KO level\n")
cat("  \u2713 Bracket order matches my_comparisons via z_order factor-level sort\n")
cat("  \u2713 stat_pvalue_manual() with p.adj.signif (BH-corrected)\n")
cat("  \u2713 Same y-axis limits across all genes\n")
cat("  \u2713 Y-axis title and axis only on Mafa (first panel)\n")
cat("  \u2713 Y-axis removed from Ins1 and Ins2 (panels 2 and 3)\n")
cat("  \u2713 IL10-KO shown on x-axis for Ins1 with ND label\n")
cat("  \u2713 Fig_1c theme applied consistently\n")
cat("\nNote:\n")
cat("  \u2022 Glp1r excluded (all CTRL samples ND)\n")
cat("  \u2022 Ins1 IL10-KO undetectable (0/4) \u2014 shown with ND label\n")
cat("  \u2022 Statistics for Ins1: only CTRL vs DSS-IL10-KO\n")
cat("\nStatistical approach:\n")
cat("  \u2022 Test: Wilcoxon rank-sum (independent samples, between-group)\n")
cat("  \u2022 Correction: Benjamini-Hochberg within each gene\n")
cat("  \u2022 Plot labels: p.adj.signif (BH-corrected, consistent with CSV)\n")
cat("  \u2022 Significance: *** p<0.001, ** p<0.01, * p<0.05, ns\n")
cat("  \u2022 Bracket order: CTRL/IL10-KO \u2192 CTRL/DSS-IL10-KO \u2192 IL10-KO/DSS-IL10-KO\n")
cat(rep("=", 60), "\n", sep = "")
