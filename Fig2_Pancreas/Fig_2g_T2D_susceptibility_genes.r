# ============================================================
# Pancreatic T2D Susceptibility Gene Expression Analysis - Fig 2g
# COMBINED PLOT (2 columns) with Fig_1c theme
# Revised for Diabetologia Journal Submission
# Log2(Fold-change) relative to CTRL baseline
# Genes: Zfand6, Srr (T2D susceptibility genes)
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Statistical method: Wilcoxon rank-sum, BH correction within gene
# Plot labels (p.adj.signif) consistent with exported CSV
# Y-AXIS: Shared scale across both panels
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

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define T2D susceptibility genes and all pairwise comparisons
stress_genes <- c("Zfand6", "Srr")

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
  filter(Gene %in% stress_genes) %>%
  filter(!is.na(Log2_Fold_Change))

dat$Group <- factor(dat$Group, levels = z_order)
dat$Gene  <- factor(dat$Gene,  levels = stress_genes)

cat("\n=== DATA STRUCTURE FOR T2D SUSCEPTIBILITY GENES ===\n")
print(head(dat, 15))
cat("\nTotal observations:", nrow(dat), "\n")
cat("\nSample sizes per gene per group:\n")
print(table(dat$Gene, dat$Group))

# Real data for statistics
dat_real <- dat %>% filter(!is.na(Log2_Fold_Change))

# ============================================================
# 5) Calculate SHARED y-axis limits across both genes
# ============================================================
cat("\n=== CALCULATING SHARED Y-AXIS LIMITS ===\n")

y_min_global <- min(dat_real$Log2_Fold_Change, na.rm = TRUE)
y_max_global <- max(dat_real$Log2_Fold_Change, na.rm = TRUE)
y_buffer     <- (y_max_global - y_min_global) * 0.15
y_shared     <- c(y_min_global - y_buffer, (y_max_global + y_buffer) * 2.2)

cat("  Shared y-axis limits: [",
    round(y_shared[1], 2), ", ",
    round(y_shared[2], 2), "]\n", sep = "")

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

for (gene in stress_genes) {
  dat_gene <- dat_real %>%
    filter(Gene == gene) %>%
    mutate(Group = droplevels(Group))
  
  if (nrow(dat_gene) > 0) {
    posthoc_gene <- dat_gene %>%
      wilcox_test(Log2_Fold_Change ~ Group,
                  comparisons = my_comparisons) %>%
      mutate(Gene = gene)
    posthoc_result_list[[gene]] <- posthoc_gene
  }
}

posthoc_result <- bind_rows(posthoc_result_list) %>%
  group_by(Gene) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(
    p            = signif(p, 3),
    p.adj        = signif(p.adj, 3),
    Test_Method  = "Wilcoxon rank-sum",
    p_adj_method = "Benjamini-Hochberg"
  ) %>%
  ungroup()

cat("\n=== Statistical Results (Wilcoxon, BH-corrected within gene) ===\n")
print(posthoc_result %>%
        select(Gene, group1, group2, n1, n2, statistic,
               Test_Method, p, p_adj_method, p.adj, p.adj.signif))

write.csv(
  posthoc_result %>%
    select(Gene, group1, group2, n1, n2, statistic,
           Test_Method, p, p_adj_method, p.adj, p.adj.signif),
  file      = "Fig_2g_T2D_susceptibility_genes_Posthoc.csv",
  row.names = FALSE
)
cat("✓ Statistical results exported to Fig_2g_T2D_susceptibility_genes_Posthoc.csv\n")

# ============================================================
# 7) GENERATE INDIVIDUAL PLOTS (stored in plot_list)
#    Significance brackets use p.adj.signif from posthoc_result
#    via stat_pvalue_manual() — consistent with CSV export
#    Shared y-axis scale across both panels
# ============================================================
cat("\n=== GENERATING PLOTS ===\n")

plot_list <- list()

y_range <- diff(y_shared)  # Derived from shared scale

for (gene_name in stress_genes) {
  cat("\nCreating plot for:", gene_name, "\n")
  
  dat_gene <- dat %>% filter(Gene == gene_name)
  
  if (nrow(dat_gene) == 0) {
    cat("  WARNING: No data for", gene_name, "— skipping\n")
    next
  }
  
  y_title <- expression(paste("Log"[2], "(Fold change)"))
  
  gene_stats <- posthoc_result %>% filter(Gene == gene_name)
  n_brackets <- nrow(gene_stats)
  
  if (n_brackets > 0) {
    bracket_start  <- y_max_global + y_range * 0.08   # Anchored to global max
    bracket_height <- y_range * 0.14
    
    gene_stats <- gene_stats %>%
      mutate(
        group1 = factor(group1, levels = z_order),
        group2 = factor(group2, levels = z_order)
      ) %>%
      arrange(group1, group2) %>%
      mutate(
        group1     = as.character(group1),
        group2     = as.character(group2),
        y.position = bracket_start + (seq_len(n()) - 1) * bracket_height
      )
  }
  
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
    ylim(y_shared[1], y_shared[2]) +                  # Shared scale
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
      axis.text.y      = element_text(size = 8, colour = "black",
                                      angle = 90, hjust = 0.5),
      axis.line        = element_line(colour = "black", linewidth = 0.75 / .pt),
      axis.ticks       = element_line(colour = "black", linewidth = 0.75 / .pt),
      plot.title       = element_text(size = 8, face = "italic", hjust = 0.5),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
  
  if (n_brackets > 0) {
    p <- p + stat_pvalue_manual(
      gene_stats,
      label        = "p.adj.signif",
      tip.length   = 0,
      size         = 2,
      bracket.size = 0.3
    )
  }
  
  plot_list[[gene_name]] <- p
  print(p)
}

# ============================================================
# 7b) COMBINED FIGURE — 2 columns, 7.2 × 5.7 cm
#     Y-axis title retained on left panel (Zfand6) only
# ============================================================
cat("\n=== GENERATING COMBINED FIGURE ===\n")

# Step 7b — replace the plot_combined line with:
plot_combined <- plot_list[["Zfand6"]] +
  (plot_list[["Srr"]] + theme(
    axis.title.y = element_blank(),   # suppress title
    axis.text.y  = element_blank(),   # suppress tick labels
    axis.line.y  = element_blank(),   # suppress axis line
    axis.ticks.y = element_blank()    # suppress ticks
  )) +
  plot_layout(ncol = 2)


ggsave(
  filename = "Fig_2g_T2D_susceptibility_genes_Combined.png",
  plot     = plot_combined,
  scale    = 1, width = 7.52, height = 5.7,
  units    = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename = "Fig_2g_T2D_susceptibility_genes_Combined.svg",
  plot     = plot_combined,
  scale    = 1, width = 7.52, height = 5.7,
  units    = "cm", limitsize = FALSE
)

print(plot_combined)
cat("✓ Saved: Fig_2g_T2D_susceptibility_genes_Combined.png (600 dpi, 7.2×5.7 cm)\n")
cat("✓ Saved: Fig_2g_T2D_susceptibility_genes_Combined.svg (vector format)\n")

# ============================================================
# 8) Summary statistics
# ============================================================
cat("\n=== CALCULATING SUMMARY STATISTICS ===\n")

summary_stats <- dat_real %>%
  group_by(Gene, Group) %>%
  summarise(
    n           = n(),
    mean_Log2FC = round(mean(Log2_Fold_Change, na.rm = TRUE), 3),
    sd_Log2FC   = round(sd(Log2_Fold_Change,   na.rm = TRUE), 3),
    sem_Log2FC  = round(sd_Log2FC / sqrt(n), 3),
    .groups     = "drop"
  )

print(summary_stats)

write.csv(
  summary_stats,
  file      = "Fig_2g_T2D_susceptibility_genes_Summary.csv",
  row.names = FALSE
)
cat("✓ Summary statistics exported to Fig_2g_T2D_susceptibility_genes_Summary.csv\n")

# ============================================================
# 9) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE — T2D SUSCEPTIBILITY GENES\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat("  Combined Figure:\n")
cat("    - Fig_2g_T2D_susceptibility_genes_Combined.png (600 dpi, 7.2×5.7 cm)\n")
cat("    - Fig_2g_T2D_susceptibility_genes_Combined.svg (vector format)\n")
cat("  Data tables:\n")
cat("    - Fig_2g_T2D_susceptibility_genes_Posthoc.csv (Wilcoxon, BH-corrected within gene)\n")
cat("    - Fig_2g_T2D_susceptibility_genes_Summary.csv\n")
cat("\nKey features:\n")
cat("  ✓ Statistics run BEFORE plot (plot labels match CSV exactly)\n")
cat("  ✓ droplevels() guards against empty factor level errors\n")
cat("  ✓ Bracket order matches my_comparisons via z_order factor-level sort\n")
cat("  ✓ stat_pvalue_manual() with p.adj.signif (BH-corrected)\n")
cat("  ✓ Shared y-axis scale across both panels\n")
cat("  ✓ Bracket positions anchored to global y_max (identical across panels)\n")
cat("  ✓ Y-axis title on left panel (Zfand6); suppressed on right (Srr)\n")
cat("  ✓ Combined 2-column layout: 7.2×5.7 cm\n")
cat("\nStatistical approach:\n")
cat("  • Test: Wilcoxon rank-sum (independent samples, between-group)\n")
cat("  • Correction: Benjamini-Hochberg within each gene\n")
cat("  • Significance: ***p<0.001, **p<0.01, *p<0.05, ns\n")
cat("  • Bracket order: CTRL/IL10-KO → CTRL/DSS-IL10-KO → IL10-KO/DSS-IL10-KO\n")
cat(rep("=", 60), "\n", sep = "")
