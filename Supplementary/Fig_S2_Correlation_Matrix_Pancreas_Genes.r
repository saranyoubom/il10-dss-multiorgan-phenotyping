# ============================================================
# Pancreas Gene Correlation Analysis - Fig_S2
# Pearson Correlation Matrix between Genes
# Diabetologia Journal Submission (Supplementary)
# Eight genes: Mafa, Ins1, Ins2, Insr, Tbc1d4, Glut1, Zfand6, Srr
# Data: Relative gene expression across all groups
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# 1) Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(corrplot)  # For additional correlation utilities
library(Hmisc)     # For correlation p-values
library(scales)    # For colour scales

# 2) Define genes for analysis
genes_figS2 <- c("Mafa", "Ins1", "Ins2", "Insr", "Tbc1d4", "Glut1", "Zfand6", "Srr")

# ============================================================
# 3) Import and prepare data
# ============================================================
cat("\n=== LOADING PANCREAS GENE EXPRESSION DATA ===\n")
dat_raw <- read.csv("Fig_S2_Pancreas_mRNA_Expression.csv", stringsAsFactors = FALSE)

# Filter for the 8 selected genes
dat <- dat_raw %>%
  filter(Gene %in% genes_figS2) %>%
  select(Group, Gene, Relative_expression) %>%
  filter(!is.na(Relative_expression))

# Convert ND to NA and then to numeric
dat$Relative_expression <- ifelse(dat$Relative_expression == "ND", NA, dat$Relative_expression)
dat$Relative_expression <- as.numeric(dat$Relative_expression)

# Remove NA values
dat <- dat %>% filter(!is.na(Relative_expression))

cat("Total observations (after removing ND):", nrow(dat), "\n")
cat("Genes:", paste(genes_figS2, collapse = ", "), "\n")

cat("\nSample sizes per gene:\n")
print(table(dat$Gene))

# ============================================================
# 4) Reshape data to wide format (samples × genes)
# ============================================================
cat("\n=== RESHAPING DATA FOR CORRELATION ANALYSIS ===\n")

# Create sample IDs based on group and replicate number
dat <- dat %>%
  group_by(Gene, Group) %>%
  mutate(Sample_ID = paste(Group, row_number(), sep = "_")) %>%
  ungroup()

# Pivot to wide format
dat_wide <- dat %>%
  select(Sample_ID, Gene, Relative_expression) %>%
  pivot_wider(names_from = Gene, values_from = Relative_expression)

# Remove Sample_ID column for correlation
dat_matrix <- dat_wide %>%
  select(-Sample_ID) %>%
  select(all_of(genes_figS2))  # Ensure consistent gene order

# Check for complete cases
complete_samples <- sum(complete.cases(dat_matrix))
total_samples    <- nrow(dat_matrix)

cat("\nTotal samples:", total_samples, "\n")
cat("Complete cases (no NA):", complete_samples, "\n")

# Use only complete cases for correlation
dat_matrix_complete <- dat_matrix[complete.cases(dat_matrix), ]

cat("Using", nrow(dat_matrix_complete), "complete samples for correlation analysis\n")

if (nrow(dat_matrix_complete) < 3) {
  stop("ERROR: Insufficient complete cases for correlation analysis (need at least 3)")
}

# ============================================================
# 5) Calculate Pearson correlation matrix with p-values
# ============================================================
cat("\n=== CALCULATING PEARSON CORRELATIONS ===\n")

cor_matrix <- cor(dat_matrix_complete, method = "pearson", use = "complete.obs")

cat("\nPearson correlation matrix:\n")
print(round(cor_matrix, 3))

cor_test <- rcorr(as.matrix(dat_matrix_complete), type = "pearson")

cor_r <- cor_test$r
cor_p <- cor_test$P
cor_n <- cor_test$n   # <-- ADD: n per pair (matrix)

cat("\nP-value matrix:\n")
print(round(cor_p, 4))

cat("\nSample size matrix (n per pair):\n")
print(cor_n)          # <-- ADD: print n matrix for verification

# ============================================================
# 6) Prepare data for ggplot heatmap
# ============================================================
cat("\n=== PREPARING DATA FOR HEATMAP ===\n")

cor_melted <- melt(cor_r)
colnames(cor_melted) <- c("Gene1", "Gene2", "Correlation")

p_melted <- melt(cor_p)
colnames(p_melted) <- c("Gene1", "Gene2", "P_value")

n_melted <- melt(cor_n)                               # <-- ADD
colnames(n_melted) <- c("Gene1", "Gene2", "N_pair")   # <-- ADD

plot_data <- cor_melted %>%
  left_join(p_melted, by = c("Gene1", "Gene2")) %>%
  left_join(n_melted, by = c("Gene1", "Gene2"))        # <-- ADD

plot_data <- plot_data %>%
  mutate(
    Signif = case_when(
      is.na(P_value)   ~ "",
      P_value < 0.001  ~ "***",
      P_value < 0.01   ~ "**",
      P_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    ),
    Label = ifelse(is.na(P_value),
                   sprintf("%.2f", Correlation),
                   sprintf("%.2f%s", Correlation,
                           ifelse(Signif == "ns", "", paste0("\n", Signif))))
  )

plot_data$Gene1 <- factor(plot_data$Gene1, levels = genes_figS2)
plot_data$Gene2 <- factor(plot_data$Gene2, levels = rev(genes_figS2))

cat("Heatmap data prepared\n")
cat("Range of correlations:", round(min(plot_data$Correlation), 3),
    "to", round(max(plot_data$Correlation), 3), "\n")

# ============================================================
# 9) Export correlation tables
# ============================================================
cat("\n=== EXPORTING CORRELATION TABLES ===\n")

cor_export <- data.frame(Gene = rownames(cor_r), round(cor_r, 3))
p_export   <- data.frame(Gene = rownames(cor_p), round(cor_p, 4))
n_export   <- data.frame(Gene = rownames(cor_n), cor_n)  # <-- ADD

write.csv(cor_export, "Fig_S2_Correlation_Matrix_Pancreas.csv",    row.names = FALSE)
write.csv(p_export,   "Fig_S2_PValue_Matrix_Pancreas.csv",         row.names = FALSE)
write.csv(n_export,   "Fig_S2_N_Matrix_Pancreas.csv",              row.names = FALSE)  # <-- ADD

# Pairwise correlations — now includes N_pair
pairwise_cor <- plot_data %>%
  filter(as.character(Gene1) < as.character(Gene2)) %>%
  select(Gene1, Gene2, N_pair, Correlation, P_value, Signif) %>%  # <-- N_pair added
  arrange(desc(abs(Correlation)))

write.csv(pairwise_cor, "Fig_S2_Pairwise_Correlations_Pancreas.csv", row.names = FALSE)

cat("\nPairwise correlations (sorted by absolute value):\n")
print(as.data.frame(pairwise_cor))

cat("\nFiles saved:\n")
cat("  - Fig_S2_Correlation_Matrix_Pancreas.csv    (r matrix)\n")
cat("  - Fig_S2_PValue_Matrix_Pancreas.csv         (p-value matrix)\n")
cat("  - Fig_S2_N_Matrix_Pancreas.csv              (n per pair matrix)\n")  # <-- ADD
cat("  - Fig_S2_Pairwise_Correlations_Pancreas.csv (long format with n)\n")


# ============================================================
# 10) Summary statistics
# ============================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("CORRELATION ANALYSIS SUMMARY - PANCREAS GENES\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nNumber of samples used:", nrow(dat_matrix_complete), "\n")
cat("Number of genes:", ncol(dat_matrix_complete), "\n")
cat("Correlation method: Pearson\n")

sig_correlations <- pairwise_cor %>%
  summarise(
    Total    = n(),
    Sig_p001 = sum(P_value < 0.001, na.rm = TRUE),
    Sig_p01  = sum(P_value < 0.01,  na.rm = TRUE),
    Sig_p05  = sum(P_value < 0.05,  na.rm = TRUE),
    NonSig   = sum(P_value >= 0.05, na.rm = TRUE)
  )

cat("\nSignificance summary (", nrow(pairwise_cor), " pairwise comparisons):\n", sep = "")
cat("  p < 0.001 (***): ", sig_correlations$Sig_p001, "\n", sep = "")
cat("  p < 0.01  (**):  ", sig_correlations$Sig_p01,  "\n", sep = "")
cat("  p < 0.05  (*):   ", sig_correlations$Sig_p05,  "\n", sep = "")
cat("  p ≥ 0.05  (ns):  ", sig_correlations$NonSig,   "\n", sep = "")

cat("\nStrongest positive correlation:\n")
print(pairwise_cor %>% filter(Correlation == max(Correlation)))

cat("\nStrongest negative correlation:\n")
print(pairwise_cor %>% filter(Correlation == min(Correlation)))

cat("\n", rep("=", 70), "\n", sep = "")
cat("FILES GENERATED\n")
cat(rep("=", 70), "\n", sep = "")
cat("Figures:\n")
cat("  - Fig_S2_Correlation_Matrix_Pancreas.png (600 dpi, 9×7 cm)\n")
cat("  - Fig_S2_Correlation_Matrix_Pancreas.svg (vector format)\n")
cat("\nData tables:\n")
cat("  - Fig_S2_Correlation_Matrix_Pancreas.csv    (correlation coefficients)\n")
cat("  - Fig_S2_PValue_Matrix_Pancreas.csv         (p-values)\n")
cat("  - Fig_S2_Pairwise_Correlations_Pancreas.csv (pairwise results)\n")
cat("\nInterpretation:\n")
cat("  • Blue  = negative correlation\n")
cat("  • White = no correlation (r = 0)\n")
cat("  • Red   = positive correlation\n")
cat("  • Significance: ***p<0.001, **p<0.01, *p<0.05, ns = not significant\n")
cat("\nNote:\n")
cat("  • ND (not detected) values were excluded from analysis\n")
cat("  • Only complete cases (samples with data for all genes) were used\n")
cat(rep("=", 70), "\n", sep = "")
