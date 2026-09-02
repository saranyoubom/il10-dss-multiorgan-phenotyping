# ============================================================
# Spleen Gene Correlation Analysis - Fig_3d
# Pearson Correlation Matrix between Genes
# Diabetologia Journal Submission
# Five genes: Mafa, Insr, Glut1, Zfand6, Tcf19
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
library(corrplot)
library(Hmisc)
library(scales)

# 2) Define genes for analysis
genes_fig3d <- c("Mafa", "Insr", "Glut1", "Zfand6", "Tcf19")

# ============================================================
# 3) Import and prepare data
# ============================================================
cat("\n=== LOADING SPLEEN GENE EXPRESSION DATA ===\n")
dat_raw <- read.csv("Fig_3c_Data_Spleen_GeneExpr.csv", stringsAsFactors = FALSE)

dat <- dat_raw %>%
  filter(Gene %in% genes_fig3d) %>%
  select(Group, Gene, Relative_expression) %>%
  filter(!is.na(Relative_expression))

dat$Relative_expression <- as.numeric(dat$Relative_expression)
dat <- dat %>% filter(!is.na(Relative_expression))

cat("Total observations (clean data):", nrow(dat), "\n")
cat("Genes:", paste(genes_fig3d, collapse = ", "), "\n")
cat("\nSample sizes per gene:\n")
print(table(dat$Gene))

# ============================================================
# 4) Reshape data to wide format (samples × genes)
# ============================================================
cat("\n=== RESHAPING DATA FOR CORRELATION ANALYSIS ===\n")

dat <- dat %>%
  group_by(Gene, Group) %>%
  mutate(Sample_ID = paste(Group, row_number(), sep = "_")) %>%
  ungroup()

dat_wide <- dat %>%
  select(Sample_ID, Gene, Relative_expression) %>%
  pivot_wider(names_from = Gene, values_from = Relative_expression)

dat_matrix <- dat_wide %>%
  select(-Sample_ID) %>%
  select(all_of(genes_fig3d))

complete_samples <- sum(complete.cases(dat_matrix))
total_samples    <- nrow(dat_matrix)

cat("\nTotal samples:", total_samples, "\n")
cat("Complete cases (no NA):", complete_samples, "\n")

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
cor_n <- cor_test$n   # n per pair

cat("\nP-value matrix:\n")
print(round(cor_p, 4))

cat("\nSample size matrix (n per pair):\n")
print(cor_n)

# ============================================================
# 6) Prepare data for ggplot heatmap
# ============================================================
cat("\n=== PREPARING DATA FOR HEATMAP ===\n")

cor_melted <- melt(cor_r)
colnames(cor_melted) <- c("Gene1", "Gene2", "Correlation")

p_melted <- melt(cor_p)
colnames(p_melted) <- c("Gene1", "Gene2", "P_value")

n_melted <- melt(cor_n)
colnames(n_melted) <- c("Gene1", "Gene2", "N_pair")

plot_data <- cor_melted %>%
  left_join(p_melted, by = c("Gene1", "Gene2")) %>%
  left_join(n_melted, by = c("Gene1", "Gene2"))

plot_data <- plot_data %>%
  mutate(
    Signif = case_when(
      is.na(P_value)  ~ "",
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    Label = case_when(
      is.na(P_value) ~ sprintf("%.2f", Correlation),            # diagonal
      Signif == "ns" ~ sprintf("%.2f", Correlation),            # ns: r only
      TRUE           ~ sprintf("%s\n%.2f", Signif, Correlation) # sig: symbol on top
    )
  )

plot_data$Gene1 <- factor(plot_data$Gene1, levels = genes_fig3d)
plot_data$Gene2 <- factor(plot_data$Gene2, levels = rev(genes_fig3d))

cat("Heatmap data prepared\n")
cat("Range of correlations:", round(min(plot_data$Correlation), 3),
    "to", round(max(plot_data$Correlation), 3), "\n")

# ============================================================
# 7) Create correlation heatmap
# ============================================================
cat("\n=== GENERATING CORRELATION HEATMAP ===\n")

col_palette <- scale_fill_gradient2(
  low      = "#7883ba",
  mid      = "white",
  high     = "#ed7374",
  midpoint = 0,
  limits   = c(-1, 1),
  breaks   = seq(-1, 1, 0.5),
  name     = "Pearson\ncorrelation"
)

p_heatmap <- ggplot(plot_data, aes(x = Gene1, y = Gene2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 2.5, color = "black", lineheight = 0.75) +
  col_palette +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(size = 8, colour = "black",
                                   angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
    axis.text.y     = element_text(size = 8, colour = "black",
                                   angle = 0, hjust = 1, face = "italic"),
    axis.ticks      = element_blank(),
    axis.line       = element_blank(),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    plot.margin     = margin(0, 0, 0, 0, "mm")
  ) +
  coord_fixed()

print(p_heatmap)

# ============================================================
# 8) Save plot
# ============================================================
cat("\n=== SAVING PLOT FILES ===\n")

ggsave(
  filename = "Fig_3d_Correlation_Matrix_Spleen.png",
  plot     = p_heatmap,
  scale    = 1, width = 7, height = 6,
  units    = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename = "Fig_3d_Correlation_Matrix_Spleen.svg",
  plot     = p_heatmap,
  scale    = 1, width = 7, height = 6,
  units    = "cm", limitsize = FALSE
)

cat("✓ Plots saved\n")

# ============================================================
# 9) Export correlation tables
# ============================================================
cat("\n=== EXPORTING CORRELATION TABLES ===\n")

cor_export <- data.frame(Gene = rownames(cor_r), round(cor_r, 3))
p_export   <- data.frame(Gene = rownames(cor_p), round(cor_p, 4))
n_export   <- data.frame(Gene = rownames(cor_n), cor_n)

write.csv(cor_export, "Fig_3d_Correlation_Matrix_Spleen.csv",    row.names = FALSE)
write.csv(p_export,   "Fig_3d_PValue_Matrix_Spleen.csv",         row.names = FALSE)
write.csv(n_export,   "Fig_3d_N_Matrix_Spleen.csv",              row.names = FALSE)

pairwise_cor <- plot_data %>%
  filter(as.character(Gene1) < as.character(Gene2)) %>%
  select(Gene1, Gene2, N_pair, Correlation, P_value, Signif) %>%
  arrange(desc(abs(Correlation)))

write.csv(pairwise_cor, "Fig_3d_Pairwise_Correlations_Spleen.csv", row.names = FALSE)

cat("\nPairwise correlations (sorted by absolute value):\n")
print(as.data.frame(pairwise_cor))

# ============================================================
# 10) Summary statistics
# ============================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("CORRELATION ANALYSIS SUMMARY - SPLEEN GENES\n")
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
print(as.data.frame(pairwise_cor %>% filter(Correlation == max(Correlation))))

cat("\nStrongest negative correlation:\n")
print(as.data.frame(pairwise_cor %>% filter(Correlation == min(Correlation))))

cat("\n", rep("=", 70), "\n", sep = "")
cat("FILES GENERATED\n")
cat(rep("=", 70), "\n", sep = "")
cat("Figures:\n")
cat("  - Fig_3d_Correlation_Matrix_Spleen.png (600 dpi, 7×6 cm)\n")
cat("  - Fig_3d_Correlation_Matrix_Spleen.svg (vector format)\n")
cat("\nData tables:\n")
cat("  - Fig_3d_Correlation_Matrix_Spleen.csv    (correlation coefficients)\n")
cat("  - Fig_3d_PValue_Matrix_Spleen.csv         (p-values)\n")
cat("  - Fig_3d_N_Matrix_Spleen.csv              (n per pair)\n")
cat("  - Fig_3d_Pairwise_Correlations_Spleen.csv (pairwise results with n)\n")
cat("\nInterpretation:\n")
cat("  • Blue  = negative correlation\n")
cat("  • White = no correlation (r = 0)\n")
cat("  • Red   = positive correlation\n")
cat("  • Significance: ***p<0.001, **p<0.01, *p<0.05, ns = not significant\n")
cat("\nNote:\n")
cat("  • All 5 genes had complete data (n=12 per gene, 3 groups × 4 replicates)\n")
cat("  • All samples included in correlation analysis\n")
cat(rep("=", 70), "\n", sep = "")
