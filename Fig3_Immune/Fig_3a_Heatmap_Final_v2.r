# ============================================================
# Cytokine Expression Heatmap - Annotated with Clusters
# Diabetologia Journal Submission - FINAL VERSION
# Features: Descriptive cluster labels, Diabetologia colour palette,
#           Top legend, cell value annotations, SVG export
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/05-Lab members/Ruby/IL10/Paper submission/Figures/Fig3")

# 1) Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(svglite)  # For SVG export

cat("\n=== CYTOKINE HEATMAP - DIABETOLOGIA STYLE ===\n")
cat("Packages loaded successfully\n")

# 2) Read cytokine data
cat("\n=== LOADING DATA ===\n")
dat <- read_csv("Fig_3a_Data_Cytokine.csv")

cat("\nData loaded:\n")
cat("  Cytokines:", length(unique(dat$Cytokines)), "\n")
cat("  Groups:", paste(unique(dat$Group), collapse=", "), "\n")
cat("  Total observations:", nrow(dat), "\n")

# Set group order
dat$Group <- factor(dat$Group, 
                    levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))

# ============================================================
# 3) Calculate Log2(Relative intensity)
# ============================================================
cat("\n=== CALCULATING LOG2(RELATIVE INTENSITY) ===\n")

dat_log2 <- dat %>%
  mutate(Log2_RelativeIntensity = log2(Relative_intensity))

cat("Log2 range: [", round(min(dat_log2$Log2_RelativeIntensity), 2), 
    ",", round(max(dat_log2$Log2_RelativeIntensity), 2), "]\n")

# ============================================================
# 4) Prepare matrix for pheatmap
# ============================================================
cat("\n=== PREPARING MATRIX ===\n")

# Wide format: rows = groups, columns = cytokines
mat_log2 <- dat_log2 %>%
  select(Cytokines, Group, Log2_RelativeIntensity) %>%
  pivot_wider(names_from = Cytokines,
              values_from = Log2_RelativeIntensity) %>%
  arrange(Group)

# Convert to matrix
mat <- as.matrix(mat_log2[, -1])
rownames(mat) <- mat_log2$Group

cat("Matrix: ", nrow(mat), "groups ×", ncol(mat), "cytokines\n")

# ============================================================
# 5) Perform hierarchical clustering
# ============================================================
cat("\n=== PERFORMING CLUSTERING ===\n")

hc_cols <- hclust(dist(t(mat), method = "euclidean"), method = "complete")
cytokine_order <- colnames(mat)[hc_cols$order]

# Cut tree to get clusters
n_clusters <- 4
clusters <- cutree(hc_cols, k = n_clusters)

cat("Number of clusters:", n_clusters, "\n")
cat("Cytokines per cluster:\n")
print(table(clusters))

# ============================================================
# 6) Define descriptive cluster labels
# ============================================================
cat("\n=== ASSIGNING DESCRIPTIVE CLUSTER LABELS ===\n")

# Descriptive labels for biological interpretation
cluster_labels <- c(
  "DSS-normalised",        # Cluster 1
  "Persistently down",     # Cluster 2
  "Strongly upregulated",  # Cluster 3
  "Mild/unchanged"         # Cluster 4
)

cat("Cluster labels:\n")
for (i in 1:n_clusters) {
  cat("  Cluster", i, ":", cluster_labels[i], "\n")
}

# Create annotation with descriptive labels
annotation_col <- data.frame(
  Cluster = factor(
    cluster_labels[clusters[colnames(mat)]],
    levels = cluster_labels
  )
)
rownames(annotation_col) <- colnames(mat)

# ============================================================
# 7) Diabetologia colour palette
# ============================================================
cat("\n=== USING DIABETOLOGIA COLOUR PALETTE ===\n")

# Diabetologia colours (from official colour palette PDF)
# Blues: R:138 G:169 B:214, R:96 G:108 B:158
# Reds: R:221 G:83 B:82, R:193 G:88 B:88
# Greens: R:109 G:165 B:103, R:152 G:201 B:139
# Purples: R:117 G:90 B:145, R:146 G:112 B:165
# Yellows: R:244 G:231 B:100, R:251 G:241 B:173

# Heatmap diverging colour scale (Blue-White-Red)
# Using Diabetologia blues and reds
breaks <- seq(-2, 2, length.out = 101)
colors <- colorRampPalette(c(
  rgb(96, 108, 158, maxColorValue=255),    # Dark Diabetologia blue
  rgb(138, 169, 214, maxColorValue=255),   # Medium Diabetologia blue
  rgb(202, 218, 232, maxColorValue=255),   # Light Diabetologia blue
  "#FFFFFF",                                # White (no change)
  rgb(244, 184, 176, maxColorValue=255),   # Light Diabetologia red
  rgb(221, 83, 82, maxColorValue=255),     # Medium Diabetologia red
  rgb(193, 88, 88, maxColorValue=255)      # Dark Diabetologia red
))(100)

# Cluster annotation colours (Diabetologia palette)
cluster_colors <- list(
  Cluster = c(
    "DSS-normalised" = rgb(244, 231, 100, maxColorValue=255),      # Diabetologia yellow
    "Persistently down" = rgb(138, 169, 214, maxColorValue=255),   # Diabetologia blue
    "Strongly upregulated" = rgb(193, 88, 88, maxColorValue=255),  # Diabetologia red
    "Mild/unchanged" = rgb(152, 201, 139, maxColorValue=255)       # Diabetologia green
  )
)

cat("Heatmap: Blue-White-Red diverging scale (Diabetologia colours)\n")
cat("Clusters: Yellow, Blue, Red, Green (Diabetologia colours)\n")

# Create cluster assignment table with descriptive labels
cluster_df <- data.frame(
  Cytokines = names(clusters),
  Cluster_Number = clusters,
  Cluster_Label = cluster_labels[clusters],
  Order = match(names(clusters), cytokine_order)
) %>%
  arrange(Order)

# ============================================================
# 8) Create annotated heatmap - PNG VERSION
# ============================================================
cat("\n=== GENERATING HEATMAP (PNG) ===\n")

png("Fig_3a_Pheatmap_With_Clusters.png", 
    width = 12, height = 5, units = "in", res = 600)

pheatmap(
  mat,
  color = colors,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  clustering_distance_cols = "euclidean",
  scale = "none",
  fontsize_row = 9,
  fontsize_col = 7,
  display_numbers = round(mat, 2),
  number_color = "black",
  fontsize_number = 4.5,
  legend = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  main = "",
  border_color = "grey90",
  annotation_col = annotation_col,
  annotation_colors = cluster_colors,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  cellwidth = NA,
  cellheight = 30
)

# Add custom legend title
grid.text("Log2(Relative intensity)", 
          x = 0.92, y = 0.85, 
          gp = gpar(fontsize = 10, fontface = "bold"))

dev.off()
cat("✓ PNG saved: Fig_3a_Pheatmap_With_Clusters.png (600 dpi)\n")

# ============================================================
# 9) Create annotated heatmap - PDF VERSION
# ============================================================
cat("\n=== GENERATING HEATMAP (PDF) ===\n")

pdf("Fig_3a_Pheatmap_With_Clusters.pdf", 
    width = 12, height = 5)

pheatmap(
  mat,
  color = colors,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  clustering_distance_cols = "euclidean",
  scale = "none",
  fontsize_row = 9,
  fontsize_col = 7,
  display_numbers = round(mat, 2),
  number_color = "black",
  fontsize_number = 4.5,
  legend = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  main = "",
  border_color = "grey90",
  annotation_col = annotation_col,
  annotation_colors = cluster_colors,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  cellwidth = NA,
  cellheight = 30
)

# Add custom legend title
grid.text("Log2(Relative intensity)", 
          x = 0.92, y = 0.85, 
          gp = gpar(fontsize = 10, fontface = "bold"))

dev.off()
cat("✓ PDF saved: Fig_3a_Pheatmap_With_Clusters.pdf (vector)\n")

# ============================================================
# 10) Create annotated heatmap - SVG VERSION
# ============================================================
cat("\n=== GENERATING HEATMAP (SVG) ===\n")

svglite("Fig_3a_Pheatmap_With_Clusters.svg", 
        width = 12, height = 5)

pheatmap(
  mat,
  color = colors,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  clustering_distance_cols = "euclidean",
  scale = "none",
  fontsize_row = 9,
  fontsize_col = 7,
  display_numbers = round(mat, 2),
  number_color = "black",
  fontsize_number = 4.5,
  legend = TRUE,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  main = "",
  border_color = "grey90",
  annotation_col = annotation_col,
  annotation_colors = cluster_colors,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  cellwidth = NA,
  cellheight = 30
)

# Add custom legend title
grid.text("Log2(Relative intensity)", 
          x = 0.92, y = 0.85, 
          gp = gpar(fontsize = 10, fontface = "bold"))

dev.off()
cat("✓ SVG saved: Fig_3a_Pheatmap_With_Clusters.svg (vector, editable)\n")

# ============================================================
# 11) Export essential data tables
# ============================================================
cat("\n=== EXPORTING DATA TABLES ===\n")

# 1. Cluster assignments with descriptive labels
write.csv(cluster_df, "Fig_3a_Cluster_Assignments.csv", row.names = FALSE)
cat("✓ Cluster assignments exported\n")

# 2. Full log2 data (long format)
log2_data <- dat_log2 %>%
  select(Cytokines, Group, Relative_intensity, Log2_RelativeIntensity) %>%
  mutate(Log2_RelativeIntensity = round(Log2_RelativeIntensity, 3),
         Relative_intensity = round(Relative_intensity, 3)) %>%
  arrange(Cytokines, Group)

write.csv(log2_data, "Fig_3a_Log2_Data.csv", row.names = FALSE)
cat("✓ Log2 data exported\n")

# 3. Comparison table (wide format, with cluster labels, in clustered order)
comparison_table <- dat_log2 %>%
  select(Cytokines, Group, Log2_RelativeIntensity) %>%
  pivot_wider(names_from = Group, values_from = Log2_RelativeIntensity) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  left_join(cluster_df %>% select(Cytokines, Cluster_Number, Cluster_Label), 
            by = "Cytokines") %>%
  arrange(match(Cytokines, cytokine_order)) %>%
  select(Cytokines, Cluster_Label, CTRL, `IL10-KO`, `DSS-IL10-KO`)

write.csv(comparison_table, "Fig_3a_Comparison_Table.csv", row.names = FALSE)
cat("✓ Comparison table exported (clustered order)\n")

# 4. Cluster characterization with descriptive labels
cluster_characterization <- dat_log2 %>%
  left_join(cluster_df %>% select(Cytokines, Cluster_Label), by = "Cytokines") %>%
  group_by(Cluster_Label, Group) %>%
  summarise(
    n_cytokines = n(),
    mean_log2 = round(mean(Log2_RelativeIntensity), 3),
    median_log2 = round(median(Log2_RelativeIntensity), 3),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = Group, 
              values_from = c(n_cytokines, mean_log2, median_log2)) %>%
  select(Cluster_Label, starts_with("n_cytokines_CTRL"), 
         starts_with("mean_log2"), starts_with("median_log2"))

write.csv(cluster_characterization, 
          "Fig_3a_Cluster_Characterization.csv", 
          row.names = FALSE)
cat("✓ Cluster characterization exported\n")

# 5. Top changed cytokines in DSS-IL10-KO
top_changes <- dat_log2 %>%
  filter(Group == "DSS-IL10-KO") %>%
  left_join(cluster_df %>% select(Cytokines, Cluster_Label), by = "Cytokines") %>%
  arrange(desc(abs(Log2_RelativeIntensity))) %>%
  head(20) %>%
  select(Cytokines, Cluster_Label, Relative_intensity, Log2_RelativeIntensity) %>%
  mutate(
    Relative_intensity = round(Relative_intensity, 3),
    Log2_RelativeIntensity = round(Log2_RelativeIntensity, 3),
    Fold_vs_CTRL = round(2^Log2_RelativeIntensity, 2),
    Direction = ifelse(Log2_RelativeIntensity > 0, "Higher", "Lower")
  )

write.csv(top_changes, "Fig_3a_Top_Changed_DSS.csv", row.names = FALSE)
cat("✓ Top changed cytokines exported\n")

# ============================================================
# 12) Summary statistics
# ============================================================
cat("\n=== SUMMARY STATISTICS ===\n")

# Group summary
group_summary <- dat_log2 %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean_log2 = round(mean(Log2_RelativeIntensity), 3),
    median_log2 = round(median(Log2_RelativeIntensity), 3),
    sd_log2 = round(sd(Log2_RelativeIntensity), 3),
    min_log2 = round(min(Log2_RelativeIntensity), 2),
    max_log2 = round(max(Log2_RelativeIntensity), 2),
    n_higher_1.5fold = sum(Log2_RelativeIntensity > log2(1.5)),
    n_lower_1.5fold = sum(Log2_RelativeIntensity < log2(1/1.5)),
    .groups = 'drop'
  )

cat("\nGroup summary:\n")
print(group_summary)

write.csv(group_summary, "Fig_3a_Group_Summary.csv", row.names = FALSE)
cat("\n✓ Group summary exported\n")

# ============================================================
# 13) Final summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("HEATMAP GENERATION COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")

cat("\nGenerated files:\n")
cat("\nHeatmaps (3 formats):\n")
cat("  1. Fig_3a_Pheatmap_With_Clusters.png (600 dpi)\n")
cat("  2. Fig_3a_Pheatmap_With_Clusters.pdf (vector)\n")
cat("  3. Fig_3a_Pheatmap_With_Clusters.svg (vector, editable)\n")

cat("\nData tables (6 files):\n")
cat("  4. Fig_3a_Cluster_Assignments.csv\n")
cat("  5. Fig_3a_Log2_Data.csv\n")
cat("  6. Fig_3a_Comparison_Table.csv\n")
cat("  7. Fig_3a_Cluster_Characterization.csv\n")
cat("  8. Fig_3a_Top_Changed_DSS.csv\n")
cat("  9. Fig_3a_Group_Summary.csv\n")

cat("\nTotal: 9 files (3 images + 6 data tables)\n")

cat("\nHeatmap features:\n")
cat("  • Layout: Cytokines (columns) × Groups (rows)\n")
cat("  • Clustering: Complete linkage, Euclidean distance\n")
cat("  • Number of clusters:", n_clusters, "\n")
cat("  • Cluster labels: Descriptive (biological interpretation)\n")
cat("  • Colour scale: Diabetologia blue-white-red (-2 to +2)\n")
cat("  • Cluster colours: Diabetologia palette\n")
cat("  • Cell annotations: Log2 values displayed (2 decimals)\n")
cat("  • Legend: Top position with custom title\n")

cat("\nCluster summary (descriptive labels):\n")
for (i in 1:n_clusters) {
  cytokines_in_cluster <- cluster_df %>% 
    filter(Cluster_Number == i) %>% 
    pull(Cytokines)
  cat("  ", cluster_labels[i], ": n =", length(cytokines_in_cluster), "cytokines\n")
}

cat("\nTop 10 most changed in DSS-IL10-KO:\n")
print(head(top_changes, 10))

cat("\n", rep("=", 60), "\n", sep = "")
cat("READY FOR DIABETOLOGIA SUBMISSION\n")
cat(rep("=", 60), "\n\n", sep = "")
