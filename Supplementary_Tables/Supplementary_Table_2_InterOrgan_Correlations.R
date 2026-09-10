# ============================================================
# Supplementary Table 2 - Benjamini-Hochberg-corrected and
# group-stratified reanalysis of inter-organ gene expression
# correlations
# Scientific Reports Rev1 rebuttal (Reviewer 2 point 4 /
# Reviewer 3's principal statistical objection)
#
# Related to Fig. 5 and the Results section "Cross-organ
# transcriptional associations reveal organ-specific
# co-expression patterns."
#
# Method:
#   Pooled analysis: for each pair of organs, Pearson
#   correlations are computed gene-by-gene across all three
#   groups combined (N <= 12 per pair), matching the approach
#   used for the original Fig. 5 panels.
#   Stratified analysis: the same gene-pair correlations are
#   recalculated separately within each of the three groups
#   (n <= 4 per pair).
#   Benjamini-Hochberg (BH) correction is applied across all
#   287 pooled inter-organ gene-pair tests (the stratified
#   tests are reported at nominal/raw p only, as in the
#   manuscript).
#
# IMPORTANT - animal pairing across organs:
#   The per-organ CSVs used here do not carry a shared animal
#   ID; within each Group x Gene stratum, replicate rows are
#   paired across organs by ROW ORDER (i.e. the 1st CTRL
#   replicate measured for one organ's gene is paired with the
#   1st CTRL replicate measured for another organ's gene, and
#   so on), matching the convention already used to build the
#   pooled Fig. 5 dataset. This pairing is NOT verified against
#   a shared per-animal identifier and should be confirmed
#   against original laboratory records (necropsy/replicate
#   logs) before the "within-animal" framing in the manuscript
#   text is relied upon; the correlations themselves are
#   correct as group-level associations regardless.
#
# Output:
#   Rebuttal_Pooled_InterOrgan_Correlations_BHcorrected.csv
#   Rebuttal_Stratified_InterOrgan_Correlations_PerGroup.csv
# ============================================================

# 0) Clear workspace and set working directory
rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# 1) Load required libraries
required_packages <- c("dplyr", "tidyr", "readr", "purrr", "tibble")
missing_packages  <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)

# ============================================================
# 2) Import per-organ gene expression data
#    Only the five main-text genes per organ (matching the
#    manuscript figure legends) are retained.
# ============================================================

files <- list(
  Colon    = "../Fig1_Intestinal/Fig_1f_Data_LargeIntestine_GeneExpr.csv",
  Pancreas = "../Fig2_Pancreas/Fig_2ef_Pancreas_FoldChange_Data.csv",
  Spleen   = "../Fig3_Immune/Fig_3c_Data_Spleen_GeneExpr.csv",
  Liver    = "../Fig4_Organs/Fig_4b_Data_Liver_GeneExpr.csv",
  Kidney   = "../Fig4_Organs/Fig_4e_Data_Kidney_GeneExpr.csv"
)

# genes actually plotted in the main-text figure for each organ
figure_genes <- list(
  Colon    = c("Mafa", "Glp1r", "Insr", "Zfand6", "Tcf19"),
  Pancreas = c("Mafa", "Ins1", "Ins2", "Insr", "Tbc1d4", "Glut1", "Zfand6", "Srr"),
  Spleen   = c("Mafa", "Insr", "Glut1", "Zfand6", "Tcf19"),
  Liver    = c("Tbc1d4", "Glut1", "Zfand6", "Tcf19", "Ghrhr"),
  Kidney   = c("Mafa", "Glut1", "Zfand6", "Insr")
)

load_organ <- function(path, organ) {
  # imap_dfr() calls .f(value, name); for a named list `value` is the
  # list element (path) and `name` is the list name (organ) - keep this
  # argument order in sync with the imap_dfr() call below.
  df <- read_csv(path, show_col_types = FALSE)
  names(df) <- trimws(names(df))
  df <- df %>%
    mutate(Gene = recode(Gene, "SRR" = "Srr")) %>%
    select(Group, Gene, Relative_expression) %>%
    filter(Gene %in% figure_genes[[organ]]) %>%
    mutate(Organ = organ) %>%
    group_by(Group, Gene) %>%
    mutate(Sample_ID = row_number() - 1) %>%
    ungroup()
  df
}

all_long <- imap_dfr(files, load_organ) %>%
  mutate(Organ_Gene = paste(Organ, Gene, sep = "_"))

cat("Total observations loaded:", nrow(all_long), "\n")

# ============================================================
# 3) Reshape to wide (one row per Group x Sample_ID) and
#    compute every inter-organ gene-pair Pearson correlation
# ============================================================

wide_for <- function(df_subset) {
  df_subset %>%
    group_by(Group, Sample_ID, Organ_Gene) %>%
    summarise(Relative_expression = mean(Relative_expression), .groups = "drop") %>%
    pivot_wider(names_from = Organ_Gene, values_from = Relative_expression)
}

pairwise_corr <- function(wide) {
  cols <- setdiff(names(wide), c("Group", "Sample_ID"))
  combos <- combn(cols, 2, simplify = FALSE)
  map_dfr(combos, function(cc) {
    c1 <- cc[1]; c2 <- cc[2]
    o1 <- strsplit(c1, "_")[[1]][1]
    o2 <- strsplit(c2, "_")[[1]][1]
    if (o1 == o2) return(NULL)  # inter-organ pairs only
    pair_df <- wide[, c(c1, c2)]
    pair_df <- pair_df[complete.cases(pair_df), ]
    n <- nrow(pair_df)
    if (n < 3) {
      r <- NA_real_; p <- NA_real_
    } else {
      ct <- suppressWarnings(cor.test(pair_df[[c1]], pair_df[[c2]], method = "pearson"))
      r <- unname(ct$estimate); p <- ct$p.value
    }
    tibble(
      Organ1 = o1, Organ2 = o2, Col1 = c1, Col2 = c2,
      Gene1 = sub("^[^_]+_", "", c1), Gene2 = sub("^[^_]+_", "", c2),
      N = n, Correlation = r, P_value = p
    )
  })
}

# --- POOLED (all groups combined, matches current manuscript Fig. 5 approach) ---
wide_pooled <- wide_for(all_long)
pooled <- pairwise_corr(wide_pooled) %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value) %>%
  mutate(
    rank        = row_number(),
    BH_threshold = rank / n() * 0.05,
    p_adj_BH    = p.adjust(P_value, method = "BH"),
    Sig_raw     = P_value < 0.05,
    Sig_BH      = p_adj_BH < 0.05
  )

cat("\n=== POOLED (N<=12 across CTRL/IL10-KO/DSS-IL10-KO combined) ===\n")
cat("Total inter-organ gene-pair tests:", nrow(pooled), "\n")
cat("Significant at raw p<0.05:        ", sum(pooled$Sig_raw), "\n")
cat("Significant after BH correction:  ", sum(pooled$Sig_BH), "\n\n")

cat("By organ pair (raw-significant / BH-significant / total tested):\n")
pooled %>%
  group_by(Organ1, Organ2) %>%
  summarise(total = n(), sig_raw = sum(Sig_raw), sig_bh = sum(Sig_BH), .groups = "drop") %>%
  print(n = Inf)

# --- STRATIFIED (within each group separately, n<=4) ---
strat <- all_long %>%
  group_split(Group) %>%
  map_dfr(function(sub) {
    grp <- unique(sub$Group)
    wide_g <- wide_for(sub)
    pairwise_corr(wide_g) %>%
      filter(!is.na(P_value)) %>%
      mutate(Group = grp, Sig_raw = P_value < 0.05)
  })

cat("\n=== STRATIFIED BY GROUP (within-group, n<=4 per group) ===\n")
strat %>%
  group_by(Group) %>%
  summarise(total = n(), sig_raw = sum(Sig_raw),
            pct = round(100 * sig_raw / total, 1), .groups = "drop") %>%
  print(n = Inf)

# Cross-check: of the pooled raw-significant pairs, how many are ALSO
# raw-significant in at least one individual group stratum?
pooled_sig <- pooled %>% filter(Sig_raw)
strat_sig  <- strat  %>% filter(Sig_raw)
pooled_sig_keys <- paste(pooled_sig$Col1, pooled_sig$Col2)
strat_sig_keys  <- unique(paste(strat_sig$Col1, strat_sig$Col2))
overlap <- intersect(pooled_sig_keys, strat_sig_keys)
cat(sprintf(
  "\nOf %d pooled raw-significant pairs, %d are also nominally significant (p<0.05) within at least one individual group.\n",
  length(pooled_sig_keys), length(overlap)
))

# ============================================================
# 4) Export
#    Written back into Fig5_InterOrgan/, alongside the rest of
#    that figure's source data, not into this script's own
#    folder.
# ============================================================
out_pooled <- "../Fig5_InterOrgan/Rebuttal_Pooled_InterOrgan_Correlations_BHcorrected.csv"
out_strat  <- "../Fig5_InterOrgan/Rebuttal_Stratified_InterOrgan_Correlations_PerGroup.csv"
write_csv(pooled, out_pooled)
write_csv(strat,  out_strat)

cat("\nSaved:\n")
cat(" -", out_pooled, "\n")
cat(" -", out_strat, "\n")

# ============================================================
# END OF SCRIPT
# ============================================================
