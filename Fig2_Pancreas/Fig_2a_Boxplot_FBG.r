# ============================================================
# Fasting Blood Glucose (FBG) Analysis - Enhanced Version
# Paired Within-Subject Comparisons
# Revised for Diabetologia Journal Submission
# Normality Testing for Information Only
# VISUAL APPEARANCE MATCHED TO FIG 1B
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
# Weeks: 0, 6, 12
# ============================================================

# STUDY DESIGN NOTE:
# This is a PAIRED design - same mice measured repeatedly over time
# Has ID column to track individuals
# Uses PAIRED t-test consistently throughout (matching plot)

# 0) Clear workspace and set working directory
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/05-Lab members/Ruby/IL10/Paper submission/Figures/Fig2")

# 1) Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(ggpubr)
library(ggsci)
library(rstatix)
library(tidyr)

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define comparison weeks
my_comparison_weeks <- list(
  c("0", "6"),   # Baseline to Week 6 (first DSS endpoint)
  c("6", "12"),  # Week 6 to 12 (long-term from second DSS)
  c("0", "12")   # Baseline to study end
)

# ============================================================
# 4) Import FBG data
# ============================================================
cat("\n=== LOADING FASTING BLOOD GLUCOSE DATA ===\n")
dat <- read_csv("Fig_2a_Data_FBG.csv")
dat <- dat %>% filter(Replicate %in% c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5"))

# Filter for Weeks 0, 6, 12 only
dat <- dat %>%
  filter(Week_num %in% c(0, 6, 12)) %>%
  filter(!is.na(FBG))

# Convert variables to factors
dat$ID       <- as.factor(dat$ID)
dat$Group    <- as.factor(dat$Group)
dat$Week_num <- as.factor(dat$Week_num)

# Set the order of groups
z_order <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
dat$Group <- factor(dat$Group, levels = z_order)

# Define group labels
Group.labs <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
names(Group.labs) <- c("CTRL", "IL10-KO", "DSS-IL10-KO")

# Check data structure
cat("\n=== DATA STRUCTURE FOR FBG ===\n")
print(head(dat, 20))
cat("\nTotal observations:", nrow(dat), "\n")
cat("\nSample sizes per group per week:\n")
print(table(dat$Group, dat$Week_num))

# ============================================================
# 4A) NORMALITY TESTING (Informational Only)
# Tests are run for transparency; test selection is NOT adaptive
# Paired t-test is used consistently throughout (see section 5)
# ============================================================
cat("\n=== NORMALITY TESTING (Shapiro-Wilk) — Informational Only ===\n")
cat("Testing paired differences for each group and comparison\n")
cat("(Paired t-test is used throughout regardless of normality result)\n\n")

normality_results <- data.frame()

for (grp in z_order) {
  for (comp in my_comparison_weeks) {
    week1 <- comp[1]
    week2 <- comp[2]
    
    data_wide <- dat %>%
      filter(Group == grp, Week_num %in% c(week1, week2)) %>%
      select(ID, Week_num, FBG) %>%
      pivot_wider(names_from = Week_num, values_from = FBG, names_prefix = "Week_")
    
    if (ncol(data_wide) >= 3) {
      differences <- data_wide[[paste0("Week_", week2)]] - data_wide[[paste0("Week_", week1)]]
      differences <- differences[!is.na(differences)]
      
      if (length(differences) >= 3) {
        shapiro_result <- shapiro.test(differences)
        is_normal      <- shapiro_result$p.value > 0.05
        
        normality_results <- rbind(normality_results, data.frame(
          Group      = grp,
          Comparison = paste0("Week ", week1, " vs ", week2),
          n          = length(differences),
          Shapiro_W  = round(shapiro_result$statistic, 4),
          Shapiro_p  = round(shapiro_result$p.value, 4),
          Normal     = is_normal,
          Note       = "Informational only — paired t-test used throughout"
        ))
        
        cat(sprintf("%s, Week %s vs %s: n=%d, W = %.4f, p = %.4f → %s\n",
                    grp, week1, week2, length(differences),
                    shapiro_result$statistic,
                    shapiro_result$p.value,
                    ifelse(is_normal, "Normal", "Non-normal")))
      }
    }
  }
  cat("\n")
}

# Export normality results
write.csv(
  normality_results,
  file      = "Fig_2a_Normality_Tests.csv",
  row.names = FALSE
)
cat("✓ Normality test results exported to Fig_2a_Normality_Tests.csv\n")

# ============================================================
# 5) PERFORM PAIRED T-TESTS (consistent method, no BH correction)
# ============================================================
cat("\n=== PERFORMING PAIRED T-TESTS (raw p-values, no correction) ===\n")

all_test_results <- data.frame()

for (grp in z_order) {
  for (comp in my_comparison_weeks) {
    week1 <- comp[1]
    week2 <- comp[2]
    
    dat_subset <- dat %>%
      filter(Group == grp, Week_num %in% c(week1, week2))
    
    test_result <- compare_means(
      formula     = FBG ~ Week_num,
      data        = dat_subset,
      method      = "t.test",   # Paired t-test used consistently
      paired      = TRUE,
      comparisons = list(c(week1, week2))
    ) %>%
      mutate(
        Group       = grp,
        Test_Method = "Paired t-test"
      )
    
    all_test_results <- rbind(all_test_results, test_result)
  }
}

# Format results — raw p-values only, no BH correction
all_test_results <- all_test_results %>%
  mutate(p = signif(p, 3)) %>%
  rename(Variable = .y.) %>%
  select(Group, Variable, group1, group2, Test_Method, p, p.signif)

cat("\n=== Statistical Test Results (raw p-values) ===\n")
print(all_test_results)

# Export test results
write.csv(
  all_test_results,
  file      = "Fig_2a_Statistical_Tests.csv",
  row.names = FALSE
)
cat("✓ Statistical test results exported to Fig_2a_Statistical_Tests.csv\n")

# ============================================================
# 6) Create boxplot with within-subject comparisons
# VISUAL APPEARANCE MATCHED TO FIG 1B
# ============================================================
cat("\n=== GENERATING BOXPLOT ===\n")

combined_plot <- ggboxplot(
  dat, x = "Week_num", y = "FBG",
  color  = "Group", fill = "Group",
  alpha  = 0.5, add = "jitter",
  add.params = list(color = "Group", shape = 21,
                    fill = "white", size = 0.75),
  palette    = diabetologia_colors,
  outlier.shape = NA
) +
  scale_y_continuous(limits = c(50, 450), breaks = c(100, 200, 300, 400)) +
  labs(
    y = "Fasting blood glucose (mg/dL)",
    x = "Week"
  ) +
  facet_wrap(vars(Group),
             ncol      = 3,
             labeller  = labeller(Group = Group.labs)) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 8),
    strip.background = element_blank(),
    axis.title       = element_text(size = 8),
    axis.text.x      = element_text(size = 8, colour = "black"),
    axis.text.y      = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
    axis.line        = element_line(colour = "black", linewidth = 0.75 / .pt),
    axis.ticks       = element_line(colour = "black", linewidth = 0.75 / .pt),
    panel.spacing    = unit(0.3, "cm"),
    plot.margin      = margin(3.8, 3.8, 3.8, 3.8, "mm")
  ) +
  stat_compare_means(
    comparisons  = my_comparison_weeks,
    method       = "t.test",   # Paired t-test — consistent with exported CSV
    paired       = TRUE,
    tip.length   = 0.0,
    step.increase = 0.39,
    label        = "p.signif",
    size         = 2.5,
    bracket.size = 0.3
  )

print(combined_plot)

# Save combined plot in PNG and SVG
ggsave(
  filename  = "Fig_2a_Boxplot_FBG.png",
  plot      = combined_plot,
  scale = 1, width = 8.5, height = 5.7,
  units = "cm", dpi = 600, limitsize = FALSE
)

ggsave(
  filename  = "Fig_2a_Boxplot_FBG.svg",
  plot      = combined_plot,
  scale = 1, width = 8.5, height = 5.7,
  units = "cm", limitsize = FALSE
)

cat("✓ Boxplot saved in PNG and SVG formats\n")

# ============================================================
# 7) Summary statistics by Group and Week
# ============================================================
cat("\n=== CALCULATING SUMMARY STATISTICS ===\n")

summary_stats <- dat %>%
  group_by(Group, Week_num) %>%
  summarise(
    n       = n(),
    mean    = round(mean(FBG, na.rm = TRUE), 2),
    sd      = round(sd(FBG, na.rm = TRUE), 2),
    sem     = round(sd(FBG, na.rm = TRUE) / sqrt(n()), 2),
    median  = round(median(FBG, na.rm = TRUE), 2),
    min     = round(min(FBG, na.rm = TRUE), 2),
    max     = round(max(FBG, na.rm = TRUE), 2),
    .groups = 'drop'
  ) %>%
  mutate(Mean_SEM = paste0(mean, " \u00b1 ", sem))

cat("\n=== Summary Statistics (Mean \u00b1 SEM) ===\n")
print(summary_stats)

write.csv(
  summary_stats,
  file      = "Fig_2a_Summary_FBG.csv",
  row.names = FALSE
)
cat("✓ Summary statistics exported to Fig_2a_Summary_FBG.csv\n")

# ============================================================
# 8) Calculate FBG changes (Δ FBG)
# ============================================================
cat("\n=== CALCULATING FBG CHANGES ===\n")

fbg_changes_list <- list()

for (comp in my_comparison_weeks) {
  week1     <- comp[1]
  week2     <- comp[2]
  comp_name <- paste0("Delta_", week1, "_", week2)
  pct_name  <- paste0("Pct_",   week1, "_", week2)
  
  fbg_changes_temp <- dat %>%
    filter(Week_num %in% c(week1, week2)) %>%
    select(ID, Group, Week_num, FBG) %>%
    pivot_wider(names_from = Week_num, values_from = FBG, names_prefix = "Week_") %>%
    mutate(
      !!comp_name := .data[[paste0("Week_", week2)]] - .data[[paste0("Week_", week1)]],
      !!pct_name  := round(((.data[[paste0("Week_", week2)]] - .data[[paste0("Week_", week1)]]) /
                              .data[[paste0("Week_", week1)]]) * 100, 2)
    ) %>%
    select(ID, Group, all_of(c(comp_name, pct_name)))
  
  fbg_changes_list[[paste0(week1, "_", week2)]] <- fbg_changes_temp
}

fbg_changes <- fbg_changes_list[[1]]
for (i in 2:length(fbg_changes_list)) {
  fbg_changes <- fbg_changes %>%
    full_join(fbg_changes_list[[i]], by = c("ID", "Group"))
}

write.csv(
  fbg_changes,
  file      = "Fig_2a_FBG_Changes.csv",
  row.names = FALSE
)
cat("✓ FBG changes exported to Fig_2a_FBG_Changes.csv\n")

# ============================================================
# 9) Summary of FBG changes by Group
# ============================================================
cat("\n=== SUMMARY OF FBG CHANGES BY GROUP ===\n")

fbg_change_summary <- fbg_changes %>%
  group_by(Group) %>%
  summarise(
    n               = n(),
    Mean_Delta_0_6  = round(mean(Delta_0_6,  na.rm = TRUE), 2),
    SEM_Delta_0_6   = round(sd(Delta_0_6,    na.rm = TRUE) / sqrt(n()), 2),
    Mean_Delta_6_12 = round(mean(Delta_6_12, na.rm = TRUE), 2),
    SEM_Delta_6_12  = round(sd(Delta_6_12,   na.rm = TRUE) / sqrt(n()), 2),
    Mean_Delta_0_12 = round(mean(Delta_0_12, na.rm = TRUE), 2),
    SEM_Delta_0_12  = round(sd(Delta_0_12,   na.rm = TRUE) / sqrt(n()), 2),
    .groups = 'drop'
  )

cat("\nSummary of FBG changes:\n")
print(fbg_change_summary)

write.csv(
  fbg_change_summary,
  file      = "Fig_2a_Summary_FBG_Changes.csv",
  row.names = FALSE
)
cat("✓ FBG change summary exported to Fig_2a_Summary_FBG_Changes.csv\n")

# ============================================================
# 10) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat(" 1. Fig_2a_Normality_Tests.csv    (Shapiro-Wilk — informational only)\n")
cat(" 2. Fig_2a_Statistical_Tests.csv  (paired t-tests, raw p-values)\n")
cat(" 3. Fig_2a_Boxplot_FBG.png        (600 dpi)\n")
cat(" 4. Fig_2a_Boxplot_FBG.svg        (vector format)\n")
cat(" 5. Fig_2a_Summary_FBG.csv        (descriptive statistics by week)\n")
cat(" 6. Fig_2a_FBG_Changes.csv        (individual-level Δ FBG)\n")
cat(" 7. Fig_2a_Summary_FBG_Changes.csv (group-level Δ FBG)\n")
cat("\nFigure dimensions: 8.5 cm (width) × 5.7 cm (height)\n")
cat("\nVISUAL APPEARANCE:\n")
cat(" - Matched to Fig 1b for consistency\n")
cat(" - Uses Diabetologia colour palette\n")
cat(" - Boxplots with 50% transparency (alpha = 0.5)\n")
cat(" - Jittered points: white fill, shape 21\n")
cat(" - Statistical brackets: size 0.3\n")
cat(" - Significance labels: size 2.5\n")
cat("\nStatistical approach:\n")
cat(" - Design: PAIRED (within-subject, repeated measures)\n")
cat(" - Normality testing: Shapiro-Wilk on paired differences (informational)\n")
cat(" - Test used throughout: Paired t-test (consistent in plot and CSV)\n")
cat(" - Multiple comparison correction: NONE (raw p-values reported)\n")
cat("\nComparisons per group:\n")
cat(" - Week 0 vs 6\n")
cat(" - Week 6 vs 12\n")
cat(" - Week 0 vs 12\n")
cat("   Total: 9 comparisons (3 groups × 3 time comparisons)\n")
cat("\nSignificance levels (raw p-values):\n")
cat(" *** p<0.001\n")
cat("  ** p<0.01\n")
cat("   * p<0.05\n")
cat("  ns  not significant\n")
cat(rep("=", 60), "\n", sep = "")
