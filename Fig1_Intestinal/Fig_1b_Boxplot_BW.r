# ============================================================
# Body Weight Analysis - Within-Subject Comparisons
# Revised for Diabetologia Journal Submission
# Paired t-test throughout - uniform across panel
# Comparing Week 0 vs Week 1, and Week 6 vs Week 7
# Three groups: CTRL, IL10-KO, DSS-IL10-KO
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

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define comparison weeks (post-DSS treatment)
my_comparison_weeks <- list(
  c("0", "1"),  # First DSS treatment effect
  c("6", "7")   # Second DSS treatment effect
)

# ============================================================
# 4) Import body weight data
# ============================================================
cat("\n=== LOADING BODY WEIGHT DATA ===\n")
dat <- read_csv("Fig_1b_Data_BW.csv")
dat <- dat %>% filter(Replicate %in% c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5"))

dat <- dat %>%
  filter(Week_num %in% c(0, 1, 6, 7)) %>%
  filter(!is.na(BW))

dat$ID       <- as.factor(dat$ID)
dat$Group    <- factor(dat$Group,
                       levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))
dat$Week_num <- as.factor(dat$Week_num)

z_order    <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
Group.labs <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
names(Group.labs) <- z_order

cat("\n=== DATA STRUCTURE FOR BODY WEIGHT ===\n")
print(head(dat, 20))
cat("\nTotal observations:", nrow(dat), "\n")
cat("\nSample sizes per group per week:\n")
print(table(dat$Group, dat$Week_num))

# ============================================================
# 4A) NORMALITY TESTING (informational only)
# ============================================================
cat("\n=== NORMALITY TESTING (Shapiro-Wilk, informational only) ===\n")

normality_results <- data.frame()

for (grp in z_order) {
  for (comp in my_comparison_weeks) {
    week1 <- comp[1]
    week2 <- comp[2]
    
    data_wide <- dat %>%
      filter(Group == grp, Week_num %in% c(week1, week2)) %>%
      select(ID, Week_num, BW) %>%
      pivot_wider(names_from   = Week_num,
                  values_from  = BW,
                  names_prefix = "Week_")
    
    if (ncol(data_wide) >= 3) {
      differences <- data_wide[[paste0("Week_", week2)]] -
        data_wide[[paste0("Week_", week1)]]
      differences <- differences[!is.na(differences)]
      
      if (length(differences) >= 3) {
        shapiro_result <- shapiro.test(differences)
        
        normality_results <- rbind(normality_results, data.frame(
          Group      = grp,
          Comparison = paste0("Week ", week1, " vs ", week2),
          n          = length(differences),
          Shapiro_W  = round(shapiro_result$statistic, 4),
          Shapiro_p  = round(shapiro_result$p.value, 4),
          Normal     = shapiro_result$p.value > 0.05,
          Test_Used  = "Paired t-test (uniform across panel)"
        ))
        
        cat(sprintf("%s, Week %s vs %s: W = %.4f, p = %.4f (%s)\n",
                    grp, week1, week2,
                    shapiro_result$statistic,
                    shapiro_result$p.value,
                    ifelse(shapiro_result$p.value > 0.05,
                           "Normal", "Non-normal")))
      }
    }
  }
}

write.csv(normality_results,
          file      = "Fig_1b_Normality_Tests.csv",
          row.names = FALSE)
cat("\n✓ Normality results exported (informational only)\n")

# ============================================================
# 5) STATISTICAL TESTS - PAIRED T-TEST THROUGHOUT
# ============================================================
cat("\n=== PERFORMING PAIRED T-TESTS (uniform across panel) ===\n")

all_test_results <- data.frame()

for (grp in z_order) {
  for (comp in my_comparison_weeks) {
    week1 <- comp[1]
    week2 <- comp[2]
    
    dat_subset <- dat %>%
      filter(Group == grp, Week_num %in% c(week1, week2))
    
    test_result <- compare_means(
      formula     = BW ~ Week_num,
      data        = dat_subset,
      method      = "t.test",      # Uniform paired t-test throughout
      paired      = TRUE,
      comparisons = list(c(week1, week2))
    ) %>%
      mutate(
        Group       = grp,
        Test_Method = "paired t-test"
      )
    
    all_test_results <- rbind(all_test_results, test_result)
  }
}

# Significance labels from uncorrected p values
all_test_results <- all_test_results %>%
  mutate(
    p        = signif(p, 3),
    p.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  ) %>%
  rename(Variable = .y.) %>%
  select(Group, Variable, group1, group2, Test_Method, p, p.signif)

cat("\n=== Statistical Test Results (uncorrected, paired t-test) ===\n")
print(all_test_results)

write.csv(all_test_results,
          file      = "Fig_1b_Statistical_Tests.csv",
          row.names = FALSE)
cat("\n✓ Statistical test results exported to Fig_1b_Statistical_Tests.csv\n")

# ============================================================
# 6) PLOT - PAIRED T-TEST THROUGHOUT (consistent with CSV)
# ============================================================
cat("\n=== GENERATING BOXPLOT ===\n")

combined_plot <- ggboxplot(
  dat, x = "Week_num", y = "BW",
  color  = "Group", fill = "Group",
  alpha  = 0.5, add = "jitter",
  add.params = list(color = "Group", shape = 21,
                    fill  = "white", size  = 0.75),
  palette       = diabetologia_colors,
  outlier.shape = NA
) +
  lims(y = c(20, 35)) +
  labs(y = "Body weight (g)", x = "Week") +
  facet_wrap(vars(Group), ncol = 3,
             labeller = labeller(Group = Group.labs)) +
  theme_classic() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 8),
    strip.background = element_blank(),
    axis.title       = element_text(size = 8),
    axis.text.x      = element_text(size = 8, colour = "black"),
    axis.text.y      = element_text(size = 8, colour = "black",
                                    angle = 90, hjust = 0.5),
    axis.line        = element_line(colour = "black",
                                    linewidth = 0.75 / .pt),
    axis.ticks       = element_line(colour = "black",
                                    linewidth = 0.75 / .pt),
    #panel.spacing    = unit(0.3, "cm"),
    #plot.margin      = margin(3.8, 3.8, 3.8, 3.8, "mm")
  ) +
  stat_compare_means(
    comparisons   = my_comparison_weeks,
    method        = "t.test",    # Consistent with CSV export
    paired        = TRUE,
    tip.length    = 0.0,
    step.increase = 0.22,
    label         = "p.signif", # Uncorrected p values
    size          = 2.5,
    bracket.size  = 0.3
  )

print(combined_plot)

ggsave(filename = "Fig_1b_Boxplot_BW.png",
       plot = combined_plot,
       scale = 1, width = 8.5, height = 5.5,
       units = "cm", dpi = 600, limitsize = FALSE)

ggsave(filename = "Fig_1b_Boxplot_BW.svg",
       plot = combined_plot,
       scale = 1, width = 8.5, height = 5.5,
       units = "cm", limitsize = FALSE)

cat("✓ Boxplot saved in PNG and SVG formats\n")

# ============================================================
# 7) Summary statistics
# ============================================================
summary_stats <- dat %>%
  group_by(Group, Week_num) %>%
  summarise(
    n      = n(),
    mean   = round(mean(BW, na.rm = TRUE), 2),
    sd     = round(sd(BW,   na.rm = TRUE), 2),
    sem    = round(sd(BW,   na.rm = TRUE) / sqrt(n()), 2),
    median = round(median(BW, na.rm = TRUE), 2),
    min    = round(min(BW,  na.rm = TRUE), 2),
    max    = round(max(BW,  na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(Mean_SEM = paste0(mean, " ± ", sem))

print(summary_stats)

write.csv(summary_stats,
          file      = "Fig_1b_Summary_BW.csv",
          row.names = FALSE)
cat("✓ Summary statistics exported to Fig_1b_Summary_BW.csv\n")

# ============================================================
# 8) Body weight changes
# ============================================================
bw_changes <- dat %>%
  select(ID, Group, Week_num, BW) %>%
  pivot_wider(names_from   = Week_num,
              values_from  = BW,
              names_prefix = "Week_") %>%
  mutate(
    Delta_BW_0_1       = Week_1 - Week_0,
    Delta_BW_6_7       = Week_7 - Week_6,
    PercentChange_0_1  = round((Delta_BW_0_1 / Week_0) * 100, 2),
    PercentChange_6_7  = round((Delta_BW_6_7 / Week_6) * 100, 2)
  )

write.csv(bw_changes,
          file      = "Fig_1b_BW_Changes.csv",
          row.names = FALSE)
cat("✓ Body weight changes exported to Fig_1b_BW_Changes.csv\n")

# ============================================================
# 9) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat(" 1. Fig_1b_Normality_Tests.csv   (informational only)\n")
cat(" 2. Fig_1b_Statistical_Tests.csv (paired t-test, uncorrected p)\n")
cat(" 3. Fig_1b_Boxplot_BW.png        (600 dpi)\n")
cat(" 4. Fig_1b_Boxplot_BW.svg        (vector format)\n")
cat(" 5. Fig_1b_Summary_BW.csv        (descriptive statistics)\n")
cat(" 6. Fig_1b_BW_Changes.csv        (individual delta values)\n")
cat("\nStatistical approach:\n")
cat(" - Paired t-test throughout (uniform across panel)\n")
cat(" - Uncorrected p values reported\n")
cat(" - Normality tests retained for transparency (informational)\n")
cat(" - Comparisons: Week 0 vs 1 (first DSS effect),")
cat("   Week 6 vs 7 (second DSS effect)\n")
cat(rep("=", 60), "\n", sep = "")
