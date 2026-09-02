# ============================================================
# FITC-Dextran Intestinal Permeability Analysis - UNPAIRED
# Revised for Diabetologia Journal Submission
# Unpaired t-test throughout - uniform across panel
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
library(patchwork)

# 2) Define Diabetologia colour palette
diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

# 3) Define comparison weeks by group
my_comparison_weeks_CTRL   <- list(c("0", "12"))
my_comparison_weeks_IL10KO <- list(c("0", "12"))
my_comparison_weeks_DSS    <- list(
  c("0", "1"),
  c("0", "7"),
  c("0", "12"),
  c("7", "12")
)

# ============================================================
# 4) Import FITC-dextran data
# ============================================================
cat("\n=== LOADING FITC-DEXTRAN DATA ===\n")
dat <- read_csv("Fig_1c_Data_FITC_dextran.csv")
dat <- dat %>% filter(Replicate %in% c("Rep1","Rep2","Rep3","Rep4","Rep5"))

dat <- dat %>%
  filter(!is.na(FITC_Conc))

dat$Group    <- factor(dat$Group,
                       levels = c("CTRL", "IL10-KO", "DSS-IL10-KO"))
dat$Week_num <- as.factor(dat$Week_num)

z_order    <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
Group.labs <- c("CTRL", "IL10-KO", "DSS-IL10-KO")
names(Group.labs) <- z_order

cat("\n=== DATA STRUCTURE FOR FITC-DEXTRAN ===\n")
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
  weeks_in_group <- unique(dat$Week_num[dat$Group == grp])
  
  for (week in weeks_in_group) {
    data_subset <- dat %>%
      filter(Group == grp, Week_num == week) %>%
      pull(FITC_Conc)
    
    if (length(data_subset) >= 3) {
      shapiro_result <- shapiro.test(data_subset)
      
      normality_results <- rbind(normality_results, data.frame(
        Group     = grp,
        Week      = week,
        n         = length(data_subset),
        Shapiro_W = round(shapiro_result$statistic, 4),
        Shapiro_p = round(shapiro_result$p.value, 4),
        Normal    = shapiro_result$p.value > 0.05,
        Test_Used = "Unpaired t-test (uniform across panel)"
      ))
      
      cat(sprintf("%s, Week %s: W = %.4f, p = %.4f (%s)\n",
                  grp, week,
                  shapiro_result$statistic,
                  shapiro_result$p.value,
                  ifelse(shapiro_result$p.value > 0.05,
                         "Normal", "Non-normal")))
    }
  }
}

write.csv(normality_results,
          file      = "Fig_1c_Normality_Results.csv",
          row.names = FALSE)
cat("\n✓ Normality results exported (informational only)\n")

# ============================================================
# 5) STATISTICAL TESTS - UNPAIRED T-TEST THROUGHOUT
# ============================================================
cat("\n=== PERFORMING UNPAIRED T-TESTS (uniform across panel) ===\n")

# Build full comparison list with group labels
all_comparisons <- list(
  list(group = "CTRL",        comps = my_comparison_weeks_CTRL),
  list(group = "IL10-KO",     comps = my_comparison_weeks_IL10KO),
  list(group = "DSS-IL10-KO", comps = my_comparison_weeks_DSS)
)

all_test_results <- data.frame()

for (entry in all_comparisons) {
  grp <- entry$group
  for (comp in entry$comps) {
    week1 <- comp[1]
    week2 <- comp[2]
    
    dat_subset <- dat %>%
      filter(Group == grp, Week_num %in% c(week1, week2))
    
    test_result <- compare_means(
      formula     = FITC_Conc ~ Week_num,
      data        = dat_subset,
      method      = "t.test",      # Uniform unpaired t-test throughout
      paired      = FALSE,
      comparisons = list(c(week1, week2))
    ) %>%
      mutate(
        Group       = grp,
        Test_Method = "unpaired t-test"
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

cat("\n=== Statistical Test Results (uncorrected, unpaired t-test) ===\n")
print(all_test_results)

write.csv(all_test_results,
          file      = "Fig_1c_Statistical_Tests.csv",
          row.names = FALSE)
cat("\n✓ Statistical test results exported to Fig_1c_Statistical_Tests.csv\n")

# ============================================================
# 6) PLOTS - UNPAIRED T-TEST THROUGHOUT (consistent with CSV)
# ============================================================
cat("\n=== GENERATING BOXPLOTS ===\n")

dat_CTRL   <- dat %>% filter(Group == "CTRL")
dat_IL10KO <- dat %>% filter(Group == "IL10-KO")
dat_DSS    <- dat %>% filter(Group == "DSS-IL10-KO")

# Shared theme elements
theme_fitc <- theme_classic() %+replace% theme(
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
                                  linewidth = 0.75 / .pt)
)

# Helper function for boxplot base
make_fitc_plot <- function(data) {
  ggboxplot(
    data, x = "Week_num", y = "FITC_Conc",
    color  = "Group", fill = "Group",
    alpha  = 0.5, add = "jitter",
    add.params = list(color = "Group", shape = 21,
                      fill  = "white", size  = 0.75),
    palette       = diabetologia_colors,
    outlier.shape = NA
  ) +
    lims(y = c(-0.0058, 5)) +
    labs(y = "FITC-dextran (µg/mL)", x = "Week") +
    facet_wrap(vars(Group), ncol = 1,
               labeller = labeller(Group = Group.labs))
}

#--------------CTRL plot--------------
plot_CTRL <- make_fitc_plot(dat_CTRL) +
  theme_fitc +
  theme(axis.title.x = element_blank()) +
  stat_compare_means(
    comparisons   = my_comparison_weeks_CTRL,
    method        = "t.test",
    paired        = FALSE,
    tip.length    = 0.0,
    step.increase = 0.3,
    label         = "p.signif",
    size          = 2.5,
    bracket.size  = 0.3
  )

#--------------IL10-KO plot--------------
plot_IL10KO <- make_fitc_plot(dat_IL10KO) +
  theme_fitc +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  stat_compare_means(
    comparisons   = my_comparison_weeks_IL10KO,
    method        = "t.test",
    paired        = FALSE,
    tip.length    = 0.0,
    step.increase = 0.3,
    label         = "p.signif",
    size          = 2.5,
    bracket.size  = 0.3
  )

#--------------DSS-IL10-KO plot--------------
plot_DSS <- make_fitc_plot(dat_DSS) +
  theme_fitc +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  stat_compare_means(
    comparisons   = my_comparison_weeks_DSS,
    method        = "t.test",
    paired        = FALSE,
    tip.length    = 0.0,
    step.increase = 0.3,
    label         = "p.signif",
    size          = 2.5,
    bracket.size  = 0.3
  )

# Combine plots
combined_plot <- plot_CTRL + plot_IL10KO + plot_DSS +
  plot_layout(ncol = 3, widths = c(1.5, 1.5, 1.5)) &
  theme(plot.margin = unit(c(2.8, 2.8, 2.8, 2.8), "pt"))

print(combined_plot)

ggsave(filename = "Fig_1c_Boxplot_FITC_dextran.png",
       plot = combined_plot,
       scale = 1, width = 8.5, height = 5.5,
       units = "cm", dpi = 600, limitsize = FALSE)

ggsave(filename = "Fig_1c_Boxplot_FITC_dextran.svg",
       plot = combined_plot,
       scale = 1, width = 8.5, height = 5.5,
       units = "cm", limitsize = FALSE)

cat("✓ Boxplots saved in PNG and SVG formats\n")

# ============================================================
# 7) Summary statistics
# ============================================================
summary_stats <- dat %>%
  group_by(Group, Week_num) %>%
  summarise(
    n      = n(),
    mean   = round(mean(FITC_Conc, na.rm = TRUE), 2),
    sd     = round(sd(FITC_Conc,   na.rm = TRUE), 2),
    sem    = round(sd(FITC_Conc,   na.rm = TRUE) / sqrt(n()), 2),
    median = round(median(FITC_Conc, na.rm = TRUE), 2),
    min    = round(min(FITC_Conc,  na.rm = TRUE), 2),
    max    = round(max(FITC_Conc,  na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(Mean_SEM = paste0(mean, " ± ", sem))

print(summary_stats)

write.csv(summary_stats,
          file      = "Fig_1c_Summary_FITC_dextran.csv",
          row.names = FALSE)
cat("✓ Summary statistics exported to Fig_1c_Summary_FITC_dextran.csv\n")

# ============================================================
# 8) Summary
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nGenerated files:\n")
cat(" 1. Fig_1c_Normality_Results.csv   (informational only)\n")
cat(" 2. Fig_1c_Statistical_Tests.csv   (unpaired t-test, uncorrected p)\n")
cat(" 3. Fig_1c_Boxplot_FITC_dextran.png (600 dpi)\n")
cat(" 4. Fig_1c_Boxplot_FITC_dextran.svg (vector format)\n")
cat(" 5. Fig_1c_Summary_FITC_dextran.csv (descriptive statistics)\n")
cat("\nStatistical approach:\n")
cat(" - Design: UNPAIRED (independent samples at each timepoint)\n")
cat(" - Unpaired t-test throughout (uniform across panel)\n")
cat(" - Uncorrected p values reported\n")
cat(" - Normality tests retained for transparency (informational)\n")
cat("\nComparisons:\n")
cat(" - CTRL:        Week 0 vs 12\n")
cat(" - IL10-KO:     Week 0 vs 12\n")
cat(" - DSS-IL10-KO: Week 0 vs 1, 0 vs 7, 0 vs 12, 7 vs 12\n")
cat(rep("=", 60), "\n", sep = "")
