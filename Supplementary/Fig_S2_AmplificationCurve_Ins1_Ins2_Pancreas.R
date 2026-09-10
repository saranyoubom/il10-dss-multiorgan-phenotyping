# ============================================================
# Supplementary Figure 2 - RT-qPCR amplification curves for
# Ins1 and Ins2 in pancreas, all three groups
#
# Purpose: direct instrument-level evidence that Ins1 shows no
# amplification (RFU stays at baseline through all 40 cycles)
# in IL10-KO pancreas, consistent with the Cq being reported
# as Undetermined by the instrument and with "undetectable"
# throughout the manuscript. Ins2 amplifies normally in all
# three groups on the same qPCR plates.
#
# Source: raw CFX Maestro baseline-subtracted amplification
# (RFU vs Cycle) export, four biological replicates per
# gene x group, extracted directly from the original
# instrument output files (not re-analysed or smoothed).
#
# Related to: Discussion, pancreatic beta-cell paragraph
# ("Ins1 was undetectable in IL10-KO mice"); R2 Comment 7's
# request for technical/biological reasons for Ins1 absence.
# ============================================================

rm(list = ls())
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

library(dplyr)
library(readr)
library(ggplot2)

diabetologia_colors <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)
)

dat <- read_csv("Fig_S2_AmplificationCurve_Ins1_Ins2_Pancreas_Data.csv", show_col_types = FALSE) %>%
  mutate(
    Group = factor(Group, levels = c("CTRL", "IL10-KO", "DSS-IL10-KO")),
    Gene  = factor(Gene, levels = c("Ins1", "Ins2")),
    Curve = paste(Gene, Group, Replicate, sep = "_")
  )

cat("Loaded", nrow(dat), "rows across", n_distinct(dat$Curve), "individual curves\n")

p <- ggplot(dat, aes(x = Cycle, y = RFU, color = Group, group = Curve)) +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  scale_color_manual(values = diabetologia_colors) +
  scale_x_continuous(breaks = seq(0, 40, 10)) +
  labs(
    x = "Cycle",
    y = "RFU (baseline-subtracted)",
    color = "Group",
    title = "Amplification curves: pancreatic Ins1 and Ins2"
  ) +
  theme_bw(base_family = "Arial", base_size = 11) +
  theme(
    strip.text = element_text(face = "italic", size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave("Fig_S2_AmplificationCurve_Ins1_Ins2_Pancreas.png", p, width = 6.2, height = 7.5, dpi = 300)
ggsave("Fig_S2_AmplificationCurve_Ins1_Ins2_Pancreas.svg", p, width = 6.2, height = 7.5)

# Cq determination summary (from the instrument's own Quantification Cq
# Results export; "Undetermined" = no threshold crossing within 40 cycles)
cq_status <- tribble(
  ~Group,         ~Gene, ~Cq_Determined, ~Total,
  "CTRL",         "Ins1", 2, 4,
  "IL10-KO",      "Ins1", 0, 4,
  "DSS-IL10-KO",  "Ins1", 3, 4,
  "CTRL",         "Ins2", 4, 4,
  "IL10-KO",      "Ins2", 4, 4,
  "DSS-IL10-KO",  "Ins2", 4, 4
)
write_csv(cq_status, "Fig_S2_AmplificationCurve_CqDetermination_Summary.csv")

cat("\nSaved:\n")
cat(" - Fig_S2_AmplificationCurve_Ins1_Ins2_Pancreas.png\n")
cat(" - Fig_S2_AmplificationCurve_Ins1_Ins2_Pancreas.svg\n")
cat(" - Fig_S2_AmplificationCurve_CqDetermination_Summary.csv\n")

# ============================================================
# END OF SCRIPT
# ============================================================
