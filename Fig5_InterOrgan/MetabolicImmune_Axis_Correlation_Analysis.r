# ==============================================================================
# Inter-Organ Correlation Analysis - GUT-PANCREAS-LIVER-KIDNEY-SPLEEN AXIS
# DSS-induced IL10-knockout mice
# Diabetologia colour theme - Metabolic and immune organs
# ==============================================================================

rm(list = ls())

# Set working directory
# Set the working directory to this script's own folder before running
# (e.g. in RStudio: Session > Set Working Directory > To Source File Location)

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

library(tidyverse)   # Data manipulation and visualisation
library(Hmisc)       # Correlation with p-values
library(ggpubr)      # stat_cor
library(reshape2)    # Data reshaping
library(ggraph)      # Network graphs
library(igraph)      # Network analysis
library(svglite)     # SVG export
library(ggrepel)     # Text labels
library(gridExtra)   # Multiple plots

cat("=============================================================\n")
cat("Gut-Pancreas-Liver-Kidney-Spleen Axis Correlation Analysis\n")
cat("Diabetologia Format\n")
cat("=============================================================\n\n")

# ==============================================================================
# OFFICIAL DIABETOLOGIA COLOUR PALETTE
# ==============================================================================

# Experimental group colours
diabetologia_groups <- c(
  "CTRL"        = rgb(120, 131, 186, maxColorValue = 255),  # Blue
  "IL10-KO"     = rgb(204, 133, 177, maxColorValue = 255),  # Mauve
  "DSS-IL10-KO" = rgb(237, 115, 116, maxColorValue = 255)   # Red
)

# Metabolic and immune organ colours (5 organs)
diabetologia_organs <- c(
  "Kidney"           = rgb(120, 131, 186, maxColorValue = 255),  # Blue
  "Colon"  = rgb(152, 201, 139, maxColorValue = 255),  # Green
  "Colon"            = rgb(152, 201, 139, maxColorValue = 255),  # Green (renamed from Colon)
  "Liver"            = rgb(242, 197, 167, maxColorValue = 255),  # Peach
  "Pancreas"         = rgb(237, 115, 116, maxColorValue = 255),  # Red
  "Spleen"           = rgb(204, 133, 177, maxColorValue = 255)   # Mauve
)

# Correlation strength colours
correlation_strength_colours <- c(
  "Strong"   = rgb(237, 115, 116, maxColorValue = 255),  # Strong red
  "Moderate" = rgb(242, 197, 167, maxColorValue = 255),  # Orange/peach
  "Weak"     = rgb(251, 241, 173, maxColorValue = 255)   # Light yellow
)

# Significance colours
significance_colours <- c(
  "p<0.001" = rgb(193, 88, 88, maxColorValue = 255),   # Dark red
  "p<0.01"  = rgb(237, 115, 116, maxColorValue = 255), # Red
  "p<0.05"   = rgb(244, 166, 151, maxColorValue = 255), # Light red
  "p>0.05"  = rgb(214, 214, 214, maxColorValue = 255)  # Grey
)

# Correlation direction colours
correlation_direction <- c(
  "positive" = rgb(237, 115, 116, maxColorValue = 255),  # Red
  "negative" = rgb(120, 131, 186, maxColorValue = 255)   # Blue
)

# ==============================================================================
# DIABETOLOGIA THEME
# ==============================================================================

theme_diabetologia <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "right",
      strip.text.x = element_text(size = base_size, face = "plain"),
      strip.background = element_blank(),
      axis.title.y = element_text(size = base_size, colour = "black"),
      axis.title.x = element_text(size = base_size, colour = "black"),
      axis.text.x = element_text(size = base_size, colour = "black"),
      axis.text.y = element_text(size = base_size, colour = "black"),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      axis.line = element_line(colour = "black", linewidth = 0.75 / .pt),
      axis.ticks = element_line(colour = "black", linewidth = 0.75 / .pt),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

theme_set(theme_diabetologia())

# Helper function to save plots
save_plot <- function(plot_obj, filename_base, width = 8, height = 6, dpi = 600) {
  ggsave(paste0(filename_base, ".png"), plot_obj,
         width = width, height = height, units = "cm", dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".svg"), plot_obj,
         width = width, height = height, units = "cm", device = "svg", bg = "white")
  cat(" ✓", filename_base, ".png and .svg\n")
}

# ==============================================================================
# 1. DATA LOADING AND PREPARATION - METABOLIC & IMMUNE ORGANS
# ==============================================================================

cat("Step 1: Loading metabolic and immune organ data...\n")

# Load gut, pancreas, liver, kidney, and spleen data
data_files <- list(
  Pancreas = "Fig_2ef_Pancreas_FoldChange_Data.csv",
  Kidney = "Fig_5b_Data_Kidney_GeneExpr.csv",
  Liver = "Fig_4b_Data_Liver_GeneExpr.csv",
  Colon = "Fig_1f_Data_Colon_GeneExpr.csv",
  Spleen = "Fig_3c_Data_Spleen_GeneExpr.csv"
)

loaded_data <- list()
for (organ_name in names(data_files)) {
  file_path <- data_files[[organ_name]]
  if (file.exists(file_path)) {
    tryCatch({
      loaded_data[[organ_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
      cat(paste(" ✓", organ_name, "\n"))
    }, error = function(e) {
      cat(paste(" ✗", organ_name, "- Error\n"))
    })
  } else {
    cat(paste(" ✗", organ_name, "- File not found\n"))
  }
}

cat(paste("\nOrgans loaded:", length(loaded_data), "\n\n"))

# Standardise data function
standardise_data <- function(df, organ_name) {
  if ("Tissue" %in% colnames(df)) {
    df <- df %>% rename(Organ = Tissue)
  } else {
    df$Organ <- organ_name
  }

  if ("Gene" %in% colnames(df)) {
    df$Gene <- gsub("SRR", "Srr", df$Gene)
  }

  df <- df %>% select(Organ, Gene, Group, Relative_expression)
  return(df)
}

# Apply standardisation
all_organs_list <- lapply(names(loaded_data), function(organ_name) {
  standardise_data(loaded_data[[organ_name]], organ_name)
})
all_organs_data <- bind_rows(all_organs_list)

# Create sample IDs
all_organs_data <- all_organs_data %>%
  group_by(Organ, Group, Gene) %>%
  mutate(Sample_ID = row_number()) %>%
  ungroup() %>%
  mutate(Unique_Sample_ID = paste(Group, Sample_ID, sep = "_"))

cat("Step 2: Data preparation complete\n")
cat(paste(" Organs:", paste(unique(all_organs_data$Organ), collapse = ", "), "\n\n"))

# ==============================================================================
# 2. CREATE WIDE FORMAT FOR CORRELATION
# ==============================================================================

cat("Step 3: Creating wide format...\n")

wide_data <- all_organs_data %>%
  mutate(Organ_Gene = paste(Organ, Gene, sep = "_")) %>%
  select(Unique_Sample_ID, Group, Organ_Gene, Relative_expression) %>%
  pivot_wider(
    names_from = Organ_Gene,
    values_from = Relative_expression,
    values_fn = mean
  )

cat(paste(" Samples:", nrow(wide_data), "\n\n"))

# ==============================================================================
# 3. CORRELATION ANALYSIS
# ==============================================================================

cat("Step 4: Calculating correlations...\n")

calc_pairwise_correlations <- function(data_matrix) {
  numeric_cols <- data_matrix %>%
    select(where(is.numeric)) %>%
    as.data.frame()

  n_vars <- ncol(numeric_cols)
  var_names <- colnames(numeric_cols)

  cor_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
  p_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
  n_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)

  rownames(cor_matrix) <- colnames(cor_matrix) <- var_names
  rownames(p_matrix) <- colnames(p_matrix) <- var_names
  rownames(n_matrix) <- colnames(n_matrix) <- var_names

  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      complete_idx <- complete.cases(numeric_cols[, c(i, j)])
      x <- as.numeric(numeric_cols[complete_idx, i])
      y <- as.numeric(numeric_cols[complete_idx, j])
      n_obs <- sum(complete_idx)

      if (n_obs >= 3) {
        if (i == j) {
          cor_matrix[i, j] <- 1
          p_matrix[i, j] <- 0
          n_matrix[i, j] <- n_obs
        } else {
          cor_test <- cor.test(x, y, method = "pearson")
          cor_matrix[i, j] <- cor_test$estimate
          p_matrix[i, j] <- cor_test$p.value
          n_matrix[i, j] <- n_obs
        }
      }
    }
  }

  return(list(r = cor_matrix, p = p_matrix, n = n_matrix))
}

extract_inter_organ_correlations <- function(cor_matrix, p_matrix, n_matrix) {
  all_col_names <- colnames(cor_matrix)
  organ_names <- unique(sapply(strsplit(all_col_names, "_"), function(x) x[1]))

  results_list <- list()

  for (i in 1:(length(organ_names)-1)) {
    for (j in (i+1):length(organ_names)) {
      organ1 <- organ_names[i]
      organ2 <- organ_names[j]

      organ1_cols <- grep(paste0("^", organ1, "_"), colnames(cor_matrix), value = TRUE)
      organ2_cols <- grep(paste0("^", organ2, "_"), colnames(cor_matrix), value = TRUE)

      if (length(organ1_cols) > 0 & length(organ2_cols) > 0) {
        cor_sub <- cor_matrix[organ1_cols, organ2_cols, drop = FALSE]
        p_sub <- p_matrix[organ1_cols, organ2_cols, drop = FALSE]
        n_sub <- n_matrix[organ1_cols, organ2_cols, drop = FALSE]

        cor_long <- melt(cor_sub)
        p_long <- melt(p_sub)
        n_long <- melt(n_sub)

        result_df <- data.frame(
          Organ1 = organ1,
          Organ2 = organ2,
          Gene_Organ1 = cor_long$Var1,
          Gene_Organ2 = cor_long$Var2,
          Correlation = cor_long$value,
          P_value = p_long$value,
          N = n_long$value
        )

        result_df$Gene1 <- gsub(paste0("^", organ1, "_"), "", result_df$Gene_Organ1)
        result_df$Gene2 <- gsub(paste0("^", organ2, "_"), "", result_df$Gene_Organ2)

        results_list[[paste(organ1, organ2, sep = "_vs_")]] <- result_df
      }
    }
  }

  return(bind_rows(results_list, .id = "Organ_Pair"))
}

# Run correlation analysis
all_groups_cor <- calc_pairwise_correlations(wide_data)
all_groups_inter_organ <- extract_inter_organ_correlations(
  all_groups_cor$r,
  all_groups_cor$p,
  all_groups_cor$n
)

all_groups_inter_organ$Group <- "All_Groups"
all_groups_inter_organ$Analysis_Note <- "Metabolic-immune axis (n=12)"

# Add significance labels
all_inter_organ_correlations <- all_groups_inter_organ %>%
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "p<0.001",
      P_value < 0.01 ~ "p<0.01",
      P_value < 0.05 ~ "p<0.05",
      TRUE ~ "p>0.05"
    ),
    Significant = P_value < 0.05
  )

cat(paste(" Total correlations:", nrow(all_inter_organ_correlations), "\n\n"))

# ==============================================================================
# 4. EXPORT CSV RESULTS
# ==============================================================================

cat("Step 5: Exporting CSV files...\n")

write.csv(all_inter_organ_correlations,
          "MetabolicImmune_Axis_Correlation_Results.csv",
          row.names = FALSE)

significant_correlations <- all_inter_organ_correlations %>%
  filter(Significant == TRUE) %>%
  arrange(P_value)

write.csv(significant_correlations,
          "MetabolicImmune_Axis_Significant_Correlations.csv",
          row.names = FALSE)

cat(" ✓ MetabolicImmune_Axis_Correlation_Results.csv\n")
cat(" ✓ MetabolicImmune_Axis_Significant_Correlations.csv\n\n")

# ==============================================================================
# 5. VISUALISATIONS - METABOLIC-IMMUNE AXIS FOCUS
# ==============================================================================

cat("Step 6: Generating visualisations...\n\n")

# --- Fig 7a: BAR PLOT SUMMARY (REVISED) ---
cat(" Creating Fig 7a: Bar plot summary...\n")

summary_data <- significant_correlations %>%
  group_by(Organ1, Organ2) %>%
  summarise(
    Total = n(),
    Strong = sum(abs(Correlation) > 0.7),
    Moderate = sum(abs(Correlation) > 0.4 & abs(Correlation) <= 0.7),
    #####Weak = sum(abs(Correlation) <= 0.4),
    .groups = "drop"
  ) %>%
  mutate(
    Organ1 = gsub("Colon", "Colon", Organ1),
    Organ2 = gsub("Colon", "Colon", Organ2),
    Organ_Pair = paste(Organ1, "vs\n", Organ2)
  ) %>%
  pivot_longer(cols = c(Strong, Moderate
                        #, Weak
                        ),
               names_to = "Strength",
               values_to = "Count") %>%
  mutate(Strength = factor(Strength, levels = c("Strong", "Moderate"
                                                #, "Weak"
                                                )))

p_barplot <- ggplot(summary_data, aes(x = reorder(Organ_Pair, Total),
                                       y = Count,
                                       fill = Strength)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = ifelse(Count > 0, Count, "")),
            position = position_stack(vjust = 0.5),
            size = 7 / .pt,
            colour = "black") +
  scale_fill_manual(values = correlation_strength_colours,
                    name = "Correlation\nStrength",
                    labels = c("Strong\n(|r|>0.7)",
                              "Moderate\n(|r|=0.4-0.7)",
                              "Weak\n(|r|<0.4)")) +
  coord_flip() +
  labs(
    #title = "Metabolic-Immune Axis: Correlations by Organ Pair",
    #subtitle = "All groups combined (n=12, p<0.05)",
    x = "Organ Pair",
    y = "Number of Significant Correlations"
  ) +
  theme_diabetologia() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8, colour = "black", face = "plain", hjust = 1),
        legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.4, "cm"),
        axis.title.y = element_blank())

save_plot(p_barplot, "Fig_5a_MetabolicImmune_Barplot_Summary", 
          width = 8, height = 8)

# --- Fig 7b: LOLLIPOP CHART (REVISED) ---
cat(" Creating Fig 7b: Lollipop chart...\n")

top_correlations <- significant_correlations %>%
  mutate(Abs_Correlation = abs(Correlation)) %>%
  arrange(desc(Abs_Correlation)) %>%
  head(30) %>%
  mutate(
    Gene_Pair = paste(Gene1, Gene2, sep = "-"),
    Organ1 = gsub("Colon", "Colon", Organ1),
    Organ2 = gsub("Colon", "Colon", Organ2),
    Organ_Pair_Label = paste(Organ1, "vs", Organ2)
  )

if (nrow(top_correlations) > 0) {
  p_lollipop <- ggplot(top_correlations, aes(x = reorder(Gene_Pair, Abs_Correlation),
                                               y = Correlation,
                                               colour = Significance)) +
    geom_segment(aes(xend = Gene_Pair, y = 0, yend = Correlation),
                 linewidth = 0.75 / .pt) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "solid",
               colour = "black", linewidth = 0.75 / .pt) +
    scale_colour_manual(values = significance_colours,
                        name = "Significance\nLevel") +
    coord_flip() +
    facet_wrap(~Organ_Pair_Label, scales = "free_y", ncol = 2) +
    labs(
      #title = "Top 30 Metabolic-Immune Axis Gene Correlations",
      #subtitle = "Gut-Pancreas-Liver-Kidney-Spleen axis (n=12)",
      x = "Gene Pair",
      y = "Pearson Correlation Coefficient (r)"
    ) +
    theme_diabetologia() +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_text(size = 8, face = "plain", 
                                  hjust = 1, vjust = 0.5, color = "black"),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 8, face = "bold")
    )

  save_plot(p_lollipop, "Fig_5b_MetabolicImmune_Lollipop_Top30", 
            width = 9, height = 14)
}

# --- Fig 7c: MULTI-PANEL SCATTER PLOTS (REVISED) ---
cat(" Creating Fig 7c: Multi-panel scatter plots...\n")

plot_gene_pair_correlation <- function(organ1_data, organ2_data,
                                        organ1_name, organ2_name,
                                        gene1, gene2,
                                        axis_title,
                                        group_filter = NULL) {
  df1 <- organ1_data %>%
    filter(Gene == gene1) %>%
    group_by(Group) %>%
    mutate(Sample_ID = row_number()) %>%
    ungroup() %>%
    select(Group, Sample_ID, Relative_expression) %>%
    rename(!!paste0(organ1_name, "_", gene1) := Relative_expression)

  df2 <- organ2_data %>%
    filter(Gene == gene2) %>%
    group_by(Group) %>%
    mutate(Sample_ID = row_number()) %>%
    ungroup() %>%
    select(Group, Sample_ID, Relative_expression) %>%
    rename(!!paste0(organ2_name, "_", gene2) := Relative_expression)

  plot_data <- inner_join(df1, df2, by = c("Group", "Sample_ID"))

  if (!is.null(group_filter)) {
    plot_data <- plot_data %>% filter(Group == group_filter)
  }

  plot_data <- plot_data %>% drop_na()

  col1 <- paste0(organ1_name, "_", gene1)
  col2 <- paste0(organ2_name, "_", gene2)

  if (nrow(plot_data) < 3) {
    return(NULL)
  }

  # Display organ name as "Colon" instead of "Colon"
  organ1_display <- ifelse(organ1_name == "Colon", "Colon", organ1_name)
  organ2_display <- ifelse(organ2_name == "Colon", "Colon", organ2_name)

  p <- ggplot(plot_data, aes_string(x = col1, y = col2, colour = "Group")) +
    geom_point(size = 1.5, alpha = 0.7, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey50",
                fill = "grey80", linetype = "dashed", linewidth = 0.75 / .pt) +
    scale_colour_manual(values = diabetologia_groups) +
    stat_cor(method = "pearson",
             label.x.npc = "left",
             label.y.npc = "top",
             size = 6 / .pt,
             colour = "black") +
    labs(
      title = axis_title,
      x = paste(organ1_display, gene1),
      y = paste(organ2_display, gene2)
    ) +
    theme_diabetologia() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5),
    )+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))

  return(p)
}

# Create 6 scatter plots
tryCatch({
  p1 <- plot_gene_pair_correlation(
    loaded_data$Colon, loaded_data$Spleen,
    "Colon", "Spleen",
    "Zfand6", "Ghrhr",
    "\nGut-Immune Axis"
  )

  p2 <- plot_gene_pair_correlation(
    loaded_data$Pancreas, loaded_data$Kidney,
    "Pancreas", "Kidney",
    "Ghrhr", "Glp1r",
    "\nEndocrine-Renal Axis"
  )

  p3 <- plot_gene_pair_correlation(
    loaded_data$Kidney, loaded_data$Colon,
    "Kidney", "Colon",
    "Mafa", "Tbc1d4",
    "\nRenal-Gut Antagonism"
  )

  p4 <- plot_gene_pair_correlation(
    loaded_data$Pancreas, loaded_data$Colon,
    "Pancreas", "Colon",
    "Zfand6", "Mafa",
    "\nPancreas-Gut Antagonism"
  )

  p5 <- plot_gene_pair_correlation(
    loaded_data$Liver, loaded_data$Spleen,
    "Liver", "Spleen",
    "Tbc1d4", "Tcf19",
    "Metabolic-Immune\nCoordination"
  )

  p6 <- plot_gene_pair_correlation(
    loaded_data$Kidney, loaded_data$Spleen,
    "Kidney", "Spleen",
    "Mafa", "Ghrhr",
    "\nRenal-Immune Axis"
  )

  # Combine plots in 2 columns x 3 rows
  combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)

  # Save combined plot
  ggsave("Fig_5c_MetabolicImmune_Scatter_MultiPanel.png", combined_plot,
         width = 8, height = 15, units = "cm", dpi = 600, bg = "white")
  ggsave("Fig_5c_MetabolicImmune_Scatter_MultiPanel.svg", combined_plot,
         width = 8, height = 15, units = "cm", device = "svg", bg = "white")
  cat(" ✓ Fig_5c_MetabolicImmune_Scatter_MultiPanel.png and .svg\n")

}, error = function(e) {
  cat(" ⚠ Could not create scatter plots:\n")
  cat("  ", conditionMessage(e), "\n")
})

# --- Fig 7d: NETWORK DIAGRAM (REVISED) ---
cat(" Creating Fig 7d: Network diagram...\n")

strong_cors <- significant_correlations %>%
  filter(abs(Correlation) >= 0.7, P_value < 0.01)

if (nrow(strong_cors) > 0) {
  # Replace Colon with Colon
  strong_cors <- strong_cors %>%
    mutate(
      Organ1 = gsub("Colon", "Colon", Organ1),
      Organ2 = gsub("Colon", "Colon", Organ2)
    )

  nodes <- data.frame(
    name = unique(c(
      paste(strong_cors$Organ1, strong_cors$Gene1, sep = "_"),
      paste(strong_cors$Organ2, strong_cors$Gene2, sep = "_")
    ))
  ) %>%
    mutate(
      Organ = sub("_.*", "", name),
      Gene = sub(".*_", "", name)
    )

  edges <- strong_cors %>%
    mutate(
      from = paste(Organ1, Gene1, sep = "_"),
      to = paste(Organ2, Gene2, sep = "_"),
      weight = abs(Correlation),
      sign = ifelse(Correlation > 0, "positive", "negative")
    ) %>%
    select(from, to, weight, sign, Correlation, P_value)

  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  layout <- create_layout(g, layout = "fr")

  p_network <- ggraph(layout) +
    geom_edge_link(aes(width = weight, colour = sign), alpha = 0.6) +
    geom_node_point(aes(colour = Organ), size = 3) +
    geom_node_text(aes(label = Gene), size = 7 / .pt, repel = TRUE, 
                   #, box.padding = 0.5
                   ) +
    scale_edge_width_continuous(name = "|r|", range = c(0.2, 2)) +
    scale_edge_colour_manual(values = correlation_direction,
                             name = "Correlation") +
    scale_colour_manual(values = diabetologia_organs, name = "Organ") +
    labs(
      #title = "Metabolic-Immune Axis Gene Correlation Network",
      #subtitle = "Strong correlations (|r|≥0.7, p<0.01)"
    ) +
    theme_graph() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      legend.spacing.y = unit(0, 'cm'),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, face = "plain")
    )

  save_plot(p_network, "Fig_5d_MetabolicImmune_Network_Strong_Correlations",
            width = 11, height = 9)
}

# --- VOLCANO PLOT (unchanged) ---
cat(" Creating volcano plot...\n")

volcano_data <- all_inter_organ_correlations %>%
  filter(!is.na(Correlation) & !is.na(P_value)) %>%
  mutate(
    neg_log10_p = -log10(P_value),
    Label = ifelse(Significant & abs(Correlation) > 0.8,
                   paste(Gene1, Gene2, sep = "-"),
                   "")
  )

p_volcano <- ggplot(volcano_data, aes(x = Correlation,
                                       y = neg_log10_p,
                                       colour = Significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = rgb(237, 115, 116, maxColorValue = 255),
             linewidth = 0.75 / .pt) +
  geom_vline(xintercept = c(-0.7, 0.7), linetype = "dashed",
             colour = rgb(154, 153, 153, maxColorValue = 255),
             linewidth = 0.75 / .pt) +
  geom_text_repel(aes(label = Label),
                  size = 6 / .pt, 
                  max.overlaps = 15,
                  colour = "black") +
  scale_colour_manual(values = c("FALSE" = rgb(180, 180, 180, maxColorValue = 255),
                                 "TRUE" = rgb(237, 115, 116, maxColorValue = 255)),
                      name = "Significant\n(p<0.05)") +
  labs(
    title = "Metabolic-Immune Axis Correlation Volcano Plot",
    subtitle = "Gut-Pancreas-Liver-Kidney-Spleen axis (n=12)",
    x = "Pearson Correlation Coefficient (r)",
    y = "-log10(p-value)"
  ) +
  theme_diabetologia()

save_plot(p_volcano, "MetabolicImmune_Volcano_Plot", width = 14, height = 12)

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n=============================================================\n")
cat("METABOLIC-IMMUNE AXIS ANALYSIS COMPLETE\n")
cat("=============================================================\n\n")

cat("FILES GENERATED:\n\n")
cat("CSV Files:\n")
cat(" • MetabolicImmune_Axis_Correlation_Results.csv\n")
cat(" • MetabolicImmune_Axis_Significant_Correlations.csv\n\n")

cat("Visualisations (PNG + SVG, Diabetologia format):\n")
cat(" • Fig_5a_MetabolicImmune_Barplot_Summary (7×8 cm)\n")
cat(" • Fig_5b_MetabolicImmune_Lollipop_Top30 (10×14 cm)\n")
cat(" • Fig_5c_MetabolicImmune_Scatter_MultiPanel (14×15 cm, 2×3 panels)\n")
cat(" • Fig_5d_MetabolicImmune_Network_Strong_Correlations (10×9 cm)\n")
cat(" • MetabolicImmune_Volcano_Plot\n\n")

organs_present <- unique(c(all_inter_organ_correlations$Organ1,
                           all_inter_organ_correlations$Organ2))

cat("METABOLIC-IMMUNE AXIS FOCUS:\n")
cat(paste(" Organs:", paste(organs_present, collapse = ", "), "\n"))
cat(paste(" Total correlations:", nrow(all_inter_organ_correlations), "\n"))
cat(paste(" Significant (p<0.05):", sum(all_inter_organ_correlations$Significant), "\n"))
cat(paste(" Strong (|r|>0.7, p<0.01):", nrow(strong_cors), "\n\n"))

cat("=============================================================\n")
cat("All files saved to:", getwd(), "\n")
cat("=============================================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
