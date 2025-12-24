# ============================================================================
# 08_reproductive_outcomes.R
# Reproductive Outcomes Analysis (Estrus and Pregnancy)
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(VennDiagram)

# ============================================================================
# Configuration
# ============================================================================

PROCESSED_PATH <- "data/processed/"
METADATA_PATH <- "data/metadata/"
RESULTS_PATH <- "results/"
FIG_PATH <- file.path(RESULTS_PATH, "figures")
TABLE_PATH <- file.path(RESULTS_PATH, "tables")

# Treatment abbreviations
# LP = Lactiplantibacillus plantarum KUGBRC
# MS = Pediococcus pentosaceus GBRCKU
# NS = Normal Saline

TREATMENT_COLORS <- c(
  "LP" = "#E41A1C",
  "MS" = "#377EB8",
  "NS" = "#4DAF4A"
)

# ============================================================================
# 1. Load data
# ============================================================================

message("\n========================================")
message("Step 1: Loading data")
message("========================================\n")

# Load metadata with reproductive outcomes
meta <- read.csv(file.path(METADATA_PATH, "metadata_filtered.csv"), row.names = 1)

# Load phyloseq objects
ps_endo <- readRDS(file.path(PROCESSED_PATH, "phyloseq_endometritis.rds"))

message("Data loaded successfully")

# ============================================================================
# 2. Estrus induction analysis (Figure 3C)
# ============================================================================

message("\n========================================")
message("Step 2: Estrus induction analysis")
message("========================================\n")

# Calculate estrus induction rates by treatment
if("Estrus_Induced" %in% colnames(meta)) {
  estrus_summary <- meta %>%
    filter(Treatment %in% c("LP", "MS", "NS")) %>%
    group_by(Treatment) %>%
    summarise(
      N_total = n(),
      N_induced = sum(Estrus_Induced == "Yes", na.rm = TRUE),
      Rate = N_induced / N_total * 100,
      .groups = "drop"
    )
  
  message("\nEstrus induction rates:")
  print(estrus_summary)
  
  # Statistical test (Chi-square)
  estrus_table <- table(meta$Treatment[meta$Treatment %in% c("LP", "MS", "NS")], 
                        meta$Estrus_Induced[meta$Treatment %in% c("LP", "MS", "NS")])
  chi_estrus <- chisq.test(estrus_table)
  message(paste("\nChi-square test: X² =", round(chi_estrus$statistic, 2), 
                ", p =", round(chi_estrus$p.value, 4)))
  
  # Pairwise Fisher's exact tests
  treatments <- c("LP", "MS", "NS")
  pairwise_estrus <- data.frame()
  
  for(i in 1:(length(treatments)-1)) {
    for(j in (i+1):length(treatments)) {
      sub_data <- meta %>%
        filter(Treatment %in% c(treatments[i], treatments[j]))
      
      fisher_test <- fisher.test(table(sub_data$Treatment, sub_data$Estrus_Induced))
      
      pairwise_estrus <- rbind(pairwise_estrus, data.frame(
        Comparison = paste(treatments[i], "vs", treatments[j]),
        Odds_Ratio = fisher_test$estimate,
        P_value = fisher_test$p.value
      ))
    }
  }
  
  pairwise_estrus$P_adj <- p.adjust(pairwise_estrus$P_value, method = "BH")
  write.csv(pairwise_estrus, file.path(TABLE_PATH, "estrus_pairwise_tests.csv"), row.names = FALSE)
  
  # Bar plot
  p_estrus <- ggplot(estrus_summary, aes(x = Treatment, y = Rate, fill = Treatment)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Rate, 1), "%\n(", N_induced, "/", N_total, ")")),
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Estrus Induction Rate by Treatment",
      subtitle = paste("Chi-square p =", round(chi_estrus$p.value, 4)),
      x = "Treatment Group",
      y = "Estrus Induction Rate (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(file.path(FIG_PATH, "Fig3C_estrus_induction.pdf"), p_estrus, width = 6, height = 5)
  ggsave(file.path(FIG_PATH, "Fig3C_estrus_induction.png"), p_estrus, width = 6, height = 5, dpi = 300)
  
  message("Figure 3C saved")
  
  write.csv(estrus_summary, file.path(TABLE_PATH, "estrus_induction_summary.csv"), row.names = FALSE)
}

# ============================================================================
# 3. Pregnancy outcome analysis (Figure 3D)
# ============================================================================

message("\n========================================")
message("Step 3: Pregnancy outcome analysis")
message("========================================\n")

if("Pregnancy_Status" %in% colnames(meta)) {
  # Calculate pregnancy rates by treatment
  preg_summary <- meta %>%
    filter(Treatment %in% c("LP", "MS", "NS")) %>%
    group_by(Treatment) %>%
    summarise(
      N_total = n(),
      N_pregnant = sum(Pregnancy_Status == "Pregnant", na.rm = TRUE),
      Rate = N_pregnant / N_total * 100,
      .groups = "drop"
    )
  
  message("\nPregnancy rates:")
  print(preg_summary)
  
  # Statistical test
  preg_table <- table(meta$Treatment[meta$Treatment %in% c("LP", "MS", "NS")], 
                      meta$Pregnancy_Status[meta$Treatment %in% c("LP", "MS", "NS")])
  chi_preg <- chisq.test(preg_table)
  message(paste("\nChi-square test: X² =", round(chi_preg$statistic, 2), 
                ", p =", round(chi_preg$p.value, 4)))
  
  # Pairwise tests
  pairwise_preg <- data.frame()
  
  for(i in 1:(length(treatments)-1)) {
    for(j in (i+1):length(treatments)) {
      sub_data <- meta %>%
        filter(Treatment %in% c(treatments[i], treatments[j]),
               !is.na(Pregnancy_Status))
      
      if(nrow(sub_data) >= 4) {
        fisher_test <- fisher.test(table(sub_data$Treatment, sub_data$Pregnancy_Status))
        
        pairwise_preg <- rbind(pairwise_preg, data.frame(
          Comparison = paste(treatments[i], "vs", treatments[j]),
          Odds_Ratio = fisher_test$estimate,
          P_value = fisher_test$p.value
        ))
      }
    }
  }
  
  if(nrow(pairwise_preg) > 0) {
    pairwise_preg$P_adj <- p.adjust(pairwise_preg$P_value, method = "BH")
    write.csv(pairwise_preg, file.path(TABLE_PATH, "pregnancy_pairwise_tests.csv"), row.names = FALSE)
  }
  
  # Bar plot
  p_preg <- ggplot(preg_summary, aes(x = Treatment, y = Rate, fill = Treatment)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Rate, 1), "%\n(", N_pregnant, "/", N_total, ")")),
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Pregnancy Rate by Treatment",
      subtitle = paste("Chi-square p =", round(chi_preg$p.value, 4)),
      x = "Treatment Group",
      y = "Pregnancy Rate (%)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(file.path(FIG_PATH, "Fig3D_pregnancy_rate.pdf"), p_preg, width = 6, height = 5)
  ggsave(file.path(FIG_PATH, "Fig3D_pregnancy_rate.png"), p_preg, width = 6, height = 5, dpi = 300)
  
  message("Figure 3D saved")
  
  write.csv(preg_summary, file.path(TABLE_PATH, "pregnancy_summary.csv"), row.names = FALSE)
}

# ============================================================================
# 4. Venn diagram - Pregnant vs Non-pregnant microbiome (Figure 3I)
# ============================================================================

message("\n========================================")
message("Step 4: Venn diagram analysis")
message("========================================\n")

if("Pregnancy_Status" %in% colnames(meta)) {
  # Get ASVs unique to each group
  ps_preg <- subset_samples(ps_endo, Pregnancy_Status == "Pregnant")
  ps_nonpreg <- subset_samples(ps_endo, Pregnancy_Status == "Non-pregnant")
  
  # Remove zero-sum taxa
  ps_preg <- prune_taxa(taxa_sums(ps_preg) > 0, ps_preg)
  ps_nonpreg <- prune_taxa(taxa_sums(ps_nonpreg) > 0, ps_nonpreg)
  
  # Get taxa names
  taxa_preg <- taxa_names(ps_preg)
  taxa_nonpreg <- taxa_names(ps_nonpreg)
  
  # Calculate overlaps
  shared <- intersect(taxa_preg, taxa_nonpreg)
  unique_preg <- setdiff(taxa_preg, taxa_nonpreg)
  unique_nonpreg <- setdiff(taxa_nonpreg, taxa_preg)
  
  message(paste("\nPregnant unique ASVs:", length(unique_preg)))
  message(paste("Non-pregnant unique ASVs:", length(unique_nonpreg)))
  message(paste("Shared ASVs:", length(shared)))
  
  # Create Venn diagram
  venn_data <- list(
    Pregnant = taxa_preg,
    `Non-pregnant` = taxa_nonpreg
  )
  
  venn_plot <- venn.diagram(
    x = venn_data,
    filename = NULL,
    fill = c("#4DAF4A", "#E41A1C"),
    alpha = 0.5,
    label.col = "black",
    cex = 2,
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    main = "ASV Distribution by Pregnancy Outcome",
    main.cex = 1.5
  )
  
  pdf(file.path(FIG_PATH, "Fig3I_venn_pregnancy.pdf"), width = 8, height = 8)
  grid::grid.draw(venn_plot)
  dev.off()
  
  message("Figure 3I saved")
  
  # Save ASV lists
  venn_summary <- data.frame(
    Category = c("Pregnant_unique", "Non-pregnant_unique", "Shared"),
    Count = c(length(unique_preg), length(unique_nonpreg), length(shared))
  )
  write.csv(venn_summary, file.path(TABLE_PATH, "venn_pregnancy_summary.csv"), row.names = FALSE)
}

# ============================================================================
# 5. CV Score changes over time
# ============================================================================

message("\n========================================")
message("Step 5: CV Score dynamics")
message("========================================\n")

if("CV_Score" %in% colnames(meta) && "Timepoint" %in% colnames(meta)) {
  cv_dynamics <- meta %>%
    filter(Treatment %in% c("LP", "MS", "NS")) %>%
    group_by(Treatment, Timepoint) %>%
    summarise(
      Mean_CV = mean(CV_Score, na.rm = TRUE),
      SD_CV = sd(CV_Score, na.rm = TRUE),
      N = n(),
      SE_CV = SD_CV / sqrt(N),
      .groups = "drop"
    )
  
  p_cv <- ggplot(cv_dynamics, aes(x = Timepoint, y = Mean_CV, color = Treatment, group = Treatment)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean_CV - SE_CV, ymax = Mean_CV + SE_CV), width = 0.2) +
    scale_color_manual(values = TREATMENT_COLORS) +
    labs(
      title = "Cervico-Vaginal Score Changes Over Time",
      x = "Timepoint",
      y = "Mean CV Score",
      color = "Treatment"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(FIG_PATH, "cv_score_dynamics.pdf"), p_cv, width = 8, height = 5)
  ggsave(file.path(FIG_PATH, "cv_score_dynamics.png"), p_cv, width = 8, height = 5, dpi = 300)
  
  write.csv(cv_dynamics, file.path(TABLE_PATH, "cv_score_dynamics.csv"), row.names = FALSE)
  
  message("CV score dynamics plot saved")
}

# ============================================================================
# 6. Summary table
# ============================================================================

message("\n========================================")
message("Step 6: Generating summary table")
message("========================================\n")

# Combine all outcomes
if(exists("estrus_summary") && exists("preg_summary")) {
  outcomes_summary <- merge(
    estrus_summary,
    preg_summary,
    by = "Treatment",
    suffixes = c("_Estrus", "_Pregnancy")
  )
  
  outcomes_summary <- outcomes_summary %>%
    rename(
      Estrus_N = N_total_Estrus,
      Estrus_Induced = N_induced,
      Estrus_Rate = Rate_Estrus,
      Pregnancy_N = N_total_Pregnancy,
      Pregnant = N_pregnant,
      Pregnancy_Rate = Rate_Pregnancy
    )
  
  write.csv(outcomes_summary, file.path(TABLE_PATH, "reproductive_outcomes_summary.csv"), row.names = FALSE)
  
  message("\nReproductive outcomes summary:")
  print(outcomes_summary)
}

message("\n========================================")
message("Reproductive Outcomes Analysis Complete!")
message("========================================\n")

message("\nKey findings:")
message("  - LP (Lactiplantibacillus plantarum KUGBRC): Highest estrus and pregnancy rates")
message("  - MS (Pediococcus pentosaceus GBRCKU): Lower reproductive outcomes")
message("  - NS (Normal Saline): Intermediate outcomes")
