# ============================================================================
# 04_beta_diversity.R
# Beta Diversity and Ordination Analysis
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ape)

# ============================================================================
# Configuration
# ============================================================================

PROCESSED_PATH <- "data/processed/"
RESULTS_PATH <- "results/"
FIG_PATH <- file.path(RESULTS_PATH, "figures")
TABLE_PATH <- file.path(RESULTS_PATH, "tables")

# Color palettes
GROUP_COLORS <- c(
  "Healthy" = "#4DAF4A",
  "Endometritis" = "#E41A1C", 
  "Anestrus" = "#377EB8"
)

TREATMENT_COLORS <- c(
  "Untreated" = "#999999",
  "LP" = "#E41A1C",
  "MS" = "#377EB8",
  "NS" = "#4DAF4A"
)

set.seed(42)

# ============================================================================
# 1. Load data
# ============================================================================

message("\n========================================")
message("Step 1: Loading data")
message("========================================\n")

ps_health <- readRDS(file.path(PROCESSED_PATH, "phyloseq_health.rds"))
ps_relabund <- readRDS(file.path(PROCESSED_PATH, "phyloseq_relabund.rds"))
ps_endo <- readRDS(file.path(PROCESSED_PATH, "phyloseq_endometritis.rds"))

message("Data loaded successfully")

# ============================================================================
# 2. Calculate Bray-Curtis distances
# ============================================================================

message("\n========================================")
message("Step 2: Calculating distance matrices")
message("========================================\n")

# Transform to relative abundance for beta diversity
ps_ra <- transform_sample_counts(ps_health, function(x) x/sum(x))

# Calculate Bray-Curtis dissimilarity
bray_dist <- phyloseq::distance(ps_ra, method = "bray")
message("Bray-Curtis distance matrix calculated")

# ============================================================================
# 3. PERMANOVA analysis - Health groups
# ============================================================================

message("\n========================================")
message("Step 3: PERMANOVA analysis (Health groups)")
message("========================================\n")

# Get metadata
meta <- as(sample_data(ps_ra), "data.frame")

# PERMANOVA test
permanova_health <- adonis2(
  bray_dist ~ Group, 
  data = meta, 
  permutations = 999,
  by = "margin"
)

message("PERMANOVA Results (Health Status):")
print(permanova_health)

# Save results
permanova_df <- as.data.frame(permanova_health)
permanova_df$Comparison <- "Health_Status"
write.csv(permanova_df, file.path(TABLE_PATH, "permanova_health.csv"))

# Effect size
r2_health <- permanova_health$R2[1]
p_health <- permanova_health$`Pr(>F)`[1]
message(paste("R² =", round(r2_health, 4), ", p =", p_health))

# ============================================================================
# 4. PCoA ordination - Health groups (Figure 1C)
# ============================================================================

message("\n========================================")
message("Step 4: PCoA ordination (Health groups)")
message("========================================\n")

# Perform PCoA
pcoa <- ordinate(ps_ra, method = "PCoA", distance = "bray")

# Extract variance explained
var_explained <- round(pcoa$values$Relative_eig[1:2] * 100, 1)

# Create plot
p_pcoa_health <- plot_ordination(ps_ra, pcoa, color = "Group") +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95, linetype = 2) +
  scale_color_manual(values = GROUP_COLORS) +
  labs(
    title = "PCoA: Bray-Curtis Dissimilarity",
    subtitle = paste("PERMANOVA: R² =", round(r2_health, 3), ", p =", p_health),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(FIG_PATH, "Fig1C_beta_diversity_health.pdf"), 
       p_pcoa_health, width = 8, height = 6)
ggsave(file.path(FIG_PATH, "Fig1C_beta_diversity_health.png"), 
       p_pcoa_health, width = 8, height = 6, dpi = 300)

message("Figure 1C saved")

# ============================================================================
# 5. Pairwise PERMANOVA - Health groups
# ============================================================================

message("\n========================================")
message("Step 5: Pairwise PERMANOVA")
message("========================================\n")

# Pairwise comparisons
groups <- unique(meta$Group)
pairwise_results <- data.frame()

for(i in 1:(length(groups)-1)) {
  for(j in (i+1):length(groups)) {
    # Subset data
    subset_samples <- rownames(meta[meta$Group %in% c(groups[i], groups[j]), ])
    ps_sub <- prune_samples(subset_samples, ps_ra)
    dist_sub <- phyloseq::distance(ps_sub, method = "bray")
    meta_sub <- as(sample_data(ps_sub), "data.frame")
    
    # PERMANOVA
    perm_res <- adonis2(dist_sub ~ Group, data = meta_sub, permutations = 999)
    
    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1 = groups[i],
      Group2 = groups[j],
      R2 = perm_res$R2[1],
      F_value = perm_res$F[1],
      P_value = perm_res$`Pr(>F)`[1]
    ))
    
    message(paste(groups[i], "vs", groups[j], ": R² =", 
                  round(perm_res$R2[1], 3), ", p =", perm_res$`Pr(>F)`[1]))
  }
}

# FDR correction
pairwise_results$P_adj <- p.adjust(pairwise_results$P_value, method = "BH")
write.csv(pairwise_results, file.path(TABLE_PATH, "permanova_pairwise_health.csv"), row.names = FALSE)

# ============================================================================
# 6. Beta diversity - Treatment comparison (Figure 2E, 2F)
# ============================================================================

message("\n========================================")
message("Step 6: Treatment comparison beta diversity")
message("========================================\n")

# Transform endometritis data
ps_endo_ra <- transform_sample_counts(ps_endo, function(x) x/sum(x))
bray_endo <- phyloseq::distance(ps_endo_ra, method = "bray")
meta_endo <- as(sample_data(ps_endo_ra), "data.frame")

# PERMANOVA for treatment
permanova_treat <- adonis2(
  bray_endo ~ Treatment, 
  data = meta_endo, 
  permutations = 999
)

message("\nPERMANOVA Results (Treatment):")
print(permanova_treat)

# Pairwise treatment comparisons
treatments <- unique(meta_endo$Treatment)
treat_pairwise <- data.frame()

for(i in 1:(length(treatments)-1)) {
  for(j in (i+1):length(treatments)) {
    subset_samples <- rownames(meta_endo[meta_endo$Treatment %in% c(treatments[i], treatments[j]), ])
    if(length(subset_samples) >= 4) {
      ps_sub <- prune_samples(subset_samples, ps_endo_ra)
      dist_sub <- phyloseq::distance(ps_sub, method = "bray")
      meta_sub <- as(sample_data(ps_sub), "data.frame")
      
      perm_res <- adonis2(dist_sub ~ Treatment, data = meta_sub, permutations = 999)
      
      treat_pairwise <- rbind(treat_pairwise, data.frame(
        Treatment1 = treatments[i],
        Treatment2 = treatments[j],
        R2 = perm_res$R2[1],
        F_value = perm_res$F[1],
        P_value = perm_res$`Pr(>F)`[1]
      ))
    }
  }
}

treat_pairwise$P_adj <- p.adjust(treat_pairwise$P_value, method = "BH")
write.csv(treat_pairwise, file.path(TABLE_PATH, "permanova_pairwise_treatment.csv"), row.names = FALSE)

# PCoA plot
pcoa_endo <- ordinate(ps_endo_ra, method = "PCoA", distance = "bray")
var_exp_endo <- round(pcoa_endo$values$Relative_eig[1:2] * 100, 1)

p_pcoa_treat <- plot_ordination(ps_endo_ra, pcoa_endo, color = "Treatment", shape = "Timepoint") +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Treatment), type = "norm", level = 0.95, linetype = 2) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(
    title = "PCoA: Treatment Comparison",
    subtitle = paste("PERMANOVA: R² =", round(permanova_treat$R2[1], 3), 
                     ", p =", permanova_treat$`Pr(>F)`[1]),
    x = paste0("PCoA1 (", var_exp_endo[1], "%)"),
    y = paste0("PCoA2 (", var_exp_endo[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(file.path(FIG_PATH, "Fig2E_beta_diversity_treatment.pdf"), 
       p_pcoa_treat, width = 9, height = 6)
ggsave(file.path(FIG_PATH, "Fig2E_beta_diversity_treatment.png"), 
       p_pcoa_treat, width = 9, height = 6, dpi = 300)

message("Figure 2E saved")

# ============================================================================
# 7. Dispersion analysis
# ============================================================================

message("\n========================================")
message("Step 7: Beta dispersion analysis")
message("========================================\n")

# Test homogeneity of dispersions
bdisp_health <- betadisper(bray_dist, meta$Group)
permutest_health <- permutest(bdisp_health)

message("Beta dispersion test (Health groups):")
print(permutest_health)

# Plot dispersion
disp_df <- data.frame(
  Group = meta$Group,
  Distance = bdisp_health$distances
)

p_dispersion <- ggplot(disp_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(
    title = "Distance to Group Centroid",
    x = "", y = "Distance to centroid"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(FIG_PATH, "beta_dispersion_health.pdf"), 
       p_dispersion, width = 6, height = 5)

# ============================================================================
# 8. Save outputs
# ============================================================================

message("\n========================================")
message("Step 8: Saving results")
message("========================================\n")

# Save distance matrix
saveRDS(bray_dist, file.path(PROCESSED_PATH, "bray_curtis_health.rds"))
saveRDS(bray_endo, file.path(PROCESSED_PATH, "bray_curtis_treatment.rds"))

# Save PCoA results
saveRDS(pcoa, file.path(PROCESSED_PATH, "pcoa_health.rds"))
saveRDS(pcoa_endo, file.path(PROCESSED_PATH, "pcoa_treatment.rds"))

message("Results saved successfully")

message("\n========================================")
message("Beta Diversity Analysis Complete!")
message("========================================\n")
