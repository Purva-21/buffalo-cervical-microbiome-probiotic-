# ============================================================================
# 03_alpha_diversity.R
# Alpha Diversity Analysis
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# ============================================================================
# Configuration
# ============================================================================

PROCESSED_PATH <- "data/processed/"
RESULTS_PATH <- "results/"
FIG_PATH <- file.path(RESULTS_PATH, "figures")
TABLE_PATH <- file.path(RESULTS_PATH, "tables")

# Color palette
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

# ============================================================================
# 1. Load data
# ============================================================================

message("\n========================================")
message("Step 1: Loading data")
message("========================================\n")

ps_rare <- readRDS(file.path(PROCESSED_PATH, "phyloseq_rarefied.rds"))
ps_endo <- readRDS(file.path(PROCESSED_PATH, "phyloseq_endometritis.rds"))
ps_anestrus <- readRDS(file.path(PROCESSED_PATH, "phyloseq_anestrus.rds"))

message("Data loaded successfully")

# ============================================================================
# 2. Calculate alpha diversity indices
# ============================================================================

message("\n========================================")
message("Step 2: Calculating alpha diversity")
message("========================================\n")

# Function to calculate all alpha diversity metrics
calc_alpha_div <- function(ps_obj) {
  # Get OTU table
  otu <- as(otu_table(ps_obj), "matrix")
  if(!taxa_are_rows(ps_obj)) otu <- t(otu)
  otu <- t(otu)  # samples as rows
  
  # Calculate indices
  alpha_df <- data.frame(
    SampleID = sample_names(ps_obj),
    Observed = specnumber(otu),
    Chao1 = estimateR(otu)["S.chao1", ],
    Shannon = diversity(otu, index = "shannon"),
    Simpson = diversity(otu, index = "simpson"),
    InvSimpson = diversity(otu, index = "invsimpson"),
    Pielou = diversity(otu, index = "shannon") / log(specnumber(otu))
  )
  
  # Add metadata
  meta <- as(sample_data(ps_obj), "data.frame")
  alpha_df <- merge(alpha_df, meta, by.x = "SampleID", by.y = "row.names")
  
  return(alpha_df)
}

# Calculate for health comparison
alpha_health <- calc_alpha_div(ps_rare)
message(paste("Calculated alpha diversity for", nrow(alpha_health), "samples"))

# ============================================================================
# 3. Statistical testing - Health groups
# ============================================================================

message("\n========================================")
message("Step 3: Statistical testing (Health groups)")
message("========================================\n")

# Kruskal-Wallis tests
kw_results <- data.frame()

for(metric in c("Observed", "Chao1", "Shannon", "Simpson")) {
  # Kruskal-Wallis test
  kw_test <- kruskal.test(as.formula(paste(metric, "~ Group")), data = alpha_health)
  
  # Pairwise Wilcoxon tests
  pw_test <- pairwise.wilcox.test(alpha_health[[metric]], alpha_health$Group, 
                                   p.adjust.method = "BH")
  
  kw_results <- rbind(kw_results, data.frame(
    Metric = metric,
    Test = "Kruskal-Wallis",
    Statistic = kw_test$statistic,
    P_value = kw_test$p.value,
    Significance = ifelse(kw_test$p.value < 0.001, "***",
                          ifelse(kw_test$p.value < 0.01, "**",
                                 ifelse(kw_test$p.value < 0.05, "*", "ns")))
  ))
  
  message(paste(metric, "- Kruskal-Wallis p-value:", round(kw_test$p.value, 4)))
}

# Save results
write.csv(kw_results, file.path(TABLE_PATH, "alpha_diversity_stats_health.csv"), row.names = FALSE)

# ============================================================================
# 4. Plot alpha diversity - Health groups (Figure 1B)
# ============================================================================

message("\n========================================")
message("Step 4: Generating alpha diversity plots")
message("========================================\n")

# Chao1 plot
p_chao1 <- ggplot(alpha_health, aes(x = Group, y = Chao1, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = GROUP_COLORS) +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_health$Chao1) * 1.1) +
  labs(title = "Species Richness (Chao1)", x = "", y = "Chao1 Index") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Shannon plot
p_shannon <- ggplot(alpha_health, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = GROUP_COLORS) +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_health$Shannon) * 1.1) +
  stat_compare_means(comparisons = list(c("Healthy", "Endometritis"),
                                         c("Healthy", "Anestrus"),
                                         c("Endometritis", "Anestrus")),
                     method = "wilcox.test", label = "p.signif") +
  labs(title = "Shannon Diversity", x = "", y = "Shannon Index") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Combined plot
library(patchwork)
p_combined_health <- p_chao1 + p_shannon + 
  plot_annotation(title = "Alpha Diversity: Health Status Comparison",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

# Save
ggsave(file.path(FIG_PATH, "Fig1B_alpha_diversity_health.pdf"), 
       p_combined_health, width = 10, height = 5)
ggsave(file.path(FIG_PATH, "Fig1B_alpha_diversity_health.png"), 
       p_combined_health, width = 10, height = 5, dpi = 300)

message("Figure 1B saved")

# ============================================================================
# 5. Alpha diversity - Treatment comparison (Figure 2B)
# ============================================================================

message("\n========================================")
message("Step 5: Treatment comparison alpha diversity")
message("========================================\n")

# Calculate for endometritis treatment groups
ps_endo_rare <- rarefy_even_depth(ps_endo, sample.size = min(sample_sums(ps_endo)), 
                                   verbose = FALSE, rngseed = 42)
alpha_endo <- calc_alpha_div(ps_endo_rare)

# Statistical tests
kw_endo <- data.frame()
for(metric in c("Observed", "Chao1", "Shannon")) {
  kw_test <- kruskal.test(as.formula(paste(metric, "~ Treatment")), data = alpha_endo)
  kw_endo <- rbind(kw_endo, data.frame(
    Metric = metric,
    Chi_sq = kw_test$statistic,
    P_value = kw_test$p.value
  ))
}

write.csv(kw_endo, file.path(TABLE_PATH, "alpha_diversity_stats_treatment.csv"), row.names = FALSE)

# Plot
p_chao1_treat <- ggplot(alpha_endo, aes(x = Treatment, y = Chao1, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  facet_wrap(~Timepoint, nrow = 1) +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_endo$Chao1, na.rm = TRUE) * 1.1) +
  labs(title = "Species Richness by Treatment", x = "", y = "Chao1 Index") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_shannon_treat <- ggplot(alpha_endo, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  facet_wrap(~Timepoint, nrow = 1) +
  stat_compare_means(method = "kruskal.test", label.y = max(alpha_endo$Shannon, na.rm = TRUE) * 1.1) +
  labs(title = "Shannon Diversity by Treatment", x = "", y = "Shannon Index") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_combined_treat <- p_chao1_treat / p_shannon_treat

ggsave(file.path(FIG_PATH, "Fig2B_alpha_diversity_treatment.pdf"), 
       p_combined_treat, width = 12, height = 10)
ggsave(file.path(FIG_PATH, "Fig2B_alpha_diversity_treatment.png"), 
       p_combined_treat, width = 12, height = 10, dpi = 300)

message("Figure 2B saved")

# ============================================================================
# 6. Alpha diversity - Pregnancy outcome (Figure 3E)
# ============================================================================

message("\n========================================")
message("Step 6: Pregnancy outcome comparison")
message("========================================\n")

# Add pregnancy status to alpha diversity data
alpha_pregnancy <- alpha_endo[!is.na(alpha_endo$Pregnancy_Status), ]

if(nrow(alpha_pregnancy) > 0) {
  # Statistical test
  wilcox_preg <- wilcox.test(Chao1 ~ Pregnancy_Status, data = alpha_pregnancy)
  message(paste("Chao1 Pregnant vs Non-pregnant: p =", round(wilcox_preg$p.value, 4)))
  
  # Plot
  p_preg <- ggplot(alpha_pregnancy, aes(x = Pregnancy_Status, y = Chao1, fill = Pregnancy_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    scale_fill_manual(values = c("Pregnant" = "#4DAF4A", "Non-pregnant" = "#E41A1C")) +
    stat_compare_means(method = "wilcox.test") +
    labs(title = "Alpha Diversity by Pregnancy Outcome",
         x = "", y = "Chao1 Index") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(file.path(FIG_PATH, "Fig3E_alpha_diversity_pregnancy.pdf"), 
         p_preg, width = 6, height = 5)
  ggsave(file.path(FIG_PATH, "Fig3E_alpha_diversity_pregnancy.png"), 
         p_preg, width = 6, height = 5, dpi = 300)
  
  message("Figure 3E saved")
}

# ============================================================================
# 7. Save alpha diversity data
# ============================================================================

message("\n========================================")
message("Step 7: Saving alpha diversity data")
message("========================================\n")

write.csv(alpha_health, file.path(TABLE_PATH, "alpha_diversity_health.csv"), row.names = FALSE)
write.csv(alpha_endo, file.path(TABLE_PATH, "alpha_diversity_treatment.csv"), row.names = FALSE)

message("Alpha diversity data saved")

message("\n========================================")
message("Alpha Diversity Analysis Complete!")
message("========================================\n")
