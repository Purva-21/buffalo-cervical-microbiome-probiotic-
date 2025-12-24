# ============================================================================
# 06_simper_analysis.R
# SIMPER (Similarity Percentage) Analysis
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

# ============================================================================
# Configuration
# ============================================================================

PROCESSED_PATH <- "data/processed/"
RESULTS_PATH <- "results/"
FIG_PATH <- file.path(RESULTS_PATH, "figures")
TABLE_PATH <- file.path(RESULTS_PATH, "tables")

set.seed(42)

# ============================================================================
# 1. Load data
# ============================================================================

message("\n========================================")
message("Step 1: Loading data")
message("========================================\n")

ps_endo <- readRDS(file.path(PROCESSED_PATH, "phyloseq_endometritis.rds"))
message("Data loaded successfully")

# ============================================================================
# 2. Prepare data for SIMPER
# ============================================================================

message("\n========================================")
message("Step 2: Preparing data")
message("========================================\n")

# Aggregate to genus level
ps_genus <- tax_glom(ps_endo, taxrank = "Genus", NArm = FALSE)

# Get relative abundance
ps_ra <- transform_sample_counts(ps_genus, function(x) x/sum(x) * 100)

# Extract OTU table and metadata
otu <- as(otu_table(ps_ra), "matrix")
if(!taxa_are_rows(ps_ra)) otu <- t(otu)
otu <- t(otu)  # samples as rows

meta <- as(sample_data(ps_ra), "data.frame")

# Get genus names
tax <- as.data.frame(tax_table(ps_ra))
colnames(otu) <- tax$Genus

message(paste("Data prepared:", nrow(otu), "samples,", ncol(otu), "genera"))

# ============================================================================
# 3. SIMPER analysis - Treatment comparisons (Figure 2J)
# ============================================================================

message("\n========================================")
message("Step 3: SIMPER analysis")
message("========================================\n")

# Define treatment groups
treatments <- c("LP", "MS", "NS")
all_simper <- list()

for(treat in treatments) {
  # Compare treatment to Untreated
  samples_use <- rownames(meta[meta$Treatment %in% c(treat, "Untreated"), ])
  otu_sub <- otu[samples_use, ]
  group_sub <- meta[samples_use, "Treatment"]
  
  # Run SIMPER
  simper_res <- simper(otu_sub, group_sub, permutations = 999)
  
  # Extract results
  comparison <- paste(treat, "vs Untreated")
  simper_summary <- summary(simper_res)[[1]]
  
  simper_df <- data.frame(
    Genus = rownames(simper_summary),
    Contribution = simper_summary$average * 100,
    Cumulative = simper_summary$cumsum * 100,
    Mean_Treat = simper_summary$ava,
    Mean_Untreated = simper_summary$avb,
    P_value = simper_summary$p,
    Comparison = comparison
  )
  
  # Top contributors
  simper_df <- simper_df %>%
    arrange(desc(Contribution)) %>%
    head(20)
  
  all_simper[[comparison]] <- simper_df
  
  # Print summary
  message(paste("\n", comparison, "- Overall dissimilarity:", 
                round(sum(simper_summary$average) * 100, 1), "%"))
  message("Top 5 contributors:")
  print(head(simper_df[, c("Genus", "Contribution", "P_value")], 5))
}

# Combine all results
simper_combined <- do.call(rbind, all_simper)
write.csv(simper_combined, file.path(TABLE_PATH, "simper_results.csv"), row.names = FALSE)

# ============================================================================
# 4. SIMPER visualization
# ============================================================================

message("\n========================================")
message("Step 4: Generating SIMPER plots")
message("========================================\n")

# Bar plot of top contributors
p_simper <- ggplot(simper_combined, aes(x = reorder(Genus, Contribution), 
                                         y = Contribution, 
                                         fill = Comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_wrap(~Comparison, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "SIMPER Analysis: Key Contributors to Community Dissimilarity",
    x = "",
    y = "Contribution to Dissimilarity (%)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(FIG_PATH, "Fig2J_simper_barplot.pdf"), p_simper, width = 10, height = 12)
ggsave(file.path(FIG_PATH, "Fig2J_simper_barplot.png"), p_simper, width = 10, height = 12, dpi = 300)

# Heatmap of SIMPER results
simper_wide <- simper_combined %>%
  select(Genus, Contribution, Comparison) %>%
  pivot_wider(names_from = Comparison, values_from = Contribution, values_fill = 0)

simper_matrix <- as.matrix(simper_wide[, -1])
rownames(simper_matrix) <- simper_wide$Genus

# Order by total contribution
row_order <- order(rowSums(simper_matrix), decreasing = TRUE)
simper_matrix <- simper_matrix[row_order[1:min(25, nrow(simper_matrix))], ]

pdf(file.path(FIG_PATH, "simper_heatmap.pdf"), width = 8, height = 10)
pheatmap::pheatmap(
  simper_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "SIMPER Contribution Heatmap",
  color = colorRampPalette(c("white", "orange", "red"))(50),
  fontsize_row = 9
)
dev.off()

message("SIMPER plots saved")

# ============================================================================
# 5. Core contributors analysis
# ============================================================================

message("\n========================================")
message("Step 5: Core contributors analysis")
message("========================================\n")

# Identify genera that consistently contribute across comparisons
core_contributors <- simper_combined %>%
  group_by(Genus) %>%
  summarise(
    Mean_Contribution = mean(Contribution),
    N_comparisons = n(),
    Min_Contribution = min(Contribution),
    Max_Contribution = max(Contribution)
  ) %>%
  filter(N_comparisons == 3) %>%  # Present in all comparisons
  arrange(desc(Mean_Contribution))

message("\nCore contributors (present in all treatment comparisons):")
print(head(core_contributors, 10))

write.csv(core_contributors, file.path(TABLE_PATH, "simper_core_contributors.csv"), row.names = FALSE)

# ============================================================================
# 6. Dissimilarity summary
# ============================================================================

message("\n========================================")
message("Step 6: Dissimilarity summary")
message("========================================\n")

# Calculate overall dissimilarity for each comparison
dissim_summary <- simper_combined %>%
  group_by(Comparison) %>%
  summarise(
    Total_Dissimilarity = sum(Contribution),
    Top10_Contribution = sum(Contribution[1:10]),
    N_significant = sum(P_value < 0.05, na.rm = TRUE)
  )

message("\nDissimilarity summary:")
print(dissim_summary)

write.csv(dissim_summary, file.path(TABLE_PATH, "simper_dissimilarity_summary.csv"), row.names = FALSE)

message("\n========================================")
message("SIMPER Analysis Complete!")
message("========================================\n")
