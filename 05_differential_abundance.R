# ============================================================================
# 05_differential_abundance.R
# Differential Abundance Analysis (DESeq2 and LEfSe)
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ============================================================================
# Configuration
# ============================================================================

PROCESSED_PATH <- "data/processed/"
RESULTS_PATH <- "results/"
FIG_PATH <- file.path(RESULTS_PATH, "figures")
TABLE_PATH <- file.path(RESULTS_PATH, "tables")

# Significance thresholds
FDR_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1

# Colors
GROUP_COLORS <- c(
  "Healthy" = "#4DAF4A",
  "Endometritis" = "#E41A1C", 
  "Anestrus" = "#377EB8"
)

set.seed(42)

# ============================================================================
# 1. Load data
# ============================================================================

message("\n========================================")
message("Step 1: Loading data")
message("========================================\n")

ps_health <- readRDS(file.path(PROCESSED_PATH, "phyloseq_health.rds"))
ps_endo <- readRDS(file.path(PROCESSED_PATH, "phyloseq_endometritis.rds"))

message("Data loaded successfully")

# ============================================================================
# 2. Aggregate to genus level
# ============================================================================

message("\n========================================")
message("Step 2: Aggregating to genus level")
message("========================================\n")

# Aggregate to genus
ps_genus <- tax_glom(ps_health, taxrank = "Genus", NArm = FALSE)
message(paste("Aggregated to", ntaxa(ps_genus), "genera"))

# ============================================================================
# 3. DESeq2 analysis - Health groups (Figure 1E, 1F)
# ============================================================================

message("\n========================================")
message("Step 3: DESeq2 differential abundance")
message("========================================\n")

# Convert to DESeq2 format
diagdds <- phyloseq_to_deseq2(ps_genus, ~ Group)

# Handle zero counts
geoMeans <- apply(counts(diagdds), 1, function(row) {
  if (all(row == 0)) { 0 } else { exp(mean(log(row[row != 0]))) }
})
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)

# Run DESeq2
diagdds <- DESeq(diagdds, test = "Wald", fitType = "local")

# Get results for each comparison
comparisons <- list(
  c("Endometritis", "Healthy"),
  c("Anestrus", "Healthy"),
  c("Anestrus", "Endometritis")
)

all_results <- data.frame()

for(comp in comparisons) {
  res <- results(diagdds, contrast = c("Group", comp[1], comp[2]))
  res_df <- as.data.frame(res)
  res_df$ASV <- rownames(res_df)
  
  # Add taxonomy
  tax <- as.data.frame(tax_table(ps_genus))
  res_df <- merge(res_df, tax, by.x = "ASV", by.y = "row.names")
  
  # Add comparison info
  res_df$Comparison <- paste(comp[1], "vs", comp[2])
  res_df$Direction <- ifelse(res_df$log2FoldChange > 0, 
                             paste("Higher in", comp[1]),
                             paste("Higher in", comp[2]))
  
  # Filter significant
  res_df$Significant <- res_df$padj < FDR_THRESHOLD & 
                        abs(res_df$log2FoldChange) > LFC_THRESHOLD
  
  all_results <- rbind(all_results, res_df)
  
  # Print summary
  sig_count <- sum(res_df$Significant, na.rm = TRUE)
  message(paste(comp[1], "vs", comp[2], ":", sig_count, "significant genera"))
}

# Save all results
write.csv(all_results, file.path(TABLE_PATH, "deseq2_results_health.csv"), row.names = FALSE)

# ============================================================================
# 4. Volcano plots
# ============================================================================

message("\n========================================")
message("Step 4: Generating volcano plots")
message("========================================\n")

# Function to create volcano plot
create_volcano <- function(res_df, comparison_name) {
  res_df <- res_df[!is.na(res_df$padj), ]
  
  res_df$Label <- ifelse(res_df$Significant, res_df$Genus, NA)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "blue") +
    ggrepel::geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15) +
    labs(
      title = comparison_name,
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

# Create volcano plots for each comparison
volcano_plots <- list()
for(comp in unique(all_results$Comparison)) {
  res_sub <- all_results[all_results$Comparison == comp, ]
  volcano_plots[[comp]] <- create_volcano(res_sub, comp)
}

# Combine and save
library(patchwork)
p_volcanos <- wrap_plots(volcano_plots, ncol = 3)
ggsave(file.path(FIG_PATH, "volcano_plots_health.pdf"), p_volcanos, width = 15, height = 5)
ggsave(file.path(FIG_PATH, "volcano_plots_health.png"), p_volcanos, width = 15, height = 5, dpi = 300)

message("Volcano plots saved")

# ============================================================================
# 5. Heatmap of significant genera (Figure 1F)
# ============================================================================

message("\n========================================")
message("Step 5: Generating heatmap")
message("========================================\n")

# Get significant genera
sig_genera <- unique(all_results$Genus[all_results$Significant == TRUE])
sig_genera <- sig_genera[!is.na(sig_genera)]

if(length(sig_genera) > 0) {
  # Get abundance data
  ps_ra <- transform_sample_counts(ps_genus, function(x) x/sum(x) * 100)
  otu <- as(otu_table(ps_ra), "matrix")
  if(!taxa_are_rows(ps_ra)) otu <- t(otu)
  
  # Get genus names
  tax <- as.data.frame(tax_table(ps_ra))
  rownames(otu) <- tax$Genus
  
  # Subset to significant genera
  otu_sig <- otu[rownames(otu) %in% sig_genera, ]
  
  # Log transform
  otu_sig <- log10(otu_sig + 0.01)
  
  # Annotation
  meta <- as(sample_data(ps_ra), "data.frame")
  anno_col <- data.frame(
    Group = meta$Group,
    row.names = colnames(otu_sig)
  )
  
  anno_colors <- list(Group = GROUP_COLORS)
  
  # Create heatmap
  pdf(file.path(FIG_PATH, "Fig1F_heatmap_significant_genera.pdf"), width = 12, height = 8)
  pheatmap(
    otu_sig,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    scale = "row",
    clustering_distance_cols = "correlation",
    clustering_distance_rows = "correlation",
    show_colnames = FALSE,
    fontsize_row = 8,
    main = "Differentially Abundant Genera"
  )
  dev.off()
  
  message("Figure 1F heatmap saved")
}

# ============================================================================
# 6. LEfSe-style analysis (Figure 1G)
# ============================================================================

message("\n========================================")
message("Step 6: LEfSe-style analysis")
message("========================================\n")

# Calculate LDA scores (simplified approach using log2FC * -log10(p))
lefse_results <- all_results %>%
  filter(Significant == TRUE, !is.na(Genus)) %>%
  mutate(
    LDA_score = abs(log2FoldChange) * -log10(padj),
    Direction = ifelse(log2FoldChange > 0, "Positive", "Negative")
  ) %>%
  group_by(Comparison) %>%
  arrange(desc(LDA_score)) %>%
  slice_head(n = 20)

# LEfSe bar plot
p_lefse <- ggplot(lefse_results, aes(x = reorder(Genus, LDA_score), y = LDA_score, fill = Comparison)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~Comparison, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "LEfSe Analysis: Biomarkers by Group",
    x = "",
    y = "LDA Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(FIG_PATH, "Fig1G_lefse_barplot.pdf"), p_lefse, width = 10, height = 12)
ggsave(file.path(FIG_PATH, "Fig1G_lefse_barplot.png"), p_lefse, width = 10, height = 12, dpi = 300)

# Save LEfSe results
write.csv(lefse_results, file.path(TABLE_PATH, "lefse_results.csv"), row.names = FALSE)

message("Figure 1G LEfSe plot saved")

# ============================================================================
# 7. Treatment comparison DESeq2
# ============================================================================

message("\n========================================")
message("Step 7: Treatment comparison DESeq2")
message("========================================\n")

# Aggregate to genus
ps_endo_genus <- tax_glom(ps_endo, taxrank = "Genus", NArm = FALSE)

# DESeq2 for treatment comparison
diagdds_treat <- phyloseq_to_deseq2(ps_endo_genus, ~ Treatment)

geoMeans_treat <- apply(counts(diagdds_treat), 1, function(row) {
  if (all(row == 0)) { 0 } else { exp(mean(log(row[row != 0]))) }
})
diagdds_treat <- estimateSizeFactors(diagdds_treat, geoMeans = geoMeans_treat)
diagdds_treat <- DESeq(diagdds_treat, test = "Wald", fitType = "local")

# Get results for LP vs Untreated
res_lp <- results(diagdds_treat, contrast = c("Treatment", "LP", "Untreated"))
res_lp_df <- as.data.frame(res_lp)
res_lp_df$ASV <- rownames(res_lp_df)

tax_treat <- as.data.frame(tax_table(ps_endo_genus))
res_lp_df <- merge(res_lp_df, tax_treat, by.x = "ASV", by.y = "row.names")
res_lp_df$Significant <- res_lp_df$padj < FDR_THRESHOLD & 
                          !is.na(res_lp_df$padj)

write.csv(res_lp_df, file.path(TABLE_PATH, "deseq2_LP_vs_Untreated.csv"), row.names = FALSE)

# Results for MS vs Untreated
res_ms <- results(diagdds_treat, contrast = c("Treatment", "MS", "Untreated"))
res_ms_df <- as.data.frame(res_ms)
res_ms_df$ASV <- rownames(res_ms_df)
res_ms_df <- merge(res_ms_df, tax_treat, by.x = "ASV", by.y = "row.names")
res_ms_df$Significant <- res_ms_df$padj < FDR_THRESHOLD & !is.na(res_ms_df$padj)

write.csv(res_ms_df, file.path(TABLE_PATH, "deseq2_MS_vs_Untreated.csv"), row.names = FALSE)

message("Treatment DESeq2 results saved")

# ============================================================================
# 8. Summary
# ============================================================================

message("\n========================================")
message("Differential Abundance Analysis Complete!")
message("========================================\n")

message("Output files:")
message("  - deseq2_results_health.csv")
message("  - deseq2_LP_vs_Untreated.csv")
message("  - deseq2_MS_vs_Untreated.csv")
message("  - lefse_results.csv")
message("  - volcano_plots_health.pdf/png")
message("  - Fig1F_heatmap_significant_genera.pdf")
message("  - Fig1G_lefse_barplot.pdf/png")
