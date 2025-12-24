# ============================================================================
# 07_correlation_analysis.R
# Pathogen-Severity Correlation Analysis
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(corrplot)

# ============================================================================
# Configuration
# ============================================================================

PROCESSED_PATH <- "data/processed/"
RESULTS_PATH <- "results/"
FIG_PATH <- file.path(RESULTS_PATH, "figures")
TABLE_PATH <- file.path(RESULTS_PATH, "tables")

# Pathogen list from literature
PATHOGENS <- c(
  "Trueperella", "Fusobacterium", "Porphyromonas", "Ureaplasma",
  "Streptococcus", "Peptoniphilus", "Clostridium", "Staphylococcus",
  "Prevotella", "Campylobacter", "Helcococcus", "Bacteroides",
  "Peptoniphilus_A", "Porphyromonas_A", "Campylobacter_B"
)

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
# 2. Prepare data
# ============================================================================

message("\n========================================")
message("Step 2: Preparing data")
message("========================================\n")

# Aggregate to genus level
ps_genus <- tax_glom(ps_endo, taxrank = "Genus", NArm = FALSE)

# Get relative abundance
ps_ra <- transform_sample_counts(ps_genus, function(x) x/sum(x) * 100)

# Extract data
otu <- as(otu_table(ps_ra), "matrix")
if(!taxa_are_rows(ps_ra)) otu <- t(otu)
otu <- t(otu)  # samples as rows

meta <- as(sample_data(ps_ra), "data.frame")

# Get genus names
tax <- as.data.frame(tax_table(ps_ra))
colnames(otu) <- tax$Genus

# Combine with metadata
abund_df <- as.data.frame(otu)
abund_df$SampleID <- rownames(abund_df)
abund_df <- merge(abund_df, meta, by.x = "SampleID", by.y = "row.names")

message(paste("Data prepared:", nrow(abund_df), "samples"))

# ============================================================================
# 3. Correlation with CV scores (Figure 3A)
# ============================================================================

message("\n========================================")
message("Step 3: Correlation with CV scores")
message("========================================\n")

# Get genera columns
genera_cols <- colnames(otu)

# Calculate Spearman correlations
cor_results <- data.frame()

for(genus in genera_cols) {
  if(genus %in% colnames(abund_df) && "CV_Score" %in% colnames(abund_df)) {
    cor_test <- cor.test(abund_df[[genus]], abund_df$CV_Score, 
                         method = "spearman", exact = FALSE)
    
    cor_results <- rbind(cor_results, data.frame(
      Genus = genus,
      Rho = cor_test$estimate,
      P_value = cor_test$p.value,
      Mean_Abundance = mean(abund_df[[genus]], na.rm = TRUE)
    ))
  }
}

# FDR correction
cor_results$FDR <- p.adjust(cor_results$P_value, method = "BH")
cor_results$Significant <- cor_results$FDR < 0.05

# Sort by significance and correlation
cor_results <- cor_results %>%
  arrange(P_value)

# Save results
write.csv(cor_results, file.path(TABLE_PATH, "correlation_cv_score.csv"), row.names = FALSE)

# Print significant correlations
sig_cors <- cor_results %>% filter(Significant)
message(paste("\nSignificant correlations with CV score:", nrow(sig_cors)))
print(sig_cors[, c("Genus", "Rho", "P_value", "FDR")])

# ============================================================================
# 4. Visualize correlations
# ============================================================================

message("\n========================================")
message("Step 4: Visualizing correlations")
message("========================================\n")

# Top correlations plot
top_cor <- cor_results %>%
  filter(Mean_Abundance > 0.1) %>%  # Filter low abundance
  arrange(desc(abs(Rho))) %>%
  head(20)

p_cor <- ggplot(top_cor, aes(x = reorder(Genus, Rho), y = Rho, fill = Rho > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("Negative", "Positive")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Correlation with CV Score",
    subtitle = "Spearman correlation coefficient",
    x = "",
    y = "Spearman's ρ",
    fill = "Direction"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(file.path(FIG_PATH, "Fig3A_correlation_cv_score.pdf"), p_cor, width = 8, height = 8)
ggsave(file.path(FIG_PATH, "Fig3A_correlation_cv_score.png"), p_cor, width = 8, height = 8, dpi = 300)

message("Figure 3A saved")

# ============================================================================
# 5. Pathogen dynamics (Figure 3B)
# ============================================================================

message("\n========================================")
message("Step 5: Pathogen dynamics analysis")
message("========================================\n")

# Get pathogen abundances
pathogens_present <- intersect(PATHOGENS, colnames(abund_df))
message(paste("Pathogens found:", length(pathogens_present)))

if(length(pathogens_present) > 0) {
  # Calculate mean abundance by treatment and timepoint
  pathogen_dynamics <- abund_df %>%
    select(all_of(c("Treatment", "Timepoint", pathogens_present))) %>%
    pivot_longer(cols = all_of(pathogens_present), 
                 names_to = "Pathogen", 
                 values_to = "Abundance") %>%
    group_by(Treatment, Timepoint, Pathogen) %>%
    summarise(
      Mean = mean(Abundance, na.rm = TRUE),
      SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Calculate fold change from Day 0
  pathogen_fc <- pathogen_dynamics %>%
    pivot_wider(names_from = Timepoint, values_from = c(Mean, SE)) %>%
    mutate(
      FC_Day7 = ifelse(Mean_Day0 > 0, Mean_Day7 / Mean_Day0, NA),
      FC_Day14 = ifelse(Mean_Day0 > 0, Mean_Day14 / Mean_Day0, NA)
    )
  
  write.csv(pathogen_dynamics, file.path(TABLE_PATH, "pathogen_dynamics.csv"), row.names = FALSE)
  write.csv(pathogen_fc, file.path(TABLE_PATH, "pathogen_fold_change.csv"), row.names = FALSE)
  
  # Heatmap of pathogen changes
  pathogen_matrix <- pathogen_dynamics %>%
    filter(Timepoint == "Day14") %>%
    select(Treatment, Pathogen, Mean) %>%
    pivot_wider(names_from = Treatment, values_from = Mean)
  
  if(nrow(pathogen_matrix) > 0) {
    pm <- as.matrix(pathogen_matrix[, -1])
    rownames(pm) <- pathogen_matrix$Pathogen
    pm <- log10(pm + 0.01)
    
    pdf(file.path(FIG_PATH, "Fig3B_pathogen_heatmap.pdf"), width = 8, height = 10)
    pheatmap::pheatmap(
      pm,
      main = "Pathogen Abundance by Treatment (Day 14)",
      color = colorRampPalette(c("blue", "white", "red"))(50),
      cluster_cols = FALSE
    )
    dev.off()
    
    message("Figure 3B saved")
  }
}

# ============================================================================
# 6. Pregnancy outcome correlations (Figure 3G, 3H)
# ============================================================================

message("\n========================================")
message("Step 6: Pregnancy outcome correlations")
message("========================================\n")

if("Pregnancy_Status" %in% colnames(abund_df)) {
  # Separate pregnant and non-pregnant
  pregnant <- abund_df %>% filter(Pregnancy_Status == "Pregnant")
  non_pregnant <- abund_df %>% filter(Pregnancy_Status == "Non-pregnant")
  
  # Correlations for pregnant animals
  if(nrow(pregnant) >= 5) {
    cor_preg <- data.frame()
    for(pathogen in pathogens_present) {
      if("CV_Score" %in% colnames(pregnant)) {
        cor_test <- cor.test(pregnant[[pathogen]], pregnant$CV_Score, 
                             method = "spearman", exact = FALSE)
        cor_preg <- rbind(cor_preg, data.frame(
          Genus = pathogen,
          Rho = cor_test$estimate,
          P_value = cor_test$p.value
        ))
      }
    }
    cor_preg$FDR <- p.adjust(cor_preg$P_value, method = "BH")
    cor_preg$Group <- "Pregnant"
    write.csv(cor_preg, file.path(TABLE_PATH, "correlation_pregnant.csv"), row.names = FALSE)
  }
  
  # Correlations for non-pregnant animals
  if(nrow(non_pregnant) >= 5) {
    cor_nonpreg <- data.frame()
    for(pathogen in pathogens_present) {
      if("CV_Score" %in% colnames(non_pregnant)) {
        cor_test <- cor.test(non_pregnant[[pathogen]], non_pregnant$CV_Score, 
                             method = "spearman", exact = FALSE)
        cor_nonpreg <- rbind(cor_nonpreg, data.frame(
          Genus = pathogen,
          Rho = cor_test$estimate,
          P_value = cor_test$p.value
        ))
      }
    }
    cor_nonpreg$FDR <- p.adjust(cor_nonpreg$P_value, method = "BH")
    cor_nonpreg$Group <- "Non-pregnant"
    write.csv(cor_nonpreg, file.path(TABLE_PATH, "correlation_non_pregnant.csv"), row.names = FALSE)
  }
}

# ============================================================================
# 7. Correlation matrix visualization
# ============================================================================

message("\n========================================")
message("Step 7: Correlation matrix")
message("========================================\n")

# Genus-genus correlations for top abundant genera
top_genera <- cor_results %>%
  arrange(desc(Mean_Abundance)) %>%
  head(20) %>%
  pull(Genus)

if(length(top_genera) >= 5) {
  cor_matrix <- cor(abund_df[, top_genera], method = "spearman", use = "pairwise.complete.obs")
  
  pdf(file.path(FIG_PATH, "correlation_matrix_genera.pdf"), width = 10, height = 10)
  corrplot(
    cor_matrix,
    method = "color",
    type = "upper",
    order = "hclust",
    tl.col = "black",
    tl.srt = 45,
    addCoef.col = "black",
    number.cex = 0.6,
    title = "Genus-Genus Correlations",
    mar = c(0, 0, 2, 0)
  )
  dev.off()
  
  message("Correlation matrix saved")
}

message("\n========================================")
message("Correlation Analysis Complete!")
message("========================================\n")
