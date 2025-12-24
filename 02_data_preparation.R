# ============================================================================
# 02_data_preparation.R
# Data Import and Preparation for Microbiome Analysis
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(tidyverse)
library(microeco)
library(phyloseq)
library(vegan)

# ============================================================================
# Configuration
# ============================================================================

# Paths
PROCESSED_PATH <- "data/processed/"
METADATA_PATH <- "data/metadata/"
OUTPUT_PATH <- "data/processed/"
RESULTS_PATH <- "results/"

# Minimum read depth for sample inclusion
MIN_READS <- 10290

# Random seed for reproducibility
set.seed(42)

# ============================================================================
# 1. Load data
# ============================================================================

message("\n========================================")
message("Step 1: Loading data files")
message("========================================\n")

# Load ASV table
asv_table <- readRDS(file.path(PROCESSED_PATH, "asv_table.rds"))
message(paste("ASV table loaded:", nrow(asv_table), "samples x", ncol(asv_table), "ASVs"))

# Load taxonomy
taxonomy <- readRDS(file.path(PROCESSED_PATH, "taxonomy.rds"))
message(paste("Taxonomy table loaded:", nrow(taxonomy), "ASVs"))

# Load metadata
metadata <- read.csv(file.path(METADATA_PATH, "sample_metadata.csv"), row.names = 1)
message(paste("Metadata loaded:", nrow(metadata), "samples"))

# ============================================================================
# 2. Prepare metadata
# ============================================================================

message("\n========================================")
message("Step 2: Preparing metadata")
message("========================================\n")

# Define group labels
# LP = Lactiplantibacillus plantarum KUGBRC
# MS = Pediococcus pentosaceus GBRCKU  
# NS = Normal Saline

# Ensure consistent factor levels
metadata$Group <- factor(metadata$Group, 
                         levels = c("Healthy", "Endometritis", "Anestrus"))

metadata$Treatment <- factor(metadata$Treatment,
                             levels = c("Untreated", "LP", "MS", "NS"))

metadata$Timepoint <- factor(metadata$Timepoint,
                             levels = c("Day0", "Day7", "Day14"))

# Print group summary
message("\nGroup distribution:")
print(table(metadata$Group))

message("\nTreatment distribution:")
print(table(metadata$Treatment))

# ============================================================================
# 3. Quality filtering
# ============================================================================

message("\n========================================")
message("Step 3: Quality filtering samples")
message("========================================\n")

# Calculate read depths
read_depths <- rowSums(asv_table)

message(paste("Read depth range:", min(read_depths), "-", max(read_depths)))
message(paste("Median read depth:", median(read_depths)))

# Filter samples by minimum reads
samples_pass <- names(read_depths[read_depths >= MIN_READS])
samples_fail <- names(read_depths[read_depths < MIN_READS])

message(paste("\nSamples passing filter (>=", MIN_READS, "reads):", length(samples_pass)))
message(paste("Samples failing filter:", length(samples_fail)))

if(length(samples_fail) > 0) {
  message("Excluded samples: ", paste(samples_fail, collapse = ", "))
}

# Filter data
asv_table_filt <- asv_table[samples_pass, ]
metadata_filt <- metadata[samples_pass, ]

# Remove zero-sum ASVs
asv_sums <- colSums(asv_table_filt)
asv_table_filt <- asv_table_filt[, asv_sums > 0]

message(paste("\nFiltered data: ", nrow(asv_table_filt), "samples x", ncol(asv_table_filt), "ASVs"))

# ============================================================================
# 4. Create microeco object
# ============================================================================

message("\n========================================")
message("Step 4: Creating microeco object")
message("========================================\n")

# Prepare taxonomy dataframe
tax_df <- as.data.frame(taxonomy)
colnames(tax_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter taxonomy to match ASVs
common_asvs <- intersect(colnames(asv_table_filt), rownames(tax_df))
asv_table_filt <- asv_table_filt[, common_asvs]
tax_df <- tax_df[common_asvs, ]

# Create microeco dataset
meco <- microtable$new(
  sample_table = metadata_filt,
  otu_table = as.data.frame(t(asv_table_filt)),  # microeco expects ASVs as rows
  tax_table = tax_df
)

# Print summary
message("\nMicroeco object created:")
meco$print_taxa()

# ============================================================================
# 5. Create phyloseq object (alternative format)
# ============================================================================

message("\n========================================")
message("Step 5: Creating phyloseq object")
message("========================================\n")

# Create phyloseq components
OTU <- otu_table(as.matrix(asv_table_filt), taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(tax_df))
SAMP <- sample_data(metadata_filt)

# Create phyloseq object
ps <- phyloseq(OTU, TAX, SAMP)

message("Phyloseq object created:")
print(ps)

# ============================================================================
# 6. Subset data by analysis groups
# ============================================================================

message("\n========================================")
message("Step 6: Creating analysis subsets")
message("========================================\n")

# Subset 1: Healthy vs Endometritis vs Anestrus comparison
ps_health <- subset_samples(ps, Group %in% c("Healthy", "Endometritis", "Anestrus"))
message(paste("Health comparison subset:", nsamples(ps_health), "samples"))

# Subset 2: Endometritis treatment comparison
ps_endo_treat <- subset_samples(ps, Group == "Endometritis")
message(paste("Endometritis treatment subset:", nsamples(ps_endo_treat), "samples"))

# Subset 3: Anestrus treatment comparison  
ps_anestrus_treat <- subset_samples(ps, Group == "Anestrus")
message(paste("Anestrus treatment subset:", nsamples(ps_anestrus_treat), "samples"))

# Subset 4: Pregnancy outcome comparison
ps_pregnancy <- subset_samples(ps, !is.na(Pregnancy_Status))
message(paste("Pregnancy outcome subset:", nsamples(ps_pregnancy), "samples"))

# ============================================================================
# 7. Normalize data
# ============================================================================

message("\n========================================")
message("Step 7: Data normalization")
message("========================================\n")

# Rarefaction (for alpha diversity)
min_depth <- min(sample_sums(ps_health))
message(paste("Rarefying to minimum depth:", min_depth))

set.seed(42)
ps_rare <- rarefy_even_depth(ps_health, sample.size = min_depth, verbose = FALSE)
message(paste("Rarefied samples:", nsamples(ps_rare)))

# Relative abundance transformation (for compositional analysis)
ps_relabund <- transform_sample_counts(ps_health, function(x) x / sum(x) * 100)

# CLR transformation (for differential abundance)
# Using microeco
meco$cal_abund()

# ============================================================================
# 8. Save processed objects
# ============================================================================

message("\n========================================")
message("Step 8: Saving processed data")
message("========================================\n")

# Save microeco object
saveRDS(meco, file.path(OUTPUT_PATH, "meco_object.rds"))

# Save phyloseq objects
saveRDS(ps, file.path(OUTPUT_PATH, "phyloseq_full.rds"))
saveRDS(ps_health, file.path(OUTPUT_PATH, "phyloseq_health.rds"))
saveRDS(ps_endo_treat, file.path(OUTPUT_PATH, "phyloseq_endometritis.rds"))
saveRDS(ps_anestrus_treat, file.path(OUTPUT_PATH, "phyloseq_anestrus.rds"))
saveRDS(ps_rare, file.path(OUTPUT_PATH, "phyloseq_rarefied.rds"))
saveRDS(ps_relabund, file.path(OUTPUT_PATH, "phyloseq_relabund.rds"))

# Save filtered metadata
write.csv(metadata_filt, file.path(METADATA_PATH, "metadata_filtered.csv"))

message("Saved objects:")
message("  - meco_object.rds")
message("  - phyloseq_full.rds")
message("  - phyloseq_health.rds")
message("  - phyloseq_endometritis.rds")
message("  - phyloseq_anestrus.rds")
message("  - phyloseq_rarefied.rds")
message("  - phyloseq_relabund.rds")

# ============================================================================
# 9. Generate summary statistics
# ============================================================================

message("\n========================================")
message("Step 9: Summary statistics")
message("========================================\n")

# Sample summary by group
sample_summary <- metadata_filt %>%
  group_by(Group, Treatment) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

message("Sample distribution:")
print(sample_summary)

# Taxonomic summary
tax_summary <- data.frame(
  Level = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  Total = sapply(tax_df, function(x) length(unique(na.omit(x))))
)

message("\nTaxonomic diversity:")
print(tax_summary)

# Save summaries
write.csv(sample_summary, file.path(RESULTS_PATH, "tables", "sample_summary.csv"), row.names = FALSE)
write.csv(tax_summary, file.path(RESULTS_PATH, "tables", "taxonomy_summary.csv"), row.names = FALSE)

message("\n========================================")
message("Data Preparation Complete!")
message("========================================\n")
