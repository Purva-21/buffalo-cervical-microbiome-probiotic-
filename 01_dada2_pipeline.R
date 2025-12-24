# ============================================================================
# 01_dada2_pipeline.R
# DADA2 Processing Pipeline for 16S rRNA Amplicon Sequencing
# Buffalo Cervical Microbiome Study
# ============================================================================

# Load required libraries
library(dada2)
library(Biostrings)
library(ShortRead)
library(ggplot2)

# ============================================================================
# Configuration
# ============================================================================

# Set paths
RAW_DATA_PATH <- "data/raw/"
OUTPUT_PATH <- "data/processed/"
RESULTS_PATH <- "results/"

# Create output directories if they don't exist
dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_PATH, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_PATH, "tables"), recursive = TRUE, showWarnings = FALSE)

# Sequencing parameters (Ion Torrent 530 chip, 400bp chemistry)
TRUNC_LEN <- 350        # Truncation length
MAX_EE <- 2             # Maximum expected errors
TRUNC_Q <- 2            # Truncate at first quality score <= this
MIN_LEN <- 200          # Minimum read length after filtering

# Taxonomy database path (NCBI RefSeq 16S)
TAX_DB_PATH <- "data/reference/ncbi_refseq_16S.fa.gz"

# ============================================================================
# 1. List and inspect input files
# ============================================================================

message("\n========================================")
message("Step 1: Inspecting input files")
message("========================================\n")

# List fastq files
fastq_files <- sort(list.files(RAW_DATA_PATH, pattern = ".fastq", full.names = TRUE))
message(paste("Found", length(fastq_files), "fastq files"))

# Extract sample names
sample_names <- sapply(strsplit(basename(fastq_files), "_"), `[`, 1)
message(paste("Sample names:", paste(head(sample_names), collapse = ", "), "..."))

# ============================================================================
# 2. Quality profile inspection
# ============================================================================

message("\n========================================")
message("Step 2: Generating quality profiles")
message("========================================\n")

# Plot quality profiles for first 4 samples
pdf(file.path(RESULTS_PATH, "figures", "quality_profiles.pdf"), width = 12, height = 8)
print(plotQualityProfile(fastq_files[1:4]))
dev.off()

message("Quality profiles saved to results/figures/quality_profiles.pdf")

# ============================================================================
# 3. Filter and trim
# ============================================================================

message("\n========================================")
message("Step 3: Filtering and trimming reads")
message("========================================\n")

# Create filtered file paths
filt_path <- file.path(OUTPUT_PATH, "filtered")
dir.create(filt_path, showWarnings = FALSE)
filt_files <- file.path(filt_path, paste0(sample_names, "_filt.fastq.gz"))
names(filt_files) <- sample_names

# Filter and trim
out <- filterAndTrim(
  fastq_files, filt_files,
  truncLen = TRUNC_LEN,
  maxEE = MAX_EE,
  truncQ = TRUNC_Q,
  minLen = MIN_LEN,
  maxN = 0,
  compress = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

# Save filtering stats
write.csv(out, file.path(RESULTS_PATH, "tables", "filtering_stats.csv"))

# Print summary
message("\nFiltering summary:")
message(paste("  - Reads in:", sum(out[,1])))
message(paste("  - Reads out:", sum(out[,2])))
message(paste("  - Retention rate:", round(sum(out[,2])/sum(out[,1])*100, 2), "%"))

# ============================================================================
# 4. Learn error rates
# ============================================================================

message("\n========================================")
message("Step 4: Learning error rates")
message("========================================\n")

# Learn errors
err <- learnErrors(filt_files, multithread = TRUE, verbose = TRUE)

# Plot error rates
pdf(file.path(RESULTS_PATH, "figures", "error_rates.pdf"), width = 10, height = 8)
print(plotErrors(err, nominalQ = TRUE))
dev.off()

message("Error rate plot saved to results/figures/error_rates.pdf")

# ============================================================================
# 5. Dereplication
# ============================================================================

message("\n========================================")
message("Step 5: Dereplicating sequences")
message("========================================\n")

derep <- derepFastq(filt_files, verbose = TRUE)
names(derep) <- sample_names

# ============================================================================
# 6. Sample inference (ASV calling)
# ============================================================================

message("\n========================================")
message("Step 6: Inferring ASVs")
message("========================================\n")

dada_res <- dada(derep, err = err, multithread = TRUE, verbose = TRUE)

# Print summary for first sample
message("\nASV inference summary (first sample):")
print(dada_res[[1]])

# ============================================================================
# 7. Construct sequence table
# ============================================================================

message("\n========================================")
message("Step 7: Constructing sequence table")
message("========================================\n")

seqtab <- makeSequenceTable(dada_res)

message(paste("Sequence table dimensions:", nrow(seqtab), "samples x", ncol(seqtab), "ASVs"))

# Check sequence length distribution
seq_lengths <- nchar(getSequences(seqtab))
message("\nSequence length distribution:")
print(table(seq_lengths))

# ============================================================================
# 8. Remove chimeras
# ============================================================================

message("\n========================================")
message("Step 8: Removing chimeric sequences")
message("========================================\n")

seqtab_nochim <- removeBimeraDenovo(
  seqtab, 
  method = "consensus", 
  multithread = TRUE, 
  verbose = TRUE
)

message(paste("\nASVs before chimera removal:", ncol(seqtab)))
message(paste("ASVs after chimera removal:", ncol(seqtab_nochim)))
message(paste("Chimeras removed:", ncol(seqtab) - ncol(seqtab_nochim)))
message(paste("Reads retained:", round(sum(seqtab_nochim)/sum(seqtab)*100, 2), "%"))

# ============================================================================
# 9. Track reads through pipeline
# ============================================================================

message("\n========================================")
message("Step 9: Generating pipeline summary")
message("========================================\n")

# Get counts at each step
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dada_res, getN),
  rowSums(seqtab_nochim)
)
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample_names

# Save tracking table
write.csv(track, file.path(RESULTS_PATH, "tables", "read_tracking.csv"))

message("Read tracking summary:")
print(head(track))

# ============================================================================
# 10. Assign taxonomy
# ============================================================================

message("\n========================================")
message("Step 10: Assigning taxonomy")
message("========================================\n")

# Check if taxonomy database exists
if (file.exists(TAX_DB_PATH)) {
  taxa <- assignTaxonomy(
    seqtab_nochim, 
    TAX_DB_PATH,
    multithread = TRUE,
    verbose = TRUE
  )
  
  # Remove sequence rownames for display
  taxa_print <- taxa
  rownames(taxa_print) <- NULL
  
  message("\nTaxonomy assignment summary:")
  print(head(taxa_print))
} else {
  message("WARNING: Taxonomy database not found at: ", TAX_DB_PATH)
  message("Please download the NCBI RefSeq 16S database and update the path.")
  taxa <- NULL
}

# ============================================================================
# 11. Save outputs
# ============================================================================

message("\n========================================")
message("Step 11: Saving output files")
message("========================================\n")

# Save ASV table
saveRDS(seqtab_nochim, file.path(OUTPUT_PATH, "asv_table.rds"))

# Save as CSV (with ASV IDs instead of sequences)
asv_seqs <- colnames(seqtab_nochim)
asv_ids <- paste0("ASV", seq_along(asv_seqs))
colnames(seqtab_nochim) <- asv_ids

asv_table_export <- as.data.frame(seqtab_nochim)
asv_table_export$SampleID <- rownames(asv_table_export)
write.csv(asv_table_export, file.path(OUTPUT_PATH, "asv_counts.csv"), row.names = FALSE)

# Save ASV sequences
asv_seqs_df <- data.frame(
  ASV_ID = asv_ids,
  Sequence = asv_seqs
)
write.csv(asv_seqs_df, file.path(OUTPUT_PATH, "asv_sequences.csv"), row.names = FALSE)

# Export as FASTA
asv_fasta <- DNAStringSet(asv_seqs)
names(asv_fasta) <- asv_ids
writeXStringSet(asv_fasta, file.path(OUTPUT_PATH, "asv_sequences.fasta"))

# Save taxonomy table
if (!is.null(taxa)) {
  taxa_export <- as.data.frame(taxa)
  taxa_export$ASV_ID <- asv_ids
  taxa_export$Sequence <- asv_seqs
  write.csv(taxa_export, file.path(OUTPUT_PATH, "taxonomy_table.csv"), row.names = FALSE)
  saveRDS(taxa, file.path(OUTPUT_PATH, "taxonomy.rds"))
}

# ============================================================================
# Summary
# ============================================================================

message("\n========================================")
message("DADA2 Pipeline Complete!")
message("========================================\n")

message("Output files saved to: ", OUTPUT_PATH)
message("  - asv_table.rds: ASV table (RDS format)")
message("  - asv_counts.csv: ASV abundance table")
message("  - asv_sequences.csv: ASV sequences")
message("  - asv_sequences.fasta: ASV sequences (FASTA)")
message("  - taxonomy_table.csv: Taxonomy assignments")
message("  - taxonomy.rds: Taxonomy (RDS format)")

message("\nResults saved to: ", RESULTS_PATH)
message("  - figures/quality_profiles.pdf")
message("  - figures/error_rates.pdf")
message("  - tables/filtering_stats.csv")
message("  - tables/read_tracking.csv")

message("\nPipeline completed successfully!")
