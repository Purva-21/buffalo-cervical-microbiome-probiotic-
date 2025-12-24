# buffalo-cervical-microbiome-probiotic-
# Probiotic-Driven Remodeling of the Cervical Microbiome in Postpartum *Bubalus bubalis*

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-4.3.2-blue.svg)](https://www.r-project.org/)

## Overview

This repository contains the analysis scripts and data processing pipelines for the study:

**"Microbial Refrains of Recovery: Probiotic-Driven Remodeling of the Cervical Microbiome in Postpartum *Bubalus bubalis*"**

The study evaluates the therapeutic efficacy of intrauterine probiotic administration (*Lactiplantibacillus plantarum* KUGBRC [LP] and *Pediococcus pentosaceus* GBRCKU [MS]) in buffaloes with endometritis and/or anestrus, using 16S rRNA amplicon sequencing.

## Study Design

- **Total samples**: 164 postpartum buffaloes
- **Groups**: Healthy (n=38), Endometritis (n=91), Anestrus (n=35)
- **Treatments**: 
  - LP (*Lactiplantibacillus plantarum* KUGBRC): n=44
  - MS (*Pediococcus pentosaceus* GBRCKU): n=26
  - NS (Normal Saline - control): n=52
- **Sampling timepoints**: Day 0, Day 7, Day 14

## Repository Structure

```
âââ README.md
âââ LICENSE
âââ data/
â   âââ raw/                    # Raw sequencing data (not included, see Data Availability)
â   âââ processed/              # Processed ASV tables and taxonomy
â   âââ metadata/               # Sample metadata files
âââ scripts/
â   âââ 01_dada2_pipeline.R     # DADA2 processing pipeline
â   âââ 02_data_preparation.R   # Data import and preparation
â   âââ 03_alpha_diversity.R    # Alpha diversity analysis
â   âââ 04_beta_diversity.R     # Beta diversity and ordination
â   âââ 05_differential_abundance.R  # DESeq2 and LEfSe analysis
â   âââ 06_simper_analysis.R    # SIMPER community dissimilarity
â   âââ 07_correlation_analysis.R    # Pathogen-severity correlations
â   âââ 08_reproductive_outcomes.R   # Estrus and pregnancy analysis
â   âââ utils/
â       âââ helper_functions.R  # Custom utility functions
âââ results/
â   âââ figures/                # Generated figures
â   âââ tables/                 # Statistical output tables
âââ environment/
    âââ requirements.R          # R package dependencies
```

## Installation

### Prerequisites

- R (>= 4.3.2)
- RStudio (recommended)

### Install Required Packages

```r
# Run this script to install all dependencies
source("environment/requirements.R")
```

## Quick Start

```r
# 1. Clone the repository
# git clone https://github.com/YOUR_USERNAME/buffalo-cervical-microbiome.git

# 2. Set working directory
setwd("path/to/buffalo-cervical-microbiome")

# 3. Install dependencies
source("environment/requirements.R")

# 4. Run the analysis pipeline
source("scripts/01_dada2_pipeline.R")      # Process raw reads
source("scripts/02_data_preparation.R")     # Prepare data objects
source("scripts/03_alpha_diversity.R")      # Alpha diversity
source("scripts/04_beta_diversity.R")       # Beta diversity
source("scripts/05_differential_abundance.R") # Differential abundance
source("scripts/06_simper_analysis.R")      # SIMPER analysis
source("scripts/07_correlation_analysis.R") # Correlations
source("scripts/08_reproductive_outcomes.R") # Clinical outcomes
```

## Analysis Overview

### 1. Sequence Processing (DADA2)
- Quality filtering and trimming
- Denoising and dereplication
- Chimera removal
- ASV inference
- Taxonomic assignment using NCBI RefSeq 16S database

### 2. Diversity Analysis
- **Alpha diversity**: Chao1, Shannon, Simpson indices
- **Beta diversity**: Bray-Curtis dissimilarity, PCoA ordination
- **Statistical tests**: Kruskal-Wallis, PERMANOVA (Adonis)

### 3. Differential Abundance
- DESeq2 for differential abundance testing
- LEfSe (LDA Effect Size) for biomarker discovery
- Benjamini-Hochberg FDR correction

### 4. Community Analysis
- SIMPER analysis for community dissimilarity drivers
- Venn diagrams for shared/unique taxa
- Correlation networks

### 5. Clinical Outcomes
- Estrus induction rates
- Pregnancy outcomes
- Pathogen-severity associations

## Key Findings

1. **Endometritis** buffaloes exhibit significantly higher microbial diversity and distinct community composition compared to healthy animals

2. **LP treatment** achieved:
   - Highest estrus induction rate (81.6%)
   - Highest pregnancy rate (59.0%)
   - Broadest pathogen suppression

3. **MS treatment** induced extensive community turnover (73% unique ASVs) but with less favorable reproductive outcomes

4. **Pregnant buffaloes** showed reduced microbial diversity and more homogeneous community structure

## Abbreviations

| Abbreviation | Full Name |
|-------------|-----------|
| LP | *Lactiplantibacillus plantarum* KUGBRC |
| MS | *Pediococcus pentosaceus* GBRCKU |
| NS | Normal Saline |
| ASV | Amplicon Sequence Variant |
| CV | Cervico-vaginal (score) |
| PCoA | Principal Coordinate Analysis |
| PERMANOVA | Permutational Multivariate Analysis of Variance |

## Data Availability

Raw sequencing data have been deposited in NCBI SRA under BioProject accession number: **[PENDING]**
## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Gujarat Biotechnology Research Centre (GBRC)
- Sanoda Dehgam Animal Ambulatory Hospital
- All participating dairy farms


