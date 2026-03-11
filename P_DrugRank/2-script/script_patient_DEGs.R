#----------------------------------------------------------------------------------
# Script: DEG_identification.R
# Author: Fiamma Romagnoli
#
# Purpose:
# Generates per-patient differential expression gene lists by comparing each
# TCGA tumor sample against a precomputed mean expression profile from normal
# control samples. For each tumor sample ID, the script:
#   1) Queries and prepares TCGA tumor count data (TCGAbiolinks)
#   2) Performs single-sample DEA vs normal mean profile
#   3) Maps identifiers to HGNC symbols
#   4) Writes outputs into a disease+drug-specific folder structure
#
# Configuration:
#   project - TCGA project ID (e.g., "TCGA-BRCA", "TCGA-COAD", "TCGA-READ", "TCGA-STAD")
#   disease - Disease label used in output folder naming
#   drug    - Drug label used to select cohort file and output folder
#   abb     - Cohort/controls abbreviation used in input filenames (must match on disk)
#
# Required inputs:
#   - Normal control mean expression matrix:
#       DrugRank/1-input_files/controls/Normal_control_<abb>/normal_ensid_count_<abb>.csv
#   - Patient cohort table with tumor sample IDs:
#       DrugRank/1-input_files/cohorts/<abb>_<drug>.csv   (expects column: sample_id)
#   - Functions sourced from:
#       DrugRank/2-script/DEG_functions.R
#     including: query_and_prepare_tumor(), perform_DEA(), map_to_hgnc()
#
# Outputs:
#   - Per-sample DEG results written to:
#       DrugRank/1-input_files/DEGs/<disease>_<drug>/<sample_id>/
#   - Log of samples with multiple tumor columns (or other flagged issues):
#       DrugRank/1-input_files/DEGs/<disease>_<drug>/samples_with_multiple_tumor_columns.txt
#
# Execution notes:
#   - Samples with existing output directories are skipped (idempotent behavior).
#   - Run from the project root so relative paths resolve correctly.
#   - Ensure the normal control file and cohort CSV exist for the chosen abb/drug.
#----------------------------------------------------------------------------------

#### Load Required Libraries ####
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)

source('2-script/DEG_functions.R')

### Modify based on analysis
project <- 'TCGA-BRCA'      # options: TCGA-COAD, TCGA READ (for CRC),  TCGA-STAD (for GC)
disease <- 'Breast_cancer'  # options: Colorectal_cancer, Stomach_cancer
drug <- 'fluorouracil'      # options: fluorouracil (for BRCA, CRC, GC); cyclophosphamide, docetaxel, doxorubicin_h, paclitaxel (for BRCA only); irinotecan (for CRC only).
abb <- 'BRCA'               # options: CRC_COAD, CRC_READ, GC

#### Define Output Directory ####
base_output_dir <- paste0("1-input_files/DEGs/", disease, "_" , drug, "/")
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

#### Mean Gene Expression in Normal Samples ####
mean_expression_normal_samples <- read.csv(paste0('1-input_files/controls/Normal_control_', abb, '/normal_ensid_count_', abb, '.csv'), row.names = 1)
normal_data <- data.matrix(mean_expression_normal_samples)


########################## Main Execution ##################################

#### Load file with patient data ####

patients <- read.csv(paste0("1-input_files/cohorts/", abb, "_", drug, ".csv"), header = TRUE, sep = ';')

# Select all sample ids
tumor_sample_ids <- patients$sample_id
total <- length(tumor_sample_ids)  # Get the total number of samples
counter <- 1 # Initialize counter
problematic_samples_log <- file.path(base_output_dir, "samples_with_multiple_tumor_columns.txt")

for (i in seq_along(tumor_sample_ids)) {
  id <- tumor_sample_ids[i]
  
  # Build output directory
  output_dir <- paste0(base_output_dir, id, "/")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  } else {
    cat("⏭️  Skipping", id, "- already processed.\n")
    next
  }
  
  cat(paste0("🔄 [", i, "/", length(tumor_sample_ids), "] Processing ID: ", id, "\n"))
  
  
  # Run analysis
  tumor_data <- query_and_prepare_tumor(
    project, 
    "Transcriptome Profiling", 
    "Gene Expression Quantification", 
    "STAR - Counts", 
    id
    )
  
  degs <- perform_DEA(
    tumor_data, 
    normal_data, 
    id, 
    output_dir,
    log_path = problematic_samples_log
    )
  
  map_to_hgnc(
    degs, 
    id, 
    output_dir
    )
    
  cat("✅ Analysis complete. Results saved in", output_dir, "\n\n")
  counter <- counter + 1  
}


