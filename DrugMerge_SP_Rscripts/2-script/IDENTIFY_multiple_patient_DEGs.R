#### Load Required Libraries ####
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)

source('DrugMerge_SP_Final/2-script/DEG_functions.R')

#### Define Output Directory ####
base_output_dir <- "DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/"

#### Mean Gene Expression in Normal Samples (n=41) ####
mean_expression_normal_samples <- read.csv('DrugMerge_SP_Final/1-input_files/normal_ensid_count.csv', row.names = 1)
normal_data <- data.matrix(mean_expression_normal_samples)


######################## Main Execution ##################################

#### Define project and load file with patient data ####
project <- "TCGA-COAD"
patients <- read.csv("DrugMerge_SP_Final/1-input_files/updated_patient_cohort.csv", header = TRUE, sep = ';')

# Select all sample ids
tumor_sample_ids <- patients$Sample.ID
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
  tumor_data <- query_and_prepare_tumor(project, "Transcriptome Profiling", "Gene Expression Quantification", "STAR - Counts", id)
  
  degs <- perform_DEA(tumor_data, normal_data, id, output_dir,log_path = problematic_samples_log)
  map_to_hgnc(degs,id, output_dir)
  
  
  cat("✅ Analysis complete. Results saved in", output_dir, "\n\n")
  counter <- counter + 1  # Increment counter  
}


