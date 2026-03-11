#------------------------------------------------------------------------------
# Script: Personalised-DrugRank - script_colorectalcancer.R
# Author: Fiamma Romagnoli
#
# Purpose:
# Runs patient-level drug repurposing analysis in batches for a selected cancer
# type and treatment drug. For each patient/sample it:
#   1) Loads UP/DOWN DEG gene lists
#   2) Runs single-module Enrichr enrichment (e.g., Core&Peel or MD)
#   3) Runs multiple-module Enrichr enrichment (per detected module)
#   4) Builds an unfiltered drug–module modulation table
#   5) Adds P.value / Adjusted.P.value from module-level Enrichr outputs
#   6) Filters modulation results for direction-consistent drug signatures
#
# Configuration:
#   - Cancer type (e.g., "Stomach_cancer", "Breast_cancer", "Colorectal_cancer")
#   - Drug name
#   - Enrichr databases, in our case Drug Perturbations from GEO 2014
#   - Module detection method and parameters
#   - Batch size and batch number
#
# Required inputs:
#   - DEG directories:
#       1-input_files/DEGs/<Cancer>_<drug>/<patient_dir>/
#         * <sample>_UP_gene_list.txt
#         * <sample>_DOWN_gene_list.txt
#   - Active subnetworks:
#       1-input_files/active_subnetworks/<Cancer>/
#         * Single-module file
#         * MultipleModules/ file
#   - Gene annotation map
#
# Outputs:
#   - Enrichr results per sample under:
#       3-Drugs/<Cancer_short>/<Cancer_short>_<drug>/<sample_id>/
#   - Multiple-module enrichment:
#       4-MultipleModules_output/<Cancer>/<sample_id>/<method>/<database>/
#   - Drug modulation tables (unfiltered and filtered):
#       <sample_id>/Drug-modulation-data/<cancer-specific-subfolder>/
#
# Execution:
#   - Set parameters in cfg
#   - Run from project root so relative paths resolve correctly
#   - Logs are written per batch under logs/<Cancer>/<drug>/
#
# Notes:
#   - Database naming must match folder structure exactly.
#   - Module file naming must uniquely match the selected method.
#------------------------------------------------------------------------------

library(limma)
library(enrichR)
library(readr)

source('2-script/functions.R')

# ----------------------------- Config ---------------------------------------
cfg <- list(
  cancer                  = "Colorectal_cancer",
  drug                    = "irinotecan",      # options: fluorouracil, irinotecan
  base_deg_dir            = NULL,              # filled below
  batch_size              = 5,
  batch_number            = 4,

  gene_map_path           = "1-input_files/symbol_entrezID_updated_table.txt",

  # Output / input roots
  drugs_root              = "3-Drugs/CRC",
  subnet_root             = "1-input_files/active_subnetworks/Colorectal_cancer/",
  multiple_modules_out    = "4-MultipleModules_output",

  # Databases
  db_single               = c("Drug_Perturbations_from_GEO_2014"),
  db_modulation           = c("drugperturbationGEO"),

  # Method metadata
  method_single           = "Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest",
  method_label_modulation = "Core&Peel-d0.7-j0.8",
  param_cp                = "-d0.7-f1-j0.8",

  # Filenames
  file_single             = "CRC_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest.txt",
  file_multiple_modules   = "CRC_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt"
)

cfg$base_deg_dir <- file.path("1-input_files/DEGs/", paste0("Colorectal_cancer_", cfg$drug))


# ----------------------------- Setup ----------------------------------------

all_patients <- list.dirs(cfg$base_deg_dir, full.names = TRUE, recursive = FALSE)
if (length(all_patients) == 0) stop("No patient directories found in: ", cfg$base_deg_dir)

batches <- split(all_patients, ceiling(seq_along(all_patients) / cfg$batch_size))

if (cfg$batch_number < 1 || cfg$batch_number > length(batches)) {
  stop("batch_number out of range. Got ", cfg$batch_number, " but batches are 1..", length(batches))
}

logfile <- file.path("logs", cfg$cancer, cfg$drug, sprintf("log_batch_%02d.txt", cfg$batch_number))
log_con <- start_log(logfile)
on.exit(end_log(log_con), add = TRUE)

cat("🟢 Starting", cfg$cancer, "batch", cfg$batch_number, "\n")
batch_start_time <- Sys.time()

gene_annotation_map <- read.delim(cfg$gene_map_path, header = TRUE, stringsAsFactors = FALSE)

# ----------------------------- Main -----------------------------------------
process_sample <- function(patient_path, cfg, gene_annotation_map) {
  folder_name <- basename(patient_path)
  sample_id <- extract_sample_id(folder_name)

  cat("\n🧬 Starting analysis for:", sample_id, "\n")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Start time for", sample_id, "\n")
  start_time <- Sys.time()
  
  # Load DEG lists
  up_path   <- file.path(patient_path, paste0(sample_id, "_UP_gene_list.txt"))
  down_path <- file.path(patient_path, paste0(sample_id, "_DOWN_gene_list.txt"))

  if (!file.exists(up_path) || !file.exists(down_path)) {
    cat("⚠️ Skipping", sample_id, "- missing DEG files.\n")
    return(invisible(NULL))
  }

  up <- read_gene_list(up_path)
  down <- read_gene_list(down_path)

  
  # Output directories
  output_prefix <- paste0('3-Drugs/CRC/CRC_', cfg$drug, '/', sample_id, '/')
  create_directory_GEO_CRC(paste0('3-Drugs/CRC/CRC_', cfg$drug, '/'), sample_id)
  
  
  # -------- Single-module enrichment (Core&Peel) ----------
  cat("🧪 Core&Peel (single module)\n")
  enrichr_analysis(
    cfg$subnet_root,
    cfg$file_single,
    cfg$method_single,
    "Colorectal_cancer",
    cfg$db_single,
    gene_annotation_map,
    up,
    down,
    output_prefix
  )
  cat("Finished single-module\n")
  
  cat("------------------------------------------------------------\n")
  
  # -------- Multiple-modules enrichment ----------

  # Create output directory for multiple modules
  safe_dir_create(cfg$multiple_modules_out)
    
  # 1) perform the enrichr analysis on the single modules
  
  cat('🧬Starting multiple-module enrichment analysis for Core&Peel \n')
  
  mm_file <- file.path(cfg$subnet_root, "MultipleModules", cfg$file_multiple_modules)

  for (d in cfg$db_single) {
    enrichr_analysis_multipleModules(mm_file, "Colorectal_cancer", cfg$method_label_modulation, d, sample_id)
  }
  cat("✅ Multiple-module enrichr done\n")

  cat("------------------------------------------------------------\n")
  
  # -------- Modulation table (no filtering + filtering) ----------

  ##### 2) get table without filtering #########

  cat("📝Creating table without filtering for Core&Peel\n")

  path <- "1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules"
    print(path)
    print(list.files(path, pattern="\\.txt$", full.names=TRUE))

  for (d in cfg$db_modulation) {

    enr <- read_enrichr(
      file.path(cfg$drugs_root, paste0("CRC_", cfg$drug, '/', sample_id, '/')),
      "Core&Peel",
      "Colorectal_cancer",
      d,
      param = cfg$param_cp
    )

    mm_path <- file.path(cfg$subnet_root, "MultipleModules")
    mod_table <- get_modulation_table(enr, "Core&Peel", mm_path, up, down)

    path_input <- file.path(cfg$multiple_modules_out, "Colorectal_cancer", sample_id)
    path_output <- file.path(cfg$drugs_root,  paste0("CRC_", cfg$drug, '/', sample_id),
                             "Drug-modulation-data/colorectal-cancer-cp")
    safe_dir_create(path_output)

    add_pvalue(path_input, "Core&Peel", d, mod_table, path_output, "Colorectal_cancer", param = cfg$param_cp)

    d1 <- switch(d,
      "drugperturbationGEO" = "GEO",
      "Drug_Perturbations_from_GEO_2014" = "GEO",
      d
    )

    csv_path <- file.path(
      path_output,
      paste0("Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_", d1, "_noFiltering.csv")
    )

    if (file.exists(csv_path)) {
      table_nf <- read.csv(csv_path)
      get_modulation_filter(table_nf, d, "Core&Peel", path_output, "Colorectal_cancer", param = cfg$param_cp)
    } else {
      cat("⚠️ Expected CSV not found:", csv_path, "\n")
    }
  }
  cat("Done\n")
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Finished analysis for", sample_id, "\n")
  cat("⏱️ Duration:", round(duration, 2), "minutes\n")
  cat("------------------------------------------------------------\n")
}

for (patient_path in batches[[cfg$batch_number]]) {
  tryCatch(
    process_sample(patient_path, cfg, gene_annotation_map),
    error = function(e) {
      folder_name <- basename(patient_path)
      sample_id <- extract_sample_id(folder_name)
      cat("❌ Error in sample", sample_id, ":", conditionMessage(e), "\n")
    }
  )
}

batch_end_time <- Sys.time()
batch_duration <- difftime(batch_end_time, batch_start_time, units = "hours")
cat("\n✅ Finished batch", cfg$batch_number, "at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("⏱️ Total batch duration:", round(batch_duration, 2), "hours\n")  


