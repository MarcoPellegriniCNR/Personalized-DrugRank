##################################################################################################################
#                                                   DEG Functions
##################################################################################################################

#-----------------------------------------------------------------------------------------------------------------
# Function: query_and_prepare_tumor
#
# Description:
# Queries, downloads, and prepares TCGA tumor expression data for a single
# sample using TCGAbiolinks, with automatic retry logic. If a download or
# preparation step fails, the function retries up to max_retries times before
# stopping with an error.
#
# Parameters:
#   project        - TCGA project ID (e.g., "TCGA-BRCA", "TCGA-COAD", "TCGA-STAD")
#   data_category  - GDC data category (e.g., "Transcriptome Profiling")
#   data_type      - GDC data type (e.g., "Gene Expression Quantification")
#   workflow_type  - Workflow type (e.g., "STAR - Counts")
#   barcode        - TCGA sample barcode
#   max_retries    - Maximum number of download/preparation attempts (default = 3)
#
# Output:
#   A SummarizedExperiment object containing tumor expression data.
#   Stops with an error if all retry attempts fail.
#-----------------------------------------------------------------------------------------------------------------

query_and_prepare_tumor <- function(project, data_category, data_type, workflow_type, barcode, max_retries = 3) {
  attempt <- 1
  success <- FALSE
  while (attempt <= max_retries && !success) {
    cat("🧪 Attempt", attempt, "for sample:", barcode, "\n")
    tryCatch({
      query <- GDCquery(
        project = project,
        data.category = data_category,
        data.type = data_type,
        workflow.type = workflow_type,
        barcode =  barcode
      )
      GDCdownload(query)
      data <- GDCprepare(query)
      success <- TRUE
      return(data)
    }, error = function(e) {
      cat("⚠️  Error on attempt", attempt, ":", conditionMessage(e), "\n")
      Sys.sleep(5)  # wait a bit before retrying
    })
    attempt <- attempt + 1
  }
  stop("❌ Failed to retrieve data for ", barcode, " after ", max_retries, " attempts.")
}


#-----------------------------------------------------------------------------------------------------------------
# Function: perform_DEA
#
# Description:
# Performs a simple single-sample differential expression analysis by comparing
# one tumor sample (TCGA unstranded counts) to a mean normal expression profile.
# Computes log2 fold-change (Tumor vs Normal), filters genes by |logFC| >= 1,
# splits into up/down sets, and writes the results to disk. If multiple tumor
# columns are present, only the first is used and the sample can be logged.
#
# Parameters:
#   tumor_data    - SummarizedExperiment returned by GDCprepare (expects assay "unstranded")
#   normal_data   - Numeric matrix of normal mean counts (one column), rownames = genes
#   output_prefix - Sample identifier used for output file naming
#   output_dir    - Output directory path (must exist)
#   log_path      - Optional file path to append sample IDs with multiple tumor columns
#
# Output:
#   Returns a data.frame of filtered DEGs (Gene, logFC) and writes:
#     <output_prefix>_UP_gene_logfc.txt
#     <output_prefix>_DOWN_gene_logfc.txt
#-----------------------------------------------------------------------------------------------------------------

perform_DEA <- function(tumor_data, normal_data, output_prefix, output_dir,log_path = NULL) {
  tumor_counts <- assay(tumor_data, "unstranded")
  normal_counts <- normal_data
  
  # Check if there's more than one tumor column
  if (ncol(tumor_counts) > 1) {
    cat("⚠️  Multiple tumor columns found for", output_prefix, "- selecting the first one only.\n")
    # Log to file if path provided
    if (!is.null(log_path)) {
      write(output_prefix, file = log_path, append = TRUE)
    }
    
    tumor_counts <- tumor_counts[, 1, drop = FALSE]
  }
  
  combined_data <- cbind(tumor_counts, normal_counts)
  
  if (ncol(combined_data) != 2) {
    stop("❌ Unexpected number of columns in combined data: ", ncol(combined_data))
  }
  
  colnames(combined_data) <- c("Tumor", "Normal")
  
  logFC <- log2(combined_data[, "Tumor"] + 1) - log2(combined_data[, "Normal"] + 1)
  degs <- data.frame(Gene = rownames(combined_data), logFC = logFC)
  
  # Apply logFC threshold
  threshold <- 1  # Corresponds to 2-fold change
  degs_filtered <- degs[abs(degs$logFC) >= threshold, ]
  
  upregulated_genes <- degs_filtered[degs_filtered$logFC > 0, ]
  downregulated_genes <- degs_filtered[degs_filtered$logFC < 0, ]
  
  write.table(upregulated_genes, file = paste0(output_dir, output_prefix, "_UP_gene_logfc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(downregulated_genes, file = paste0(output_dir, output_prefix, "_DOWN_gene_logfc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(degs_filtered)
}


#-----------------------------------------------------------------------------------------------------------------
# Function: map_to_hgnc
#
# Description:
# Maps Ensembl gene IDs from a DEG table to HGNC gene symbols using biomaRt.
# Removes Ensembl version suffixes (e.g., ENSG...".12"), merges HGNC symbols
# into the DEG table, drops unmapped genes, splits into up/down by logFC, and
# writes clean HGNC symbol lists to disk.
#
# Parameters:
#   degs          - data.frame of DEGs (expects columns: Gene = Ensembl IDs, logFC)
#   output_prefix - Sample identifier used for output file naming
#   output_dir    - Output directory path (must exist)
#
# Output:
#   Writes two files containing HGNC symbols (one per line):
#     <output_prefix>_UP_gene_list.txt
#     <output_prefix>_DOWN_gene_list.txt
#-----------------------------------------------------------------------------------------------------------------

map_to_hgnc <- function(degs, output_prefix, output_dir) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  degs$Gene <- sub("\\..*", "", degs$Gene)  # Removes version from ENSG IDs
  
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = degs$Gene,
    mart = ensembl
  )
  
  matched_genes <- merge(degs, gene_mapping, by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)
  matched_genes <- matched_genes[!is.na(matched_genes$hgnc_symbol), ]
  
   # Extract and clean gene lists
  upregulated_hgnc <- matched_genes[matched_genes$logFC > 0, "hgnc_symbol"]
  downregulated_hgnc <- matched_genes[matched_genes$logFC < 0, "hgnc_symbol"]

  # Filter out empty or NA entries
  upregulated_hgnc <- upregulated_hgnc[!is.na(upregulated_hgnc) & upregulated_hgnc != ""]
  downregulated_hgnc <- downregulated_hgnc[!is.na(downregulated_hgnc) & downregulated_hgnc != ""]

  # Save clean lists
  write.table(upregulated_hgnc, file = paste0(output_dir, output_prefix, "_UP_gene_list.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(downregulated_hgnc, file = paste0(output_dir, output_prefix, "_DOWN_gene_list.txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}