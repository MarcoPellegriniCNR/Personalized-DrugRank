#### Function to Query and Prepare Data ####
query_and_prepare_tumor_simple <- function(project, data_category, data_type, workflow_type, barcode) {
  query <- GDCquery(
    project = project,
    data.category = data_category,
    data.type = data_type,
    workflow.type = workflow_type,
    barcode =  barcode
  )
  GDCdownload(query)
  return(GDCprepare(query))
}

## version that deals with GDC server issue (500 error) ##
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


#### Function to Perform Differential Expression Analysis ####
perform_DEA_debug <- function(tumor_data, normal_data, output_prefix, output_dir) {
  tumor_counts <- assay(tumor_data)
  normal_counts <- normal_data
  combined_data <- cbind(tumor_counts, normal_counts)
  #colnames(combined_data) <- c("Tumor", "Normal")
  cat("Tumor shape:", dim(tumor_counts), "\n") #debug step
  cat("Normal shape:", dim(normal_counts), "\n") #debug step
  cat("Combined shape:", dim(combined_data), "\n") #debug step
  if (ncol(combined_data) == 2) {
    colnames(combined_data) <- c("Tumor", "Normal")
  } else {
    stop("Unexpected number of columns in combined data: ", ncol(combined_data))
  }
  
  logFC <- log2(combined_data[, "Tumor"] + 1) - log2(combined_data[, "Normal"] + 1)
  degs <- data.frame(Gene = rownames(combined_data), logFC = logFC)
  
  upregulated_genes <- degs[degs$logFC > 0, ]
  downregulated_genes <- degs[degs$logFC < 0, ]
  
  write.table(upregulated_genes, file = paste0(output_dir, output_prefix, "_UP_gene_logfc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(downregulated_genes, file = paste0(output_dir, output_prefix, "_DOWN_gene_logfc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(degs)
}


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
  
  upregulated_genes <- degs[degs$logFC > 0, ]
  downregulated_genes <- degs[degs$logFC < 0, ]
  
  write.table(upregulated_genes, file = paste0(output_dir, output_prefix, "_UP_gene_logfc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(downregulated_genes, file = paste0(output_dir, output_prefix, "_DOWN_gene_logfc.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(degs)
}


#### Function to Map Ensembl IDs to HGNC Symbols and generate gene lists ####
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
  
  upregulated_hgnc <- matched_genes[matched_genes$logFC > 0, "hgnc_symbol"]
  downregulated_hgnc <- matched_genes[matched_genes$logFC < 0, "hgnc_symbol"]
  
  write.table(upregulated_hgnc, file = paste0(output_dir, output_prefix, "_UP_gene_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(downregulated_hgnc, file = paste0(output_dir, output_prefix, "_DOWN_gene_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}