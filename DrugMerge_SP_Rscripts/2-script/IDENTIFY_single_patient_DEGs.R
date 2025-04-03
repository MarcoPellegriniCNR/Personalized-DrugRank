#### Load Required Libraries ####
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)


#### Define Output Directory ####
output_dir <- "D:/Progetti/CNR/DrugMerge_single_patient/1-input_files/DEGs/Colorectal_cancer_sp/"  # Specify the folder where output files will be saved
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#### Mean Gene Expression in Normal Samples (n=41) ####
mean_expression_normal_samples <- read.csv('D:/Progetti/CNR/DrugMerge_single_patient/1-input_files/Normal_control_crc/normal_ensid_count.csv', row.names = 1)
normal_data <- data.matrix(mean_expression_normal_samples)

#### Function to Query and Prepare Data ####
query_and_prepare_tumor <- function(project, data_category, data_type, workflow_type, barcode) {
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

#### Function to Perform Differential Expression Analysis ####
perform_DEA <- function(tumor_data, normal_data, output_prefix) {
  tumor_counts <- assay(tumor_data)
  normal_counts <- normal_data
  combined_data <- cbind(tumor_counts, normal_counts)
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
map_to_hgnc <- function(degs, output_prefix) {
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

#### Main Execution ####
project <- "TCGA-COAD"  # TCGA-COAD for colorectal, TCGA-PRAD for prostate
case_id <- "TCGA-A6-5667"
tumor_sample_id <- "TCGA-A6-5667-01A"

tumor_data <- query_and_prepare_tumor(project, "Transcriptome Profiling", "Gene Expression Quantification", "STAR - Counts", tumor_sample_id)

degs <- perform_DEA(tumor_data, normal_data, paste0(project,'_', case_id, '_test'))
map_to_hgnc(degs, paste0(project,'_', case_id, '_test'))

cat("Analysis complete. Results saved in", output_dir,"\n")
