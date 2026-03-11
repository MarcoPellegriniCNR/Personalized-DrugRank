#----------------------------------------------------------------------------------
# Script: TCGA population-level differential expression analysis (edgeR)
#
# Purpose:
# Performs population-level differential expression analysis between tumor
# and normal samples from a TCGA project using edgeR. The script:
#   1) Downloads and prepares TCGA RNA-seq count data (TCGAbiolinks)
#   2) Filters samples to Primary Tumor (TP) and Normal Tissue (NT)
#   3) Runs edgeR differential expression analysis (Tumor vs Normal)
#   4) Saves full DEG results and filtered up/down gene tables
#   5) Maps Ensembl gene IDs to HGNC symbols using biomaRt
#   6) Writes clean HGNC up/down gene lists for downstream analyses
#
# Configuration:
#   project - TCGA project identifier
#             (e.g., "TCGA-BRCA", "TCGA-COAD", "TCGA-READ", "TCGA-STAD")
#   abb     - Abbreviation used for output directory naming
#             (e.g., "BRCA", "CRC_COAD", "CRC_READ", "GC")
#
# Required libraries:
#   TCGAbiolinks, SummarizedExperiment, edgeR, biomaRt
#
# Outputs:
#   1-input_files/DEGs/<abb>_population/
#       <project>_population_edgeR_results.txt
#       <project>_population_UP_gene_edgeR.txt
#       <project>_population_DOWN_gene_edgeR.txt
#       <project>_population_UP_gene_list.txt
#       <project>_population_DOWN_gene_list.txt
#
# Notes:
#   - Differential expression significance threshold: FDR < 0.01
#   - Ensembl version numbers are removed before gene mapping.
#   - The resulting HGNC gene lists can be used for downstream enrichment
#     or drug-repurposing pipelines.
#----------------------------------------------------------------------------------

#### Load Required Libraries ####
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(edgeR)

### Modify based on analysis
project <- "TCGA-BRCA" # options TCGA-COAD, TCGA READ (for CRC), TCGA-STAD (for GC)
abb <- "BRCA"          # options: CRC_COAD, CRC_READ, GC
output_prefix <- paste0(project, "_population")

#### Define Output Directory ####
output_dir <- paste0("1-input_files/DEGs/", abb, "_population/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#### Function to Query and Prepare Data ####
query_and_prepare <- function(project, data_category, data_type, workflow_type) {
  query <- GDCquery(
    project = project,
    data.category = data_category,
    data.type = data_type,
    workflow.type = workflow_type,
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(query)
  return(GDCprepare(query))
}

#### Function to Perform Differential Expression Analysis Using edgeR ####
perform_edgeR_DEA_population <- function(se, output_prefix, output_dir) {
  sample_types <- colData(se)$shortLetterCode
  keep_samples <- sample_types %in% c("TP", "NT")
  counts <- assay(se)[, keep_samples]
  sample_types <- sample_types[keep_samples]
  group <- factor(ifelse(sample_types == "TP", "Tumor", "Normal"))
  
  # Clean rownames (remove Ensembl version numbers)
  rownames(counts) <- sub("\\..*", "", rownames(counts))
  
  # Create DGEList
  dge <- DGEList(counts = counts, group = group)
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~ group)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  res <- glmQLFTest(fit, coef = 2)
  
  degs <- topTags(res, n = Inf, sort.by = "PValue")$table
  degs$Gene <- rownames(degs)
  
  # Save all results
  write.table(degs, file = paste0(output_dir, output_prefix, "_edgeR_results.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # DEG filtering
  upregulated_genes <- degs[degs$logFC > 0 & degs$FDR < 0.01, ]
  downregulated_genes <- degs[degs$logFC < 0 & degs$FDR < 0.01, ]
  
  # Diagnostic messages
  cat("Number of upregulated genes:", nrow(upregulated_genes), "\n")
  cat("Number of downregulated genes:", nrow(downregulated_genes), "\n")
  
  # Save DEG tables
  write.table(upregulated_genes, file = paste0(output_dir, output_prefix, "_UP_gene_edgeR.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(downregulated_genes, file = paste0(output_dir, output_prefix, "_DOWN_gene_edgeR.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(degs)
}

#### Function to Map Ensembl IDs to HGNC Symbols ####
map_to_hgnc <- function(degs, output_prefix, output_dir) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  degs$Gene <- sub("\\..*", "", degs$Gene)
  
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = degs$Gene,
    mart = ensembl
  )
  
  matched_genes <- merge(degs, gene_mapping, by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)
  matched_genes <- matched_genes[!is.na(matched_genes$hgnc_symbol), ]
  
  upregulated_hgnc <- matched_genes[matched_genes$logFC > 0 & matched_genes$FDR < 0.01, "hgnc_symbol"]
  downregulated_hgnc <- matched_genes[matched_genes$logFC < 0 & matched_genes$FDR < 0.01, "hgnc_symbol"]
  
  # Remove empty or whitespace-only entries
  upregulated_hgnc <- upregulated_hgnc[nzchar(trimws(upregulated_hgnc))]
  downregulated_hgnc <- downregulated_hgnc[nzchar(trimws(downregulated_hgnc))]
  
  write.table(upregulated_hgnc, file = paste0(output_dir, output_prefix, "_UP_gene_list.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(downregulated_hgnc, file = paste0(output_dir, output_prefix, "_DOWN_gene_list.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  cat("HGNC-mapped gene lists saved.\n")
}

#### Main Execution ####
se <- query_and_prepare(project, "Transcriptome Profiling", "Gene Expression Quantification", "STAR - Counts")
degs <- perform_edgeR_DEA_population(se, output_prefix, output_dir)
map_to_hgnc(degs, output_prefix, output_dir)

cat("Population-level analysis complete. Results saved in:\n", output_dir, "\n")
