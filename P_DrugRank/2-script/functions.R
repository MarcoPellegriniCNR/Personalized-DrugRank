########################################################################
#                               Functions
########################################################################

#-----------------------------------------------------------------------
#                              Directories
#-----------------------------------------------------------------------

# Functions that create directories in output = 'DrugRank/3-Drugs'
create_directory_GEO_CRC <- function(output, sample_id){
  d <- c('drugperturbationGEO')
  s <- c('Colorectal_cancer')
  
  for(i in d){
    for(j in s){
      if(dir.exists(paste0(output, sample_id,'/', j,'_',i)) == FALSE){
        dir.create(paste0(output, sample_id,'/',j,'_',i), recursive = TRUE)}
    }
  }
  if(dir.exists(paste0(output, sample_id,'/Drug-modulation-data')) == FALSE){
    dir.create(paste0(output, sample_id,'/Drug-modulation-data'))}
  f <- c('colorectal-cancer')
  for (k in f) {
    if(dir.exists(paste0(output, sample_id,'/Drug-modulation-data/',k,'-cp')) == FALSE){
      dir.create(paste0(output, sample_id,'/Drug-modulation-data/',k,'-cp'))}
  }
}

create_directory_GEO_BRCA <- function(output, sample_id){
  d <- c('drugperturbationGEO')
  s <- c('Breast_cancer')
  
  for(i in d){
    for(j in s){
      if(dir.exists(paste0(output, sample_id,'/', j,'_',i)) == FALSE){
        dir.create(paste0(output, sample_id,'/',j,'_',i), recursive = TRUE)}
    }
  }
  if(dir.exists(paste0(output, sample_id,'/Drug-modulation-data')) == FALSE){
    dir.create(paste0(output, sample_id,'/Drug-modulation-data'))}
  f <- c('breast-cancer')
  for (k in f) {
    if(dir.exists(paste0(output, sample_id,'/Drug-modulation-data/',k,'-cp')) == FALSE){
      dir.create(paste0(output, sample_id,'/Drug-modulation-data/',k,'-cp'))}
 }
}

create_directory_GEO_GC <- function(output, sample_id){
  d <- c('drugperturbationGEO')
  s <- c('Stomach_cancer')
  
  for(i in d){
    for(j in s){
      if(dir.exists(paste0(output, sample_id,'/', j,'_',i)) == FALSE){
        dir.create(paste0(output, sample_id,'/',j,'_',i), recursive = TRUE)}
    }
  }
  if(dir.exists(paste0(output, sample_id,'/Drug-modulation-data')) == FALSE){
    dir.create(paste0(output, sample_id,'/Drug-modulation-data'))}
  f <- c('stomach-cancer')
  for (k in f) {
    if(dir.exists(paste0(output, sample_id,'/Drug-modulation-data/',k,'-cp')) == FALSE){
      dir.create(paste0(output, sample_id,'/Drug-modulation-data/',k,'-cp'))}
  }
}
#-------------------------------------------------------------------------

# Helper to check directory existence
safe_dir_create <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# Helper to get directory name
get_directory_name <- function(d){
  if(d == 'LINCS_L1000_Chem_Pert_up' | d == 'LINCS_L1000_Chem_Pert_down'){
    d1 <- 'L1000'
  }
  if(d == 'Old_CMAP_down' | d == 'Old_CMAP_up'){
    d1 <- 'CMAP'
  }
  if(d == 'DrugMatrix'){
    d1 <- 'drugmatrix'
  }
  if(d == 'Drug_Perturbations_from_GEO_2014'){
    d1 <- 'drugperturbationGEO'
  }
  return(d1)
}

#------------------------------------------------------------------------
#                                 Logging
#------------------------------------------------------------------------
start_log <- function(logfile) {
  safe_dir_create(dirname(logfile))
  con <- file(logfile, open = "at")   # append text connection
  sink(con, split = TRUE)             # stdout -> con
  sink(con, type = "message")         # messages/warnings -> con
  return(con)
}

end_log <- function(con = NULL) {
  # close both sinks
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number() > 0) sink()

  # close connection if provided
  if (!is.null(con)) close(con)
}

#------------------------------------------------------------------------
#                                 Enrichr
#------------------------------------------------------------------------


# Extract sample id for file naming
extract_sample_id <- function(folder_name) {
  sub(".*_(TCGA-.*)", "\\1", folder_name)
}

# Helper to read the DEG lists
read_gene_list <- function(path) {
  readr::read_lines(path, progress = FALSE) |> trimws() |> (\(x) x[nzchar(x)])()
}

#----------------------------------------------------------------------------------
# Function: gene_disease_enrichr
#
# Description:
# Computes, for each Enrichr term, how many genes in the term’s gene list overlap
# with a provided disease gene set. Adds two columns:
#   - perc: percent of term genes overlapping the disease set
#   - num : number of overlapping genes
#
# Parameters:
#   gene_database - Vector of disease-associated genes (symbols/aliases)
#   enrichr       - Enrichr results data.frame containing a 'Genes' column
#                  (semicolon-separated gene lists)
#
# Output:
#   The Enrichr data.frame with added columns 'perc' and 'num' (or NULL on error).
#----------------------------------------------------------------------------------

gene_disease_enrichr <- function(gene_database, enrichr) {
  
  # Ensure enrichr is a data frame and has a Genes column
  if (!is.data.frame(enrichr) || !"Genes" %in% colnames(enrichr)) {
    cat("❌ enrichr is not a valid data frame with a Genes column. Got class:", class(enrichr), "\n")
    return(NULL)
  }
  
  alias <- alias2SymbolTable(gene_database, species = "Hs")
  alias[which(is.na(alias))] <- gene_database[which(is.na(alias))]
  
  cat("🧪 enrichr input class:", class(enrichr), "\n")
  cat("🧪 enrichr column names:", paste(colnames(enrichr), collapse = ", "), "\n")
  
  perc <- c()
  num <- c()
  # Check if enrichr$Genes is a character vector
  for (g in enrichr$Genes) {
    genes <- unlist(strsplit(as.character(g),';'))
    genes <- alias2SymbolTable(genes, species = 'Hs')
    disease_genes <- length(intersect(genes, alias))
    p <- (disease_genes * 100) / length(genes)
    perc <- c(perc, p)
    num <- c(num, disease_genes)
  }
  
  enrichr <- cbind(enrichr, perc, num)
  return(enrichr)
}


#----------------------------------------------------------------------------------
# Function: final_enrichr_sm
#
# Description:
# Post-processes Enrichr results by annotating each term with overlap against:
#   (1) DisGeNET disease genes, (2) NCG cancer genes, and (3) sample DEGs (up∪down).
# Adds percent and count columns for each source.
#
# Parameters:
#   pathway_final - Enrichr results data.frame (must include 'Genes' column)
#   disease       - Disease label ("Stomach_cancer", "Colorectal_cancer", "Breast_cancer")
#   up            - Vector of upregulated DEG symbols
#   down          - Vector of downregulated DEG symbols
#
# Output:
#   pathway_final with added columns:
#   %DisGeNet, num_DisGeNet, %NCG, num_NCG, %DEGs, num_DEGs
#----------------------------------------------------------------------------------

final_enrichr_sm <- function(pathway_final,disease,up,down){
  
  disgenet <- read_tsv('1-input_files/curated_gene_disease_associations.tsv')
  ncg <- read_tsv('1-input_files/NCG6_cancergenes.csv')
  
  if(disease == 'Stomach_cancer'){
    pattern <- 'gastric'
  }
  if(disease == 'Colorectal_cancer'){
    pattern <- 'colorectal'
  }
    if(disease == 'Breast_cancer'){
    pattern <- 'breast'
  }
  
  #include percentage of disease genes
  disgenet <- disgenet[grep(pattern,disgenet$diseaseName,ignore.case = T),]
  pathway_final <- gene_disease_enrichr(disgenet$geneSymbol,pathway_final) 
  colnames(pathway_final)[c(10,11)] <- c('%DisGeNet','num_DisGeNet')
  
  ncg <- ncg[grep(pattern,ncg$cancer_type, ignore.case = T),]
  pathway_final <- gene_disease_enrichr(ncg$symbol,pathway_final)
  colnames(pathway_final)[c(12,13)] <- c('%NCG','num_NCG')
  
  #pathway_final <- gene_disease_enrichr(union(up$V1,down$V1),pathway_final)
  pathway_final <- gene_disease_enrichr(union(up,down),pathway_final)
  colnames(pathway_final)[c(14,15)] <- c('%DEGs','num_DEGs')
  
  return(pathway_final)
  
}


#----------------------------------------------------------------------------------
# Function: enrichr_analysis
#
# Description:
# Performs pathway enrichment analysis using Enrichr on genes extracted from
# active subnetwork results. Gene lists are derived according to the specified
# module detection method (e.g., Core&Peel or MD), converted to official gene
# symbols, and queried against one or more Enrichr databases. Significant
# results (Adjusted P-value ≤ 0.05) are optionally post-processed and saved
# to disk.
#
# Parameters:
#   input      - Path to directory containing the active subnetwork file
#   file       - Name of the active subnetwork file
#   method     - Module detection method used (e.g., "Core&Peel", "MD")
#   disease    - Disease name (used for output organization and filtering logic)
#   database   - Character vector of Enrichr database names to query
#   background - Background gene list (reserved for compatibility; not used here)
#   up         - Optional vector of upregulated genes (used for post-processing)
#   down       - Optional vector of downregulated genes (used for post-processing)
#   output     - Base output directory for saving enrichment results
#
# Details:
#   - Automatically skips databases if output already exists.
#   - Includes retry logic to handle Enrichr API rate limits (HTTP 429).
#   - Filters results by Adjusted P-value ≤ 0.05.
#   - Applies disease-specific post-processing when applicable.
#
# Output:
#   Saves filtered enrichment results as CSV files in structured
#   disease-specific directories.
#----------------------------------------------------------------------------------

enrichr_analysis <- function(input, file, method, disease, database, background, up = NULL, down = NULL, output) {
  print(output)

  # Read the active subnetwork file
  singleModule <- read.table(paste0(input, file))

  # Extract genes based on method
  if (length(grep('Core&Peel', method)) == 1) {
    genes <- singleModule$V1
  }

  if (method == 'MD') {
    genes <- union(singleModule$V1, singleModule$V2)
  }

  # Convert gene names to official symbols
  alias <- alias2SymbolTable(genes, species = "Hs")
  alias[which(is.na(alias))] <- genes[which(is.na(alias))]

  # Perform enrichment analysis for each database
  for (d in database) {
    d1 <- get_directory_name(d)
    output_directory <- paste0(output, disease, '_', d1)
    suffix <- "_perc.csv"
    out_file <- paste0(output_directory, '/enrichR_', method, '_', d, suffix)

    # Skip if output file already exists
    if (file.exists(out_file)) {
      cat("⏭️ Skipping", d, "- already exists for", method, "\n")
      next
    }

    cat("🔍 Querying", d, "for", method, "\n")

    # Retry-safe Enrichr call
    e <- NULL
    max_tries <- 5
    wait_time <- 10

    for (attempt in 1:max_tries) {
      e <- tryCatch({
        res <- enrichr(genes = alias, databases = d)
        if (is.null(res) || is.null(res[[1]])) {
          cat("⚠️ No results returned for", d, "- skipping\n")
          NULL
        } else {
          res[[1]]
        }
      }, error = function(err) {
        msg <- conditionMessage(err)
        if (grepl("Too Many Requests|429", msg)) {
          cat(sprintf("⚠️ Rate limit hit for %s. Waiting %d seconds (Attempt %d of %d)...\n",
                      d, wait_time, attempt, max_tries))
          Sys.sleep(wait_time)
          return(NULL)
        } else {
          cat("❌ Enrichr error for", d, ":", msg, "\n")
          return(NULL)
        }
      })

      if (!is.null(e)) break
    }

    # Filter and save results
    if (!is.null(e) && nrow(e) > 0 && "Adjusted.P.value" %in% names(e)) {
      e <- e[e$Adjusted.P.value <= 0.05, , drop = FALSE]

      if (nrow(e) > 0) {
        if (disease %in% c("Stomach_cancer", "Colorectal_cancer", "Breast_cancer")) {
          e <- final_enrichr_sm(e, disease, up, down)
        }

        dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
        write.csv(e, out_file, row.names = FALSE)
        cat("✅ Saved enrichment results to", out_file, "\n")
      } else {
        cat("ℹ️ No significant results for", d, "in", method, "\n")
      }
    }
  }
}


#----------------------------------------------------------------------------------
# Function: enrichr_analysis_multipleModules
#
# Description:
# Performs Enrichr enrichment analysis separately for each module detected
# in a multiple-module output file. Modules are parsed according to the
# specified method (e.g., Core&Peel or tab-delimited formats), gene symbols
# are standardized, and each module is queried independently against one or
# more Enrichr databases. Significant results are filtered and saved in a
# structured directory hierarchy.
#
# Parameters:
#   file       - Path to the multiple-module file containing gene sets
#   disease    - Disease name (used for output directory structure)
#   method     - Module detection method used (e.g., "Core&Peel", "MD")
#   database   - Character vector of Enrichr database names to query
#   sample_id  - Sample identifier used for output organization and logging
#
# Details:
#   - Automatically detects file format:
#       * Core&Peel → space-delimited modules (one module per line)
#       * Other methods → tab-delimited format
#   - Converts gene aliases to official human gene symbols.
#   - Runs enrichment independently for each module.
#   - Skips analysis if output file already exists.
#   - Filters results using raw P.value ≤ 0.05.
#   - Logs Enrichr errors to 'logs/enrichr_multipleModules_errors.csv'
#     including sample, method, module index, database, and timestamp.
#
# Output:
#   Saves enrichment results for each module and database as CSV files
#   under:
#   4-MultipleModules_output/<disease>/<sample_id>/<method>/<database>/
#----------------------------------------------------------------------------------


enrichr_analysis_multipleModules <- function(file, disease, method, database, sample_id) {
  
  # Internal helper for output directory names
  get_directory_name <- function(d) {
    if (d %in% c("LINCS_L1000_Chem_Pert_up", "LINCS_L1000_Chem_Pert_down")) return("L1000")
    if (d %in% c("Old_CMAP_up", "Old_CMAP_down")) return("CMAP")
    return(d)
  }

  # Detect delimiter based on method
  if (grepl("Core&Peel", method)) {
    sep <- " "
    cat("📄 Detected Core&Peel format. Using space as delimiter.\n")
    modules_raw <- scan(file, what = "", sep = "\n", quiet = TRUE)
    modules_list <- lapply(modules_raw, function(line) {
      g <- unlist(strsplit(line, split = sep))
      g[g != ""]
    })
  } else {
    sep <- "\t"
    cat("📄 Assuming tab-delimited format.\n")
    modules_df <- read.table(file, header = FALSE, sep = sep, stringsAsFactors = FALSE, fill = TRUE)
    modules_list <- apply(modules_df, 1, function(row) {
      g <- as.character(na.omit(row))
      g[g != ""]
    })
  }

  cat("------------------------------------------------------------\n")
  cat("🧬 Number of modules found:", length(modules_list), "\n")

  if (!is.vector(database)) {
    database <- list(database)
  }

  for (i in seq_along(modules_list)) {
    module_genes <- modules_list[[i]]

    # Gene alias conversion
    alias <- alias2SymbolTable(module_genes, species = "Hs")
    alias[is.na(alias)] <- module_genes[is.na(alias)]

    for (d in database) {
      d1 <- get_directory_name(d)
      output_dir <- file.path("4-MultipleModules_output", disease, sample_id, method, d1)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

      output_file <- file.path(output_dir, paste0("enrichR_module", i, "_", d, ".csv"))

      if (file.exists(output_file)) {
        cat("⏭️ Skipping module", i, "-", d, "- already exists\n")
        next
      }

      cat("🔍 Querying Enrichr for module", i, "-", d, "\n")

      e <- tryCatch({
        res <- enrichr(genes = alias, databases = d)
        if (is.null(res) || is.null(res[[1]])) {
          cat("⚠️ No results returned for module", i, "-", d, "\n")
          return(NULL)
        }
        res[[1]]
      }, error = function(err) {
        cat("❌ Enrichr error for module", i, "-", d, ":", conditionMessage(err), "\n")
        log_entry <- data.frame(
          sample_id = sample_id,
          method = method,
          module = i,
          database = d,
          error = gsub("[\r\n]", " ", conditionMessage(err)),
          time = as.character(Sys.time()),
          stringsAsFactors = FALSE
        )
        log_file <- paste0("logs/", disease, "/enrichr_multipleModules_errors.csv")
        write.csv(
          log_entry,
          file = log_file,
          append = file.exists(log_file),
          row.names = FALSE,
          col.names = !file.exists(log_file)
        )
        return(NULL)
      })

      # Filter by raw P.value ≤ 0.05
      if (!is.null(e) && nrow(e) > 0 && "P.value" %in% names(e)) {
        e <- e[e$P.value <= 0.05, , drop = FALSE]

        if (nrow(e) > 0) {
          write.csv(e, output_file, row.names = FALSE)
          cat("✅ Saved enrichment for module", i, "-", d, "to", output_file, "\n")
        } else {
          cat("ℹ️ No enrichment below P.value 0.05 for module", i, "-", d, "\n")
        }
      }
    }
  }
}



#----------------------------------------------------------------------------------
# Function: add_pvalue
#
# Description:
# Appends Enrichr P.value and Adjusted.P.value to a drug-by-module results table.
# For each module, the function reads the corresponding Enrichr output file,
# matches table$drugs to Enrichr$Term, and extracts statistical values.
# Missing files or unmatched terms are filled with NA.
#
# Parameters:
#   path_input  - Directory containing per-module Enrichr results
#   method      - Module detection method ("Core&Peel" or "MD")
#   d           - Database label (e.g., "GEO", "L1000", "CMAP", etc.)
#   table       - Data frame with at least 'modules' and 'drugs' columns
#   path_output - Directory for saving the updated table
#   disease     - Disease name used in output filename
#   param       - Optional Core&Peel parameter string
#
# Output:
#   Returns the updated table (with p_value and adjusted_pvalue columns)
#   and writes it to disk as a CSV file.
#----------------------------------------------------------------------------------

add_pvalue <- function(path_input,method,d,table,path_output,disease,param = NULL){
  
  if(d == 'L1000'){
    d1 <- 'L1000'
  }
  if(d == 'CMAP'){
    d1 <- 'CMAP'
  }
  if(d == 'drugmatrix'){
    d1 <- 'DrugMatrix'
  }
  if(d == 'drugperturbationGEO'){
    d1 <- 'Drug_Perturbations_from_GEO_2014'
  }
  if(method == 'Core&Peel'){
    m <- paste0(method,'-',strsplit(param,split = '-')[[1]][2],'-',strsplit(param,split = '-')[[1]][4])
  }
  if(method == 'MD'){
    m <- 'MD'
  }
  if(d == 'L1000' | d == 'CMAP'){
    path1 <- paste0(path_input,'/',m,'/',d1,'/','merge_enrich')
  }
  else{
    path1 <- paste0(path_input,'/',m,'/',d1)
  }
  p_value <- c()
  adjusted_pvalue <- c()
  
  for(j in unique(table[,'modules'])){
    if (d == "L1000") {
    file_path <- paste0(path1, '/enrichR_module', j, '_LINCS_L1000_Chem_Pert.csv')
    } else if (d == "CMAP") {
    file_path <- paste0(path1, '/enrichR_module', j, '_CMAP.csv')
    } else {
    file_path <- paste0(path1, '/enrichR_module', j, '_', d1, '.csv')
    }
    print(paste("Checking file:", file_path))  # Debugging
    
    if (file.exists(file_path)) {
        enrich <- read.csv(file_path)
        subtable <- as.matrix(table[which(table[,'modules'] == j),]) 
  
      if(class(subtable)[1] == 'matrix' & ncol(subtable) > 1){
        p_value <- c(p_value, enrich$P.value[match(subtable[,'drugs'],enrich$Term)])
        adjusted_pvalue <- c(adjusted_pvalue, enrich$Adjusted.P.value[match(subtable[,'drugs'],enrich$Term)])
        }
      
      else{
        #print(j)
        p_value <- c(p_value, enrich$P.value[match(subtable['drugs'],enrich$Term)])
        adjusted_pvalue <- c(adjusted_pvalue, enrich$Adjusted.P.value[match(subtable['drugs'],enrich$Term)])
        }
      
    }
    else{
      print(paste("File not found:", file_path))  # Debugging
      subtable <- table[which(table[,'modules'] == j),]
      if(class(subtable)[1] == 'character'){
        p_value <- c(p_value,rep(NA,1))
        adjusted_pvalue <- c(adjusted_pvalue,rep(NA,1))
        
      }
      else{
        p_value <- c(p_value,rep(NA,nrow(subtable)))
        adjusted_pvalue <- c(adjusted_pvalue,rep(NA,nrow(subtable)))
      }
    }
  }
  
  table <- cbind(table,p_value,adjusted_pvalue)
  
  if(d == 'drugperturbationGEO'){
    d1 <- 'GEO'
  }
  if(method == 'Core&Peel'){
    m2 <- paste0(method,'-coexpr025marta-r1-nl10',param)
  }
  if(method == 'MD'){
    m2 <- 'MD'
  }
  write.csv(table, file = paste0(path_output,'/',m2,'_',disease,'_',d1,'_noFiltering.csv'), row.names = F)
  return(table)
}

#----------------------------------------------------------------------------------
# Function: read_enrichr
#
# Description:
# Reads Enrichr CSV results from a disease/database folder and returns a single
# data frame. Files are filtered by method (MD or Core&Peel + param). For L1000
# and CMAP, up/down files are loaded and merged, with "-up"/"-dn" appended to Term.
#
# Parameters:
#   path_enrich - Base path containing "<disease>_<d>/" directories
#   method      - Method name ("MD" or "Core&Peel") used to filter filenames
#   disease     - Disease label used in directory naming
#   d           - Database label (e.g., "GEO", "L1000", "CMAP", "drugmatrix", "drugperturbationGEO")
#   param       - Optional Core&Peel parameter string used to filter filenames
#
# Output:
#   A data.frame of Enrichr results (combined for up/down where applicable).
#----------------------------------------------------------------------------------

read_enrichr <- function(path_enrich, method, disease, d, param = NULL) {
  
  message("Reading Enrichr files from: ", path_enrich, disease, "_", d)
  files <- list.files(paste0(path_enrich, disease, "_", d), pattern = '.csv')
  message(length(files), " CSV file(s) found.")

  # Filter based on method
  if (method == 'MD') {
    files <- files[grep('enrichR_MD_', files)]
    message("Filtering for method 'MD': ", length(files), " file(s) remaining.")
    message("Files: ", paste(files, collapse = ", "))
  }

  if (length(grep('Core&Peel', method)) >= 1) {
    files <- files[grep(paste0('enrichR_', method), files)]
    message("After filtering for method '", method, "': ", length(files), " file(s).")
    
    files <- files[grep(param, files)]
    message("After filtering with param '", param, "': ", length(files), " file(s).")
    message("Files: ", paste(files, collapse = ", "))
  }

  # Load files based on database type
  path_full <- paste0(path_enrich, disease, "_", d, "/")
  if (d == 'L1000') {
    message("Loading LINCS L1000 up/down regulated files...")
    f_d <- files[grep('LINCS_L1000_Chem_Pert_down', files)]
    f_u <- files[grep('LINCS_L1000_Chem_Pert_up', files)]
    
    enrich_down <- read.csv(paste0(path_full, f_d))
    enrich_down$Term <- paste(enrich_down$Term, '-dn', sep = '')
    
    enrich_up <- read.csv(paste0(path_full, f_u))
    enrich_up$Term <- paste(enrich_up$Term, '-up', sep = '')
    
    enrichr <- as.data.frame(rbind(enrich_down, enrich_up))
    message("LINCS L1000 data combined: ", nrow(enrichr), " total terms.")
  }

  if (d == 'CMAP') {
    message("Loading Old CMAP up/down regulated files...")
    f_d <- files[grep('Old_CMAP_down', files)]
    f_u <- files[grep('Old_CMAP_up', files)]

    enrich_down <- read.csv(paste0(path_full, f_d))
    enrich_down$Term <- paste(enrich_down$Term, '-dn', sep = '')

    enrich_up <- read.csv(paste0(path_full, f_u))
    enrich_up$Term <- paste(enrich_up$Term, '-up', sep = '')

    enrichr <- as.data.frame(rbind(enrich_down, enrich_up))
    message("Old CMAP data combined: ", nrow(enrichr), " total terms.")
  }

  if (d == 'drugmatrix' | d == 'drugperturbationGEO'| d== "Drug_Perturbations_from_GEO_2014") {
    message("Loading data from: ", files)
    enrichr <- read.csv(paste0(path_full, files))
    message("Loaded ", nrow(enrichr), " rows from ", files)
  }

  message("Enrichment reading complete.\n")
  return(enrichr)
}


#----------------------------------------------------------------------------------
# Function: get_modulation_table
#
# Description:
# Builds a drug–module modulation table by intersecting module gene sets with
# Enrichr drug signatures. For each module and each drug term, it computes:
# module size, #DEGs (up/down) in the module, #genes shared with the drug
# signature, and #shared up/down DEGs. Rows with zero overlap are removed.
#
# Parameters:
#   enrichr - Enrichr results data.frame (must include 'Term' and 'Genes' columns)
#   method  - Module method ("MD" or "Core&Peel"); used to select and parse module file
#   path    - Directory containing .txt module file(s)
#   up      - Vector of upregulated DEG symbols
#   down    - Vector of downregulated DEG symbols
#
# Output:
#   Matrix/data.frame with columns:
#   modules, drugs, module_size, DEG_up, DEG_down, drug_module_genes,
#   DE_up_drug, DE_down_drug
#----------------------------------------------------------------------------------

get_modulation_table <- function(enrichr,method,path,up,down){
  
  # list full paths, not just basenames
  files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No .txt module files found in: ", path)
  }

  # choose the right file
  if (method == "MD") {
    hits <- grep("MD", basename(files), value = TRUE)
    file_use <- files[match(hits, basename(files))][1]
  } else {
    hits <- grep("Core&Peel", basename(files), value = TRUE)
    file_use <- files[match(hits, basename(files))][1]
  }

  if (is.na(file_use) || length(file_use) == 0) {
    stop("No module file matched method='", method, "' in: ", path,
         "\nFiles seen: ", paste(basename(files), collapse = ", "))
  }

  if (!file.exists(file_use)) {
    stop("Module file does not exist (constructed path wrong?): ", file_use)
  }

  # read modules
  if (method == "MD") {
    modules <- scan(file_use, what = "", sep = "\n", quiet = TRUE)
    modules1 <- strsplit(modules, split = "\t")
  } else {
    modules <- scan(file_use, what = "", sep = "\n", quiet = TRUE)
    modules1 <- strsplit(modules, split = " ")
  }
  enrichr1 <- strsplit(as.character(enrichr$Genes), split = ';')
  names(enrichr1) <- enrichr$Term
  
  l <- lapply(modules1, function(x) lapply(enrichr1, function(y) intersect(x,y)))
  
  size <- lapply(modules1, function(x) length(x))
  size_module_gene <- lapply(l,function(x) lengths(x))
  
  DEGs_up <- lapply(modules1, function(x) length(intersect(x,up)))
  DEGs_down <- lapply(modules1, function(x) length(intersect(x,down)))
  
  DEG_up_drug <- vector(mode = 'list')
  for(j in seq(1,length(modules))){
    DEG_up_drug[[j]] <- lapply(l[[j]], function(x) length(intersect(unlist(x),up)))
  }
  
  DEG_down_drug <- vector(mode = 'list')
  for(j in seq(1,length(modules))){
    DEG_down_drug[[j]] <- lapply(l[[j]], function(x) length(intersect(unlist(x),down)))
  }
  
  active_modules <- c()
  drugs <- c()
  module_size <- c()
  drug_module_genes <- c()
  DEG_up <- c()
  DEG_down <- c()
  DE_up_drug <- c()
  DE_down_drug <- c()
  for(i in seq(1,length(modules))){
    active_modules <- c(active_modules,rep(i,length(l[[i]])))
    drugs <- c(drugs,names(l[[i]]))
    module_size <- c(module_size,rep(size[[i]],length(l[[i]])))
    DEG_up <- c(DEG_up,rep(DEGs_up[[i]],length(l[[i]])))
    DEG_down <- c(DEG_down,rep(DEGs_down[[i]],length(l[[i]])))
    drug_module_genes <- c(drug_module_genes,size_module_gene[[i]])
    DE_up_drug <- c(DE_up_drug,DEG_up_drug[[i]])
    DE_down_drug <- c(DE_down_drug,DEG_down_drug[[i]])
    
  }
  
  DE_up_drug <- unlist(DE_up_drug)
  names(DE_up_drug) <- NULL
  
  DE_down_drug <- unlist(DE_down_drug)
  names(DE_down_drug) <- NULL
  
  table <- cbind(active_modules,drugs,module_size,DEG_up,DEG_down,drug_module_genes,DE_up_drug,DE_down_drug)
  table <- table[which(table[,'drug_module_genes'] != 0),]
  colnames(table) <- c('modules','drugs','module_size','DEG_up','DEG_down','drug_module_genes','DE_up_drug','DE_down_drug')
  return(table)
}

#----------------------------------------------------------------------------------
# Function: get_modulation_filter
#
# Description:
# Filters a drug–module modulation table to keep only “direction-consistent”
# pairs: modules with only up-DEGs are matched to down-drug signatures, and
# modules with only down-DEGs are matched to up-drug signatures. Computes a
# single DE_drug column (max of DE_up_drug / DE_down_drug), drops the originals,
# reorders columns, and writes the filtered table to CSV.
#
# Parameters:
#   table1      - Modulation table (expects columns: DEG_up, DEG_down, drugs,
#                 drug_module_genes, DE_up_drug, DE_down_drug, etc.)
#   d           - Database label (e.g., "GEO", "L1000", "CMAP", "drugmatrix", "drugperturbationGEO")
#   method      - Method label ("Core&Peel" or "MD") used in output naming
#   path_output - Output directory for the filtered CSV
#   disease     - Disease name used in output filename
#   param       - Optional Core&Peel parameter string used in output naming
#
# Output:
#   Writes the filtered table to:
#   <path_output>/<method>_<disease>_<db>.csv
#----------------------------------------------------------------------------------

get_modulation_filter <- function(table1,d,method,path_output,disease,param = NULL){
  
  index_down <- which(table1[,'DEG_up'] > 0 & table1[,'DEG_down'] == 0)
  index_up <- which(table1[,'DEG_up'] == 0 & table1[,'DEG_down'] > 0)
  if(d == 'L1000' | d == 'CMAP' | d == 'drugmatrix'){
    index_drug_down <- grep('-dn',table1[,'drugs'])
    index_drug_up <- grep('-up',table1[,'drugs'])
  }
  if(d == 'drugperturbationGEO'){
    index_drug_down <- grep(' down',table1[,'drugs'])
    index_drug_up <- grep(' up',table1[,'drugs'])
  }
  
  index_down_final <- intersect(index_down,index_drug_down)
  index_up_final <- intersect(index_up,index_drug_up)
  
  index <- which(table1[,'drug_module_genes'] > 0)
  index_final <- union(index_down_final,index_up_final)
  table2 <- table1[intersect(index,index_final),]
  
  table2[,'DE_drug'] <- apply(cbind(table2[,'DE_up_drug'],table2[,'DE_down_drug']),1,function(x) ifelse(x[1] > 0,x[1], x[2]))
  
  table2[,'DE_up_drug'] <- NULL
  table2[,'DE_down_drug'] <- NULL
  table2 <- table2[,c(1,2,3,4,5,6,9,7,8)]
  
  if(d == 'L1000'){
    d1 <- 'L1000'
  }
  if(d == 'CMAP'){
    d1 <- 'CMAP'
  }
  if(d == 'drugmatrix'){
    d1 <- 'DrugMatrix'
  }
  if(d == 'drugperturbationGEO'){
    d1 <- 'GEO'
  }
  if(d == 'Drug_Perturbations_from_GEO_2014'){
    d1 <- 'GEO'
  }
  if(method == 'Core&Peel'){
    m <- paste0(method,'-coexpr025marta-r1-nl10',param)
  }
  if(method == 'MD'){
    m <- 'MD'
  }
  write.csv(table2, file = paste0(path_output,'/',m,'_',disease,'_',d1,'.csv'), row.names = F)
}