#!/usr/bin/env Rscript

######################
# Code to download functionl annotations on patient specific disease modules for stomach adenocarcicoma (STAD) aptients and Genoe othology BP annotaitions

#####################
# command line support

# Get command-line arguments.
# `trailingOnly = TRUE` ensures you only get the arguments *after* the script name
# and any R-specific options (like --args, though Rscript handles this implicitly).
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments were provided
if (length(args) != 3) {
  stop("Error: Wrong parameters provided. Usage: ./cl_My_script_stomach_pipeline_2_functional.R <param1_value: patient id> <param2_value: mode up or down > <param3_value: delay in seconds between network calls", call. = FALSE)
}

# Assign parameters to variables
# Arguments are always read as character strings, so convert if needed
case_id <- args[1]
mode <- args[2]
delay <- as.numeric(args[3])

#param2 <- as.numeric(args[2]) # Example: converting the second parameter to a number

# Now, use your parameters in your R logic
print(paste("First parameter case_id (string):", case_id))
print(paste("Second parameter mode (string):", mode))
print(paste("Third paameter dealy (integer):", delay))





#####################################
# Secod pipeline for handling moduular data


###### set your working directory ######
#setwd('/home/fiamma/DrugMerge/')
#setwd('C:/Users/Marco/Desktop/Personalized_drug_ranking_project/New_R_code_2025_04_15')
#setwd('/Users/pellegrini/Desktop/Personalized_drug_ranking_project/New_R_code_2025_04_15')

library(limma)
library(enrichR)
library(readr)

#source('DrugMerge_SP_Final/2-script/functions.R')
#source('My_functions.R')
#source('My_functions_breast.R')
#source('My_functions_stomach.R')


#---------------------------------------------------------------------------------
#                              Support functions
#----------------------------------------------------------------------------------

#################
# a function to cehck existence or gernate the directory structure

check_gen_dir <- function(rootname, disease, case_id, method, d, mode){
	
	d1 <-  paste0(d , '_',  mode)
  
  #create output directory structure 
  top_dir_name = paste0(rootname, '/' , disease,'/', case_id)
  if(dir.exists(top_dir_name) == FALSE){
    creation_successful = dir.create(top_dir_name, recursive = TRUE)
    if (creation_successful) {
      print(paste("Succesful creation of ", top_dir_name ))
     }
    else 
     {
      print(paste("UnSuccesful creation of ", top_dir_name ))
      stop()
     }
  }

  if(dir.exists(paste0(top_dir_name, '/',method)) == FALSE){
    dir.create(paste0(top_dir_name,'/',method))
  }
  if(dir.exists(paste0(top_dir_name, '/',method,'/',d1)) == FALSE){
    dir.create(paste0(top_dir_name,'/',method,'/',d1))
  }
  dir_output <- paste0(top_dir_name,'/',method,'/',d1)
  return(dir_output)
}  


###########################

enrichR_single_multipleModule <- function(rootname, genes,d,i,disease,case_id,method, mode, delay){
	
	
	# create a composite name for file naming 
	#d1 <-  paste0(d , mode)  
	dir_output = check_gen_dir(rootname, disease, case_id, method, d, mode)
  	filename = paste0(dir_output,'/enrichR_module',i, '.csv')
  
  if (!file.exists(filename)){
	  
	  alias <- alias2SymbolTable(genes, species = "Hs")
	  #if alias == NA then use the gene names in genes
	  alias[which(is.na(alias))] <- genes[which(is.na(alias))]
	  
	  Sys.sleep(delay) # introduces a delay not to oveload the server
	  e <- enrichr(genes = alias,
	               databases = d)
	  e <- e[[1]][which(e[[1]]$P.value <= 0.05),]
	  if (is.null(e) == FALSE) {
	    if (nrow(e)!=0) {
	#      write.csv(e,paste0(dir_output,'/enrichR_module',i,'_',d,'.csv'), row.names = F)
	      write.csv(e,paste0(dir_output,'/enrichR_module',i, '.csv'), row.names = F) # shorter filename
	
	      
	    }
	  }
    }	  
}


#############################


enrichr_analysis_multipleModules <- function(rootname, f,disease,method,d, up, down, mode, delay){
  
  # read the module file
  modules <- scan(f, what = '', sep = '\n')
  print(paste('Number of modules is : ', length(modules)))
  # set the seprator character
  if(length(grep('Core&Peel',method) >=1 )){
    s <- ' '
  }
  if(method == 'MD'){
    s <- '\t'
  }

  #perform interesections with up/down or both or none

  if (mode == 'none') { 
    for(i in seq(1,length(modules))){
      genes <- unlist(strsplit(modules[i], split = s))
      if (length(genes) > 1) {
          enrichR_single_multipleModule(rootname, genes,d,i,disease,case_id,method, mode, delay)
          }
    }
  } 

  if (mode == 'up') { 
    new_up <- alias2SymbolTable(up$V1, species = "Hs")
    #if alias == NA then use the gene names in genes
    new_up[which(is.na(new_up))] <- up$V1[which(is.na(new_up))]
    print('new_up')
    print(length(new_up))
    for(i in seq(1,length(modules))){
      genes <- unlist(strsplit(modules[i], split = s))
      new_genes = intersect(genes, new_up)
      print(c(i, length(genes),length(new_genes)))
      if (length(new_genes) > 1) {
        enrichR_single_multipleModule(rootname, new_genes,d,i,disease,case_id,method, mode, delay)
      }
    }
  }

  if (mode == 'down') { 
    new_down <- alias2SymbolTable(down$V1, species = "Hs")
    #if alias == NA then use the gene names in genes
    new_down[which(is.na(new_down))] <- down$V1[which(is.na(new_down))]
    print('new_down')
    print(length(new_down))
    for(i in seq(1,length(modules))){
      genes <- unlist(strsplit(modules[i], split = s))
      new_genes = intersect(genes, new_down)
      print(c(i, length(genes),length(new_genes)))
      if (length(new_genes) > 1) {
        enrichR_single_multipleModule(rootname, new_genes,d,i,disease,case_id,method, mode, delay)
      }
    }
  }



}



#---------------------------------------------------------------------------------
#                               MAIN
#----------------------------------------------------------------------------------


#background <- read.table('DrugMerge_SP_Final/1-input_files/symbol_entrezID_updated_table.txt', header = T)
background <- read.table('DrugMerge_SP_Final/1-input_files/symbol_entrezID_updated_table.txt', header = T)

###SPECIFY Drug databases###

#database <- c('Drug_Perturbations_from_GEO_2014')
#databases <- c('GO_Biological_Process_2025', 'Reactome_Pathways_2024')
#databases <- c('GO_Molecular_Function_2025')
#databases <- c('KEGG_2021_Human')


alg_name = 'Core&Peel-d0.7-j0.8'
cancer_name = 'Stomach_cancer'
cancer_ID = 'STAD'

# directorty root check/creation
root <- 'DrugMerge_SP_Final/5-MultipleModules_Functional_Enrichment'
if(dir.exists(root) == FALSE){
  dir.create(root)
}



# an old directorys tructure I am not sunign anymore 
#output <- paste0('DrugMerge_SP_Final/3-Drugs/') #, case_id, '_colorectal_cancer/')
#output <- paste0('DrugMerge_SP_Final/3-Drugs/') #, case_id, '_colorectal_cancer/')
#create_directory_breast(output, case_id)
#create_directory_stomach(output, case_id)




####-------------------------------------------------------------------------------
#                           DEG up and down forStomach cancer
#----------------------------------------------------------------------------------



# balanced at 4000 genes max for Stomach STAD  fluorouracil patients
up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Stomach_cancer_sp2_fluorouracil/', case_id, '_4000_UP_gene_list.txt'))
down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Stomach_cancer_sp2_fluorouracil/', case_id, '_4000_DOWN_gene_list.txt'))


print('Structure of up')
head(up)
tail(up)
str(up)
dim(up)

print('Structure of down')
head(down)
tail(down)
str(down)
dim(down)

 

print('Starting Modular EnrichR analysis')


# 1) perform the enrich analysis on the single modules
# Core&Peel
print('Starting multiple-module enrichment analysis for Core&Peel')
print('starting enrichr_analysis_multipleModules')

#f <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'
#f <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Breast_cancer/MultipleModules/BRCA_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'

disease_modules_filename  <-  'DrugMerge_SP_Final/1-input_files/active_subnetworks/Stomach_cancer/MultipleModules/STAD_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'

for(database in databases){
  enrichr_analysis_multipleModules(root, disease_modules_filename, cancer_name , alg_name ,database, up, down, mode, delay)
}
print('enrichr_analysis_multipleModules finished for Core&Peel')



