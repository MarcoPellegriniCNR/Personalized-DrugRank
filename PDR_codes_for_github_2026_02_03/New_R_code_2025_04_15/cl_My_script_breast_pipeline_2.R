#!/usr/bin/env Rscript

#####################
# command line support

# Get command-line arguments.
# `trailingOnly = TRUE` ensures you only get the arguments *after* the script name
# and any R-specific options (like --args, though Rscript handles this implicitly).
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments were provided
if (length(args) == 0) {
  stop("Error: No parameters provided. Usage: ./my_executable_script.R <param1_value> <param2_value>", call. = FALSE)
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
source('My_functions_breast.R')


#---------------------------------------------------------------------------------
#                               MAIN
#----------------------------------------------------------------------------------


#background <- read.table('DrugMerge_SP_Final/1-input_files/symbol_entrezID_updated_table.txt', header = T)
background <- read.table('DrugMerge_SP_Final/1-input_files/symbol_entrezID_updated_table.txt', header = T)

###SPECIFY Drug databases###

#database <- c('LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down','Old_CMAP_down','Old_CMAP_up',
#              'DrugMatrix', 'Drug_Perturbations_from_GEO_2014')

database <- c('Drug_Perturbations_from_GEO_2014')


#database <- c('LINCS_L1000_Chem_Pert_up')
#database <- c('LINCS_L1000_Chem_Pert_down')###SPECIFY mode up down both or none###

#database <- c('Old_CMAP_down')
#mode = 'up'

#database <- c('Old_CMAP_up')
#mode = 'down'



#output <- paste0('DrugMerge_SP_Final/3-Drugs/') #, case_id, '_colorectal_cancer/')
output <- paste0('DrugMerge_SP_Final/3-Drugs/') #, case_id, '_colorectal_cancer/')

create_directory_breast(output, case_id)




####-------------------------------------------------------------------------------
#                           Breast cancer modules
#----------------------------------------------------------------------------------


# balanced at 4000 genes max for BRCA fluorouracil patients
#up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients/', case_id, '_4000_UP_gene_list.txt'))
#down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients/', case_id, '_4000_DOWN_gene_list.txt'))


# balanced at 4000 genes max for BRCA Doxorubicin Hydrochloride  patients
#up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients/', case_id, '_4000_UP_gene_list.txt'))
#down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients/', case_id, '_4000_DOWN_gene_list.txt'))


# balanced at 4000 genes max for BRCA Paclitaxel  patients
#up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients/', case_id, '_4000_UP_gene_list.txt'))
#down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients/', case_id, '_4000_DOWN_gene_list.txt'))


# balanced at 4000 genes max for BRCA Docetaxel patients
#up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients/', case_id, '_4000_UP_gene_list.txt'))
#down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients/', case_id, '_4000_DOWN_gene_list.txt'))

# balanced at 4000 genes max for BRCA Cyclophosphamide patients
#up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_cyclophosphamide_patients/', case_id, '_4000_UP_gene_list.txt'))
#down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_cyclophosphamide_patients/', case_id, '_4000_DOWN_gene_list.txt'))




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


# not used here 
#input <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/'
#input <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/'

# not used here
#output <- paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/')
#output <- paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/') 

print('Starting Modular EnrichR analysis')

if(dir.exists('DrugMerge_SP_Final/4-MultipleModules_output') == FALSE){
  dir.create('DrugMerge_SP_Final/4-MultipleModules_output')
}

# 1) perform the enrich analysis on the single modules
# Core&Peel
print('Starting multiple-module enrichment analysis for Core&Peel')

print('starting enrichr_analysis_multipleModules')
f <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Breast_cancer/MultipleModules/BRCA_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'

for(d in database){
  enrichr_analysis_multipleModules(f,'Breast_cancer','Core&Peel-d0.7-j0.8',d, up, down, mode, delay)
}
print('enrichr_analysis_multipleModules finished for Core&Peel')


####################
# merging of L1000 and CMAP up and down DB results - adding a suffix to the drung names..

flag_merge_L1000 = FALSE
flag_merge_CMAP = FALSE

if (flag_merge_L1000){
  print('merging up and down files for L1000')
  p <- paste0('DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/',case_id,'/Core&Peel-d0.7-j0.8/L1000')
  get_merge_enrichr('L1000',p)
  print('merging for L1000 finished')
  }

if (flag_merge_CMAP){
  print('merging up and down files for CMAP')
  p <- paste0('DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/',case_id,'/Core&Peel-d0.7-j0.8/CMAP')
  get_merge_enrichr('CMAP',p)
  print('merging for CMAP finished')
  }


print('Analysis for Core&Peel complete')


