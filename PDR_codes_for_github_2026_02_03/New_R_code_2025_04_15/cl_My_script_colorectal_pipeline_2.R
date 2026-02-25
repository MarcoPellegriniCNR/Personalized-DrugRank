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
source('My_functions_colorectal.R')



#---------------------------------------------------------------------------------
#                               MAIN
#----------------------------------------------------------------------------------

###SPECIFY CASE ID###
#case_id <- "TCGA-A6-5667"
#case_id <- "TCGA-AZ-6600"
#case_id <- "TCGA-A6-2671"
#case_id <- "TCGA-AA-3517"
#case_id <- "TCGA-AZ-6605"

# new series
#case_id <- "TCGA-AA-3679"
#case_id  <- "TCGA-A6-2671"
#case_id <- "TCGA-AA-3869"
#case_id <- "TCGA-CM-5868"
#case_id  <- "TCGA-NH-A6GA"
#case_id  <- "TCGA-AY-4070"
#case_id <- "TCGA-AA-A01R"
#case_id <- "TCGA-AA-A02F"
#case_id <- "TCGA-F4-6806"
#case_id <- "TCGA-AA-3973"
#case_id <- "TCGA-AA-3548"

#case_id <- "TCGA-A6-5664"
#case_id <- "TCGA-A6-2683"
#case_id <- "TCGA-AA-A01D"
#case_id <- "TCGA-CM-4748"
#case_id <- "TCGA-CM-6675"
#case_id <- "TCGA-CM-6678"
#case_id <- "TCGA-G4-6297"
#case_id <- "TCGA-G4-6303"
#case_id <- "TCGA-AA-3494"

#case_id <- "TCGA-AZ-4682"

# new 6 rectal cancer patients with irinotecan

#'TCGA-F5-6571-01A', 'TCGA-DC-6156-01A', 'TCGA-AG-A02G-01A', 'TCGA-AG-3887-01A', 'TCGA-AG-3731-01A', 'TCGA-AG-3584-01A'

#case_id <- 'TCGA-F5-6571'
#case_id <-  'TCGA-DC-6156'
#case_id  <-  'TCGA-AG-A02G'
#case_id <- "TCGA-AG-3887"
#case_id <- "TCGA-AG-3731"
#case_id <- "TCGA-AG-3584"

########
# fluorouracil patients

#case_id <- "TCGA-A6-A566"
#case_id <- "TCGA-A6-A56B"
# TCGA-A6-A5ZU
# TCGA-AA-3971
# TCGA-AA-A00J
# TCGA-AA-A00Q
# TCGA-AA-A010
# TCGA-AA-A01F
# TCGA-AA-A01K
# TCGA-AA-A01T
# TCGA-AD-6965
# TCGA-AF-6136
# TCGA-AF-A56L
# TCGA-AG-3593
# TCGA-AG-3600
# TCGA-AG-3609
# TCGA-AG-4022
# TCGA-AG-A01W
# TCGA-AG-A01Y
# TCGA-AG-A020
# TCGA-AH-6544
# TCGA-AH-6643
# TCGA-AH-6903
# TCGA-AY-A8YK
# TCGA-AZ-6606
# TCGA-CA-5254
# TCGA-CA-5255
# TCGA-CA-5256
# TCGA-CA-5796
# TCGA-CA-6715
# TCGA-CA-6716
# TCGA-CA-6717
# TCGA-CK-4950
# TCGA-F4-6569
# TCGA-F4-6856
# TCGA-F5-6702
# TCGA-F5-6863
# TCGA-F5-6864
# TCGA-NH-A50U
# TCGA-NH-A50V
# TCGA-NH-A6GB
# TCGA-NH-A8F8
# TCGA-SS-A7HO

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

create_directory_colorectal(output, case_id)




####-------------------------------------------------------------------------------
#                           Colorectal cancer
#----------------------------------------------------------------------------------





# balanced at 4000 genes max for crc irinotecan  patients
#up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_bis/', case_id, '_4000_UP_gene_list.txt'))
#down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_bis/', case_id, '_4000_DOWN_gene_list.txt'))

# balanced at 4000 genes max for fluorouracil patients
up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients/', case_id, '_4000_UP_gene_list.txt'))
down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients/', case_id, '_4000_DOWN_gene_list.txt'))





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
f <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'
for(d in database){
  enrichr_analysis_multipleModules(f,'Colorectal_cancer','Core&Peel-d0.7-j0.8',d, up, down, mode, delay)
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


