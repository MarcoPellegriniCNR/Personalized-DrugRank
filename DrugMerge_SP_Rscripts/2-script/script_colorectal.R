###### set your working directory ######
setwd('/home/fiamma/DrugMerge/')
library(limma)
library(enrichR)
library(readr)

source('DrugMerge_SP_Final/2-script/functions.R')

#---------------------------------------------------------------------------------
#                               MAIN
#----------------------------------------------------------------------------------

###SPECIFY CASE ID###
case_id <- "TCGA-A6-5667"

background <- read.table('DrugMerge_SP_Final/1-input_files/symbol_entrezID_updated_table.txt', header = T)

database <- c('LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down','Old_CMAP_down','Old_CMAP_up',
              'DrugMatrix', 'Drug_Perturbations_from_GEO_2014')

output <- paste0('DrugMerge_SP_Final/3-Drugs/') #, case_id, '_colorectal_cancer/')
create_directory_colorectal(output, case_id)

####-------------------------------------------------------------------------------
#                           Colorectal cancer
#----------------------------------------------------------------------------------
up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_', case_id, '_UP_gene_list.txt'))
down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_', case_id, '_DOWN_gene_list.txt'))


input <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/'
output <- paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/') 
print('Starting EnrichR analysis')


# Core&Peel
print('Core&Peel')
file <- 'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Colorectal.txt'
enrichr_analysis(input,file,'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest','Colorectal_cancer',
                 database,background,up,down,output)
print('Finished')

# MD
print('ModuleDiscoverer')
file <- 'MD_Colorectal.txt'
enrichr_analysis(input,file,'MD','Colorectal_cancer',
                 database,background,up,down,output)
print('Finished')

# ClustEx
print('ClustEx3800')
file <- 'ClustEx3800_Colorectal.txt'
enrichr_analysis(input,file,'ClustEx3800','Colorectal_cancer',
                 database,background,up,down,output)
print('Finished')

# KPM
print('KeyPathwayMiner')
file <- 'KPM_Colorectal.txt'
enrichr_analysis(input,file,'KPM','Colorectal_cancer',
                 database,background,up,down,output)
print('Finished')

# Degas
print('DEGAS')
file <- 'Degas_Colorectal.txt'
enrichr_analysis(input,file,'Degas','Colorectal_cancer',
                 database,background,up,down,output)
print('Finished')

################################### Modulation drug ########################

if(dir.exists('DrugMerge_SP_Final/4-MultipleModules_output') == FALSE){
  dir.create('DrugMerge_SP_Final/4-MultipleModules_output')
}

# 1) perform the enrich analysis on the single modules
# Core&Peel
print('Starting multiple-module enrichment analysis for Core&Peel')

print('starting enrichr_analysis_multipleModules')
f <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'
for(d in database){
  enrichr_analysis_multipleModules(f,'Colorectal_cancer','Core&Peel-d0.7-j0.8',d)
}
print('enrichr_analysis_multipleModules finished')


print('merging up and down files for L1000')
p <- paste0('DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/',case_id,'/Core&Peel-d0.7-j0.8/L1000')
get_merge_enrichr('L1000',p)
print('merging for L1000 finished')


print('merging up and down files for CMAP')
p <- paste0('DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/',case_id,'/Core&Peel-d0.7-j0.8/CMAP')
get_merge_enrichr('CMAP',p)
print('merging for CMAP finished')


print('Analysis for Core&Peel complete')

# MD
print('Starting multiple-module enrichment analysis for MD')

print('starting enrichr_analysis_multipleModules')
f <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/MD_MultipleModules_Colorectal_cancer.txt'
for(d in database){
  enrichr_analysis_multipleModules(f,'Colorectal_cancer','MD',d)
}
print('enrichr_analysis_multipleModules finished')

print('merging up and down files for L1000')
p <- paste0('DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/', case_id, '/MD/L1000')
get_merge_enrichr('L1000',p)
print('merging for L1000 finished')

print('merging up and down files for CMAP')
p <- paste0('DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/',case_id,'MD/CMAP')
get_merge_enrichr('CMAP',p)
print('merging for CMAP finished')

print('Analysis for MD complete')

##### 2) get table without filtering #########

up <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_', case_id, '_UP_gene_list.txt'))
up <- up$V1
down <- read.table(paste0('DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_', case_id, '_DOWN_gene_list.txt'))
down <- down$V1
database1 <- c('L1000','CMAP','drugperturbationGEO', 'drugmatrix')

# Core&Peel
print("Creating table without filtering for Core&Peel")
for(d in database1){
  enrichr <- read_enrichr(paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/'),'Core&Peel','Colorectal_cancer',d, param = '-d0.7-f1-j0.8')
  path <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/'
  table <- get_modulation_table(enrichr,'Core&Peel',path,up,down)
  path_input <- paste0("DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/", case_id, "/")
  path_output <- paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp/')
  table1 <- add_pvalue(path_input,'Core&Peel',d,table,path_output,'Colorectal_cancer',param = '-d0.7-f1-j0.8')
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
  table1 <- read.csv(paste0(path_output,'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8','_','Colorectal_cancer','_',d1,'_noFiltering.csv'))
  get_modulation_filter(table1,d,'Core&Peel',path_output,'Colorectal_cancer',param = '-d0.7-f1-j0.8')
}
print("Finished")

# MD
print("Creating table without filtering for MD")
for(d in database1){
  enrichr <- read_enrichr(paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/'),'MD','Colorectal_cancer',d)
  path <- 'DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/'
  table <- get_modulation_table(enrichr,'MD',path,up,down)
  path_input <- paste0("DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/", case_id, "/")
  path_output <- paste0('DrugMerge_SP_Final/3-Drugs/', case_id, '_colorectal_cancer/Drug-modulation-data/colorectal-cancer-MD/')
  table1 <- add_pvalue(path_input,'MD',d,table,path_output,'Colorectal_cancer')
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
  table1 <- read.csv(paste0(path_output,'MD','_','Colorectal_cancer','_',d1,'_noFiltering.csv'))
  get_modulation_filter(table1,d,'MD',path_output,'Colorectal_cancer')
}
print("Finished")
