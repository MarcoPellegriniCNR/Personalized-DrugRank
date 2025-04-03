# -*- coding: cp1252 -*-
# 2025-03-19
# This file holds the main mapping to store the directory names and
# filenames of data sets.
# The  data structure is a mapping indexed by a tupe of parameters.
# index includes the patient id, the drug data set, the algorithm for building the module covering of the
# active subnetwork.
# this file holds data from the active subnetwork computation with several algorithm and pvalue as mains score function.

subnetwork_data = {}
path_GEO = {}
path_CMAP = {}
path_L1000 = {}
path_DrugMatrix = {}

########################
# paziente TCGA-AZ-6605

path_GEO['TCGA-AZ-6605'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6605_colorectal_cancer/Colorectal_cancer_drugperturbationGEO'
path_CMAP['TCGA-AZ-6605'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6605_colorectal_cancer/Colorectal_cancer_CMAP'
path_L1000['TCGA-AZ-6605'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6605_colorectal_cancer/Colorectal_cancer_L1000'
path_DrugMatrix['TCGA-AZ-6605'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6605_colorectal_cancer/Colorectal_cancer_drugmatrix'

###########
# GEO data

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'GEO', '', 'crc')] = (path_GEO['TCGA-AZ-6605'] ,\
                                                                         'enrichR_ClustEx3800_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'GEO', '', 'crc')] = (path_GEO['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'Degas', 'GEO', '', 'crc')] = (path_GEO['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Degas_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'GEO', '', 'crc')] = (path_GEO['TCGA-AZ-6605'] ,\
                                                                         'enrichR_KPM_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'MD', 'GEO', '', 'crc')] = (path_GEO['TCGA-AZ-6605'] ,\
                                                                         'enrichR_MD_Drug_Perturbations_from_GEO_2014_perc.csv')


############
# Drugmarix data




subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AZ-6605'] ,\
                                                                         'enrichR_ClustEx3800_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_DrugMatrix_perc.csv')
# non existent
##subnetwork_data[('TCGA-AZ-6605', 'Degas', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AZ-6605'] ,\
##                                                                         'enrichR_Degas_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AZ-6605'] ,\
                                                                         'enrichR_KPM_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'MD', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AZ-6605'] ,\
                                                                         'enrichR_MD_DrugMatrix_perc.csv')


################
# CMAP Data
# up

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'Degas', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                        'enrichR_Degas_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_KPM_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'MD', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_MD_Old_CMAP_up_perc.csv')

# down

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'Degas', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                        'enrichR_Degas_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_KPM_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'MD', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AZ-6605'] ,\
                                                                         'enrichR_MD_Old_CMAP_down_perc.csv')

# up down

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'Degas', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'MD', 'CMAPupdown', '', 'crc')] = ()

################
# L1000 Data
# up

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_up_perc.csv')

# missing file
##subnetwork_data[('TCGA-AZ-6605', 'Degas', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
##                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'MD', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_up_perc.csv')

# down

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'Degas', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AZ-6605', 'MD', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AZ-6605'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_down_perc.csv')

# up down

subnetwork_data[('TCGA-AZ-6605', 'clex3800', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'cpd07j08', 'L1000updown', '', 'crc')] = ()

# missing one file 
##subnetwork_data[('TCGA-AZ-6605', 'Degas', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'KPM', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-AZ-6605', 'MD', 'L1000updown', '', 'crc')] = ()





########################
# paziente TCGA-A6-2671

path_GEO['TCGA-A6-2671'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-2671_colorectal_cancer/Colorectal_cancer_drugperturbationGEO'
path_CMAP['TCGA-A6-2671'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-2671_colorectal_cancer/Colorectal_cancer_CMAP'
path_L1000['TCGA-A6-2671'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-2671_colorectal_cancer/Colorectal_cancer_L1000'
path_DrugMatrix['TCGA-A6-2671'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-2671_colorectal_cancer/Colorectal_cancer_drugmatrix'

###########
# GEO data

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-2671'] ,\
                                                                         'enrichR_ClustEx3800_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-2671'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'Degas', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-2671'] ,\
                                                                         'enrichR_Degas_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'KPM', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-2671'] ,\
                                                                         'enrichR_KPM_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'MD', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-2671'] ,\
                                                                         'enrichR_MD_Drug_Perturbations_from_GEO_2014_perc.csv')


############
# Drugmarix data




subnetwork_data[('TCGA-A6-2671', 'clex3800', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-2671'] ,\
                                                                         'enrichR_ClustEx3800_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-2671'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_DrugMatrix_perc.csv')
# non existent
##subnetwork_data[('TCGA-A6-2671', 'Degas', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-2671'] ,\
##                                                                         'enrichR_Degas_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'KPM', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-2671'] ,\
                                                                         'enrichR_KPM_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'MD', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-2671'] ,\
                                                                         'enrichR_MD_DrugMatrix_perc.csv')


################
# CMAP Data
# up

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'Degas', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                        'enrichR_Degas_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'KPM', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_KPM_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'MD', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_MD_Old_CMAP_up_perc.csv')

# down

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'Degas', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                        'enrichR_Degas_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'KPM', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_KPM_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'MD', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-2671'] ,\
                                                                         'enrichR_MD_Old_CMAP_down_perc.csv')

# up down

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'Degas', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'KPM', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'MD', 'CMAPupdown', '', 'crc')] = ()

################
# L1000 Data
# up

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_up_perc.csv')

# missing file
##subnetwork_data[('TCGA-A6-2671', 'Degas', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
##                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'KPM', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'MD', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_up_perc.csv')

# down

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'Degas', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'KPM', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-2671', 'MD', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-2671'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_down_perc.csv')

# up down

subnetwork_data[('TCGA-A6-2671', 'clex3800', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'cpd07j08', 'L1000updown', '', 'crc')] = ()

# missing one file 
##subnetwork_data[('TCGA-A6-2671', 'Degas', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'KPM', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-2671', 'MD', 'L1000updown', '', 'crc')] = ()



########################
# paziente TCGA-AA-3517

path_GEO['TCGA-AA-3517'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AA-3517_colorectal_cancer/Colorectal_cancer_drugperturbationGEO'
path_CMAP['TCGA-AA-3517'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AA-3517_colorectal_cancer/Colorectal_cancer_CMAP'
path_L1000['TCGA-AA-3517'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AA-3517_colorectal_cancer/Colorectal_cancer_L1000'
path_DrugMatrix['TCGA-AA-3517'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AA-3517_colorectal_cancer/Colorectal_cancer_drugmatrix'

###########
# GEO data

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'GEO', '', 'crc')] = (path_GEO['TCGA-AA-3517'] ,\
                                                                         'enrichR_ClustEx3800_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'GEO', '', 'crc')] = (path_GEO['TCGA-AA-3517'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'Degas', 'GEO', '', 'crc')] = (path_GEO['TCGA-AA-3517'] ,\
                                                                         'enrichR_Degas_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'KPM', 'GEO', '', 'crc')] = (path_GEO['TCGA-AA-3517'] ,\
                                                                         'enrichR_KPM_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'MD', 'GEO', '', 'crc')] = (path_GEO['TCGA-AA-3517'] ,\
                                                                         'enrichR_MD_Drug_Perturbations_from_GEO_2014_perc.csv')


############
# Drugmarix data




subnetwork_data[('TCGA-AA-3517', 'clex3800', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AA-3517'] ,\
                                                                         'enrichR_ClustEx3800_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AA-3517'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_DrugMatrix_perc.csv')
# non existent
##subnetwork_data[('TCGA-AA-3517', 'Degas', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AA-3517'] ,\
##                                                                         'enrichR_Degas_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'KPM', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AA-3517'] ,\
                                                                         'enrichR_KPM_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'MD', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-AA-3517'] ,\
                                                                         'enrichR_MD_DrugMatrix_perc.csv')


################
# CMAP Data
# up

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'Degas', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                        'enrichR_Degas_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'KPM', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_KPM_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'MD', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_MD_Old_CMAP_up_perc.csv')

# down

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'Degas', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                        'enrichR_Degas_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'KPM', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_KPM_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'MD', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-AA-3517'] ,\
                                                                         'enrichR_MD_Old_CMAP_down_perc.csv')

# up down

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'Degas', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'KPM', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'MD', 'CMAPupdown', '', 'crc')] = ()

################
# L1000 Data
# up

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_up_perc.csv')

# missing file
##subnetwork_data[('TCGA-AA-3517', 'Degas', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
##                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'KPM', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'MD', 'L1000up', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_up_perc.csv')

# down

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'Degas', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'KPM', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-AA-3517', 'MD', 'L1000down', '', 'crc')] = (path_L1000['TCGA-AA-3517'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_down_perc.csv')

# up down

subnetwork_data[('TCGA-AA-3517', 'clex3800', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'cpd07j08', 'L1000updown', '', 'crc')] = ()

# missing one file 
##subnetwork_data[('TCGA-AA-3517', 'Degas', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'KPM', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-AA-3517', 'MD', 'L1000updown', '', 'crc')] = ()



########################
# paziente TCGA-HC-7742

path_GEO['TCGA-HC-7742'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-HC-7742_prostate_cancer/prostate_cancer_drugperturbationGEO'
path_CMAP['TCGA-HC-7742'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-HC-7742_prostate_cancer/prostate_cancer_CMAP'
path_L1000['TCGA-HC-7742'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-HC-7742_prostate_cancer/prostate_cancer_L1000'
path_DrugMatrix['TCGA-HC-7742'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-HC-7742_prostate_cancer/prostate_cancer_drugmatrix'

###########
# GEO data

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'GEO', '', 'prc')] = (path_GEO['TCGA-HC-7742'] ,\
                                                                         'enrichR_ClustEx3000_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'GEO', '', 'prc')] = (path_GEO['TCGA-HC-7742'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'Degas', 'GEO', '', 'prc')] = (path_GEO['TCGA-HC-7742'] ,\
                                                                         'enrichR_Degas_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'KPM', 'GEO', '', 'prc')] = (path_GEO['TCGA-HC-7742'] ,\
                                                                         'enrichR_KPM_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'MD', 'GEO', '', 'prc')] = (path_GEO['TCGA-HC-7742'] ,\
                                                                         'enrichR_MD_Drug_Perturbations_from_GEO_2014_perc.csv')


############
# Drugmarix data




subnetwork_data[('TCGA-HC-7742', 'clex3000', 'DM', '', 'prc')] = (path_DrugMatrix['TCGA-HC-7742'] ,\
                                                                         'enrichR_ClustEx3000_DrugMatrix_perc.csv')

# non existent
##subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'DM', '', 'prc')] = (path_DrugMatrix['TCGA-HC-7742'] ,\
##                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_DrugMatrix_perc.csv')
##
##subnetwork_data[('TCGA-HC-7742', 'Degas', 'DM', '', 'prc')] = (path_DrugMatrix['TCGA-HC-7742'] ,\
##                                                                         'enrichR_Degas_DrugMatrix_perc.csv')
##
##subnetwork_data[('TCGA-HC-7742', 'KPM', 'DM', '', 'prc')] = (path_DrugMatrix['TCGA-HC-7742'] ,\
##                                                                         'enrichR_KPM_DrugMatrix_perc.csv')
##
##subnetwork_data[('TCGA-HC-7742', 'MD', 'DM', '', 'prc')] = (path_DrugMatrix['TCGA-HC-7742'] ,\
##                                                                         'enrichR_MD_DrugMatrix_perc.csv')


################
# CMAP Data
# up

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'CMAPup', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_ClustEx3000_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'CMAPup', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'Degas', 'CMAPup', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                        'enrichR_Degas_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'KPM', 'CMAPup', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_KPM_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'MD', 'CMAPup', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_MD_Old_CMAP_up_perc.csv')

# down

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'CMAPdown', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_ClustEx3000_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'CMAPdown', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'Degas', 'CMAPdown', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                        'enrichR_Degas_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'KPM', 'CMAPdown', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_KPM_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'MD', 'CMAPdown', '', 'prc')] = (path_CMAP['TCGA-HC-7742'] ,\
                                                                         'enrichR_MD_Old_CMAP_down_perc.csv')

# up down

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'CMAPupdown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'CMAPupdown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'Degas', 'CMAPupdown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'KPM', 'CMAPupdown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'MD', 'CMAPupdown', '', 'prc')] = ()

################
# L1000 Data
# up

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'L1000up', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_ClustEx3000_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'L1000up', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_up_perc.csv')

# missing file
##subnetwork_data[('TCGA-HC-7742', 'Degas', 'L1000up', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
##                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'KPM', 'L1000up', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'MD', 'L1000up', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_up_perc.csv')

# down

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'L1000down', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_ClustEx3000_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'L1000down', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'Degas', 'L1000down', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'KPM', 'L1000down', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-HC-7742', 'MD', 'L1000down', '', 'prc')] = (path_L1000['TCGA-HC-7742'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_down_perc.csv')

# up down

subnetwork_data[('TCGA-HC-7742', 'clex3000', 'L1000updown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'cpd07j08', 'L1000updown', '', 'prc')] = ()

# missing one file 
##subnetwork_data[('TCGA-HC-7742', 'Degas', 'L1000updown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'KPM', 'L1000updown', '', 'prc')] = ()

subnetwork_data[('TCGA-HC-7742', 'MD', 'L1000updown', '', 'prc')] = ()



########################
# paziente TCGA-A6-5667

path_GEO['TCGA-A6-5667'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-5667_colorectal_cancer/Colorectal_cancer_drugperturbationGEO'
path_CMAP['TCGA-A6-5667'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-5667_colorectal_cancer/Colorectal_cancer_CMAP'
path_L1000['TCGA-A6-5667'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-5667_colorectal_cancer/Colorectal_cancer_L1000'
path_DrugMatrix['TCGA-A6-5667'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-5667_colorectal_cancer/Colorectal_cancer_drugmatrix'

###########
# GEO data

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-5667'] ,\
                                                                         'enrichR_ClustEx3800_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-5667'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'Degas', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-5667'] ,\
                                                                         'enrichR_Degas_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'KPM', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-5667'] ,\
                                                                         'enrichR_KPM_Drug_Perturbations_from_GEO_2014_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'MD', 'GEO', '', 'crc')] = (path_GEO['TCGA-A6-5667'] ,\
                                                                         'enrichR_MD_Drug_Perturbations_from_GEO_2014_perc.csv')


############
# Drugmarix data




subnetwork_data[('TCGA-A6-5667', 'clex3800', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-5667'] ,\
                                                                         'enrichR_ClustEx3800_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-5667'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_DrugMatrix_perc.csv')
# non existent
##subnetwork_data[('TCGA-A6-5667', 'Degas', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-5667'] ,\
##                                                                         'enrichR_Degas_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'KPM', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-5667'] ,\
                                                                         'enrichR_KPM_DrugMatrix_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'MD', 'DM', '', 'crc')] = (path_DrugMatrix['TCGA-A6-5667'] ,\
                                                                         'enrichR_MD_DrugMatrix_perc.csv')


################
# CMAP Data
# up

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'Degas', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                        'enrichR_Degas_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'KPM', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_KPM_Old_CMAP_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'MD', 'CMAPup', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_MD_Old_CMAP_up_perc.csv')

# down

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_ClustEx3800_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'Degas', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                        'enrichR_Degas_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'KPM', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_KPM_Old_CMAP_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'MD', 'CMAPdown', '', 'crc')] = (path_CMAP['TCGA-A6-5667'] ,\
                                                                         'enrichR_MD_Old_CMAP_down_perc.csv')

# up down

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'Degas', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'KPM', 'CMAPupdown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'MD', 'CMAPupdown', '', 'crc')] = ()

################
# L1000 Data
# up

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_up_perc.csv')

# missing file
##subnetwork_data[('TCGA-A6-5667', 'Degas', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
##                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'KPM', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_up_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'MD', 'L1000up', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_up_perc.csv')

# down

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_ClustEx3800_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'Degas', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                        'enrichR_Degas_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'KPM', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_KPM_LINCS_L1000_Chem_Pert_down_perc.csv')

subnetwork_data[('TCGA-A6-5667', 'MD', 'L1000down', '', 'crc')] = (path_L1000['TCGA-A6-5667'] ,\
                                                                         'enrichR_MD_LINCS_L1000_Chem_Pert_down_perc.csv')

# up down

subnetwork_data[('TCGA-A6-5667', 'clex3800', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'cpd07j08', 'L1000updown', '', 'crc')] = ()

# missing one file 
##subnetwork_data[('TCGA-A6-5667', 'Degas', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'KPM', 'L1000updown', '', 'crc')] = ()

subnetwork_data[('TCGA-A6-5667', 'MD', 'L1000updown', '', 'crc')] = ()
