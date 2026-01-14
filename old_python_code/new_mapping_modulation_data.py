# -*- coding: cp1252 -*-
# 2025-02-24
# This file holds the main mapping to store the directory names and
# filenames of data sets.
# The  data structure is a mapping indexed by a tupe of parameters.
# index includes the patient id, the drug data set, the algorithm for building the module covering of the
# active subnetwork.
# pilot is the first patient.

drugs_to_check = {}
drugs_to_check_lable = {}

##
drugs_to_check['pilot'] = []
drugs_to_check_lable['pilot'] = 'pilot'

##
drugs_to_check['TCGA-AZ-6600']  = ['Irinotecan','Fluorouracil','Leucovorin'] # record
drugs_to_check_lable['TCGA-AZ-6600'] = 'crc_drug_used_for_TCGA-AZ-6600'


##
drugs_to_check['TCGA-HC-7742'] = ['Leuprolide', # record
                  'Leuprorelin', 'Lupron', 'Eligard', 'Fensolvi', 'Viadur', 'Lupaneta', # aliases 
                  ]
drugs_to_check_lable['TCGA-HC-7742'] = 'prc_drug_used_for_TCGA-HC-7742'

##
drugs_to_check['TCGA-A6-2671'] = ['Conatumumab', 'Leucovorin', 'Fluorouracil', 'Bevacizumab', 'Oxaliplatin'] # record
drugs_to_check_lable['TCGA-A6-2671'] = 'crc_drug_used_for_TCGA-A6-2671'

##
drugs_to_check['TCGA-AA-3517'] = ['Leucovorin','Fluorouracil','Irinotecan','Bevacizumab','Oxaliplatin'] # record
drugs_to_check_lable['TCGA-AA-3517'] = 'crc_drug_used_for_TCGA-AA-3517'


##
drugs_to_check['TCGA-AZ-6605'] = []
drugs_to_check_lable['TCGA-AZ-6605'] = 'TCGA-AZ-6605'


##
drugs_to_check['TCGA-A6-5667'] = []
drugs_to_check_lable['TCGA-A6-5667'] = 'TCGA-A6-5667'



modulation_data = {}

########################
# pilot patient

############################
# Core and peel cpd07j08

#####
# GEO

modulation_data[('pilot', 'cpd07j08', 'GEO', 'filter' , 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('pilot', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                         'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')

#####
# CMAP

modulation_data[('pilot', 'cpd07j08', 'CMAP', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('pilot', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

#####
# DrugMatrix

modulation_data[('pilot', 'cpd07j08', 'DM', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('pilot', 'cpd07j08', 'DM', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noF.csv')
# problems with long filenames in windows..
#
#modulation_data[('pilot', 'cpd07j08', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
#                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####
# L1000

modulation_data[('pilot', 'cpd07j08', 'L1000', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('pilot', 'cpd07j08', 'L1000', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('pilot', 'MD', 'GEO', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_GEO.csv')

modulation_data[('pilot', 'MD', 'GEO', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                         'MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('pilot', 'MD', 'CMAP', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_CMAP.csv')

modulation_data[('pilot', 'MD', 'CMAP', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# DrugMatrix

modulation_data[('pilot', 'MD', 'DM', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('pilot', 'MD', 'DM', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                      'MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####
# L1000

modulation_data[('pilot', 'MD', 'L1000', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_L1000.csv')

modulation_data[('pilot', 'MD', 'L1000', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_L1000_noFiltering.csv')


########################
# paziente TCGA-AZ-6600



############################
# Core and peel cpd07j08

#####
# GEO

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_GEO.csv')

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                         'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_noFiltering.csv')


modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                         'TCGA-AZ-6600_cpd07j08_GEO_noFiltering.csv')

#####
# CMAP

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_CMAP.csv')



##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_CMAP_noFiltering.csv')

#####
# DrugMatrix

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_DrugMatrix.csv')



##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_DrugMatrix_noFiltering.csv')


#####
# L1000

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_L1000.csv')



##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_L1000_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-AZ-6600', 'MD', 'GEO', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'GEO', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                         'TCGA-AZ-6600_MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-AZ-6600', 'MD', 'CMAP', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'CMAP', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# DrugMatrix

modulation_data[('TCGA-AZ-6600', 'MD', 'DM', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'DM', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                      'TCGA-AZ-6600_MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-AZ-6600', 'MD', 'L1000', 'filter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'L1000', 'nofilter', 'crc')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_L1000_noFiltering.csv')


########################
# paziente TCGA-HC-7742

path_cp ={}
path_MD ={}
path_cp['TCGA-HC-7742'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-HC-7742_prostate_cancer'
path_MD['TCGA-HC-7742'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-HC-7742_prostate_cancer'



############################
# Core and peel cpd07j08

#####
# GEO

##modulation_data[('TCGA-HC-7742', 'cpd07j08', 'GEO', 'filter')] = (path_cp['TCGA-HC-7742'], \
##                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_GEO.csv')

modulation_data[('TCGA-HC-7742', 'cpd07j08', 'GEO', 'filter', 'prc')] = (path_cp['TCGA-HC-7742'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_GEO.csv')

modulation_data[('TCGA-HC-7742', 'cpd07j08', 'GEO', 'nofilter', 'prc')] = (path_cp['TCGA-HC-7742'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_GEO_noFiltering.csv')



######
# CMAP


modulation_data[('TCGA-HC-7742', 'cpd07j08', 'CMAP', 'filter', 'prc')] = (path_cp['TCGA-HC-7742'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_CMAP.csv')

modulation_data[('TCGA-HC-7742', 'cpd07j08', 'CMAP', 'nofilter', 'prc')] = (path_cp['TCGA-HC-7742'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_CMAP_noFiltering.csv')

######
# L1000


modulation_data[('TCGA-HC-7742', 'cpd07j08', 'L1000', 'filter', 'prc')] = (path_cp['TCGA-HC-7742'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_L1000.csv')

modulation_data[('TCGA-HC-7742', 'cpd07j08', 'L1000', 'nofilter', 'prc')] = (path_cp['TCGA-HC-7742'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Prostate_cancer_L1000_noFiltering.csv')



############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-HC-7742', 'MD', 'GEO', 'filter', 'prc')] = (path_MD['TCGA-HC-7742'], \
                                                           'MD_Prostate_cancer_GEO.csv')

modulation_data[('TCGA-HC-7742', 'MD', 'GEO', 'nofilter', 'prc')] = (path_MD['TCGA-HC-7742'] , \
                                                           'MD_Prostate_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-HC-7742', 'MD', 'CMAP', 'filter', 'prc')] = (path_MD['TCGA-HC-7742'], \
                                                           'MD_Prostate_cancer_CMAP.csv')

modulation_data[('TCGA-HC-7742', 'MD', 'CMAP', 'nofilter', 'prc')] = (path_MD['TCGA-HC-7742'] , \
                                                           'MD_Prostate_cancer_CMAP_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-HC-7742', 'MD', 'L1000', 'filter', 'prc')] = (path_MD['TCGA-HC-7742'], \
                                                           'MD_Prostate_cancer_L1000.csv')

modulation_data[('TCGA-HC-7742', 'MD', 'L1000', 'nofilter', 'prc')] = (path_MD['TCGA-HC-7742'] , \
                                                           'MD_Prostate_cancer_L1000_noFiltering.csv')

#################################################################################################
########################
# paziente TCGA-A6-2671


path_cp['TCGA-A6-2671'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-2671_colorectal_cancer'
path_MD['TCGA-A6-2671'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-2671_colorectal_cancer'

######
# core and peel

#####
# GEO


modulation_data[('TCGA-A6-2671', 'cpd07j08', 'GEO', 'filter', 'crc')] = (path_cp['TCGA-A6-2671'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-A6-2671', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-A6-2671'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')



######
# CMAP


modulation_data[('TCGA-A6-2671', 'cpd07j08', 'CMAP', 'filter', 'crc')] = (path_cp['TCGA-A6-2671'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-A6-2671', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-A6-2671'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

######
# L1000


modulation_data[('TCGA-A6-2671', 'cpd07j08', 'L1000', 'filter', 'crc')] = (path_cp['TCGA-A6-2671'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-A6-2671', 'cpd07j08', 'L1000', 'nofilter', 'crc')] = (path_cp['TCGA-A6-2671'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')

#######
# Drug matrix


modulation_data[('TCGA-A6-2671', 'cpd07j08', 'DM', 'filter', 'crc')] = (path_cp['TCGA-A6-2671'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-A6-2671', 'cpd07j08', 'DM', 'nofilter', 'crc')] = (path_cp['TCGA-A6-2671'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-A6-2671', 'MD', 'GEO', 'filter', 'crc')] = (path_MD['TCGA-A6-2671'], \
                                                           'MD_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-A6-2671', 'MD', 'GEO', 'nofilter', 'crc')] = (path_MD['TCGA-A6-2671'] , \
                                                           'MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-A6-2671', 'MD', 'CMAP', 'filter', 'crc')] = (path_MD['TCGA-A6-2671'], \
                                                           'MD_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-A6-2671', 'MD', 'CMAP', 'nofilter', 'crc')] = (path_MD['TCGA-A6-2671'] , \
                                                           'MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-A6-2671', 'MD', 'L1000', 'filter', 'crc')] = (path_MD['TCGA-A6-2671'], \
                                                           'MD_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-A6-2671', 'MD', 'L1000', 'nofilter', 'crc')] = (path_MD['TCGA-A6-2671'] , \
                                                           'MD_Colorectal_cancer_L1000_noFiltering.csv')
####
# Drug marix

modulation_data[('TCGA-A6-2671', 'MD', 'DM', 'filter', 'crc')] = (path_MD['TCGA-A6-2671'], \
                                                           'MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-A6-2671', 'MD', 'DM', 'nofilter', 'crc')] = (path_MD['TCGA-A6-2671'] , \
                                                           'MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#################################################################################################
########################
# paziente TCGA-AA-3517


path_cp['TCGA-AA-3517'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AA-3517_colorectal_cancer'
path_MD['TCGA-AA-3517'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AA-3517_colorectal_cancer'



######
# core and peel

#####
# GEO


modulation_data[('TCGA-AA-3517', 'cpd07j08', 'GEO', 'filter', 'crc')] = (path_cp['TCGA-AA-3517'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AA-3517', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-AA-3517'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')



######
# CMAP


modulation_data[('TCGA-AA-3517', 'cpd07j08', 'CMAP', 'filter', 'crc')] = (path_cp['TCGA-AA-3517'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AA-3517', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-AA-3517'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

######
# L1000


modulation_data[('TCGA-AA-3517', 'cpd07j08', 'L1000', 'filter', 'crc')] = (path_cp['TCGA-AA-3517'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AA-3517', 'cpd07j08', 'L1000', 'nofilter', 'crc')] = (path_cp['TCGA-AA-3517'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')

#######
# Drug matrix


modulation_data[('TCGA-AA-3517', 'cpd07j08', 'DM', 'filter', 'crc')] = (path_cp['TCGA-AA-3517'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AA-3517', 'cpd07j08', 'DM', 'nofilter', 'crc')] = (path_cp['TCGA-AA-3517'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-AA-3517', 'MD', 'GEO', 'filter', 'crc')] = (path_MD['TCGA-AA-3517'], \
                                                           'MD_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AA-3517', 'MD', 'GEO', 'nofilter', 'crc')] = (path_MD['TCGA-AA-3517'] , \
                                                           'MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-AA-3517', 'MD', 'CMAP', 'filter', 'crc')] = (path_MD['TCGA-AA-3517'], \
                                                           'MD_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AA-3517', 'MD', 'CMAP', 'nofilter', 'crc')] = (path_MD['TCGA-AA-3517'] , \
                                                           'MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-AA-3517', 'MD', 'L1000', 'filter', 'crc')] = (path_MD['TCGA-AA-3517'], \
                                                           'MD_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AA-3517', 'MD', 'L1000', 'nofilter', 'crc')] = (path_MD['TCGA-AA-3517'] , \
                                                           'MD_Colorectal_cancer_L1000_noFiltering.csv')
####
# Drug marix

modulation_data[('TCGA-AA-3517', 'MD', 'DM', 'filter', 'crc')] = (path_MD['TCGA-AA-3517'], \
                                                           'MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AA-3517', 'MD', 'DM', 'nofilter', 'crc')] = (path_MD['TCGA-AA-3517'] , \
                                                           'MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')



##########################

########################
# paziente TCGA-AZ-6605


path_cp['TCGA-AZ-6605'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6605_colorectal_cancer'
path_MD['TCGA-AZ-6605'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6605_colorectal_cancer'



######
# core and peel

#####
# GEO


modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'GEO', 'filter', 'crc')] = (path_cp['TCGA-AZ-6605'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6605'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')



######
# CMAP


modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'CMAP', 'filter', 'crc')] = (path_cp['TCGA-AZ-6605'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6605'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

######
# L1000


modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'L1000', 'filter', 'crc')] = (path_cp['TCGA-AZ-6605'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'L1000', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6605'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')

#######
# Drug matrix


modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'DM', 'filter', 'crc')] = (path_cp['TCGA-AZ-6605'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AZ-6605', 'cpd07j08', 'DM', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6605'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-AZ-6605', 'MD', 'GEO', 'filter', 'crc')] = (path_MD['TCGA-AZ-6605'], \
                                                           'MD_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AZ-6605', 'MD', 'GEO', 'nofilter', 'crc')] = (path_MD['TCGA-AZ-6605'] , \
                                                           'MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-AZ-6605', 'MD', 'CMAP', 'filter', 'crc')] = (path_MD['TCGA-AZ-6605'], \
                                                           'MD_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AZ-6605', 'MD', 'CMAP', 'nofilter', 'crc')] = (path_MD['TCGA-AZ-6605'] , \
                                                           'MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-AZ-6605', 'MD', 'L1000', 'filter', 'crc')] = (path_MD['TCGA-AZ-6605'], \
                                                           'MD_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AZ-6605', 'MD', 'L1000', 'nofilter', 'crc')] = (path_MD['TCGA-AZ-6605'] , \
                                                           'MD_Colorectal_cancer_L1000_noFiltering.csv')
####
# Drug marix

modulation_data[('TCGA-AZ-6605', 'MD', 'DM', 'filter', 'crc')] = (path_MD['TCGA-AZ-6605'], \
                                                           'MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AZ-6605', 'MD', 'DM', 'nofilter', 'crc')] = (path_MD['TCGA-AZ-6605'] , \
                                                           'MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')




########################
# paziente TCGA-A6-5667


path_cp['TCGA-A6-5667'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-5667_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp'
path_MD['TCGA-A6-5667'] = '../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-A6-5667_colorectal_cancer/Drug-modulation-data/colorectal-cancer-MD'

######
# core and peel

#####
# GEO


modulation_data[('TCGA-A6-5667', 'cpd07j08', 'GEO', 'filter', 'crc')] = (path_cp['TCGA-A6-5667'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-A6-5667', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-A6-5667'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')



######
# CMAP


modulation_data[('TCGA-A6-5667', 'cpd07j08', 'CMAP', 'filter', 'crc')] = (path_cp['TCGA-A6-5667'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-A6-5667', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-A6-5667'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

######
# L1000


modulation_data[('TCGA-A6-5667', 'cpd07j08', 'L1000', 'filter', 'crc')] = (path_cp['TCGA-A6-5667'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-A6-5667', 'cpd07j08', 'L1000', 'nofilter', 'crc')] = (path_cp['TCGA-A6-5667'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')

#######
# Drug matrix


modulation_data[('TCGA-A6-5667', 'cpd07j08', 'DM', 'filter', 'crc')] = (path_cp['TCGA-A6-5667'], \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-A6-5667', 'cpd07j08', 'DM', 'nofilter', 'crc')] = (path_cp['TCGA-A6-5667'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-A6-5667', 'MD', 'GEO', 'filter', 'crc')] = (path_MD['TCGA-A6-5667'], \
                                                           'MD_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-A6-5667', 'MD', 'GEO', 'nofilter', 'crc')] = (path_MD['TCGA-A6-5667'] , \
                                                           'MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-A6-5667', 'MD', 'CMAP', 'filter', 'crc')] = (path_MD['TCGA-A6-5667'], \
                                                           'MD_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-A6-5667', 'MD', 'CMAP', 'nofilter', 'crc')] = (path_MD['TCGA-A6-5667'] , \
                                                           'MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-A6-5667', 'MD', 'L1000', 'filter', 'crc')] = (path_MD['TCGA-A6-5667'], \
                                                           'MD_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-A6-5667', 'MD', 'L1000', 'nofilter', 'crc')] = (path_MD['TCGA-A6-5667'] , \
                                                           'MD_Colorectal_cancer_L1000_noFiltering.csv')
####
# Drug marix

modulation_data[('TCGA-A6-5667', 'MD', 'DM', 'filter', 'crc')] = (path_MD['TCGA-A6-5667'], \
                                                           'MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-A6-5667', 'MD', 'DM', 'nofilter', 'crc')] = (path_MD['TCGA-A6-5667'] , \
                                                           'MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####################
