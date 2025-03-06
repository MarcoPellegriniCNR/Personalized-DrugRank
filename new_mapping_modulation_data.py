# -*- coding: cp1252 -*-
# 2025-02-24
# This file holds the main mapping to store the directory names and
# filenames of data sets.
# The expred data structure is mapping indexed by a tupe of parameters.
# index includes te patient id, the drug data set, the algorithm for building the module covering of the
# active subnetwork.
# pilot is the first patient.


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

