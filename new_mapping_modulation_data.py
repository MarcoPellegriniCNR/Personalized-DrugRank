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

modulation_data[('pilot', 'cpd07j08', 'GEO', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('pilot', 'cpd07j08', 'GEO', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                         'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')

#####
# CMAP

modulation_data[('pilot', 'cpd07j08', 'CMAP', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('pilot', 'cpd07j08', 'CMAP', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

#####
# DrugMatrix

modulation_data[('pilot', 'cpd07j08', 'DM', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('pilot', 'cpd07j08', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noF.csv')
# problems with long filenames in windows..
#
#modulation_data[('pilot', 'cpd07j08', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
#                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####
# L1000

modulation_data[('pilot', 'cpd07j08', 'L1000', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('pilot', 'cpd07j08', 'L1000', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
                                                           'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('pilot', 'MD', 'GEO', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_GEO.csv')

modulation_data[('pilot', 'MD', 'GEO', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                         'MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('pilot', 'MD', 'CMAP', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_CMAP.csv')

modulation_data[('pilot', 'MD', 'CMAP', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# DrugMatrix

modulation_data[('pilot', 'MD', 'DM', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('pilot', 'MD', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                      'MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####
# L1000

modulation_data[('pilot', 'MD', 'L1000', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_L1000.csv')

modulation_data[('pilot', 'MD', 'L1000', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_MD/colorectal-cancer-MD' , \
                                                           'MD_Colorectal_cancer_L1000_noFiltering.csv')


########################
# paziente TCGA-AZ-6600



############################
# Core and peel cpd07j08

#####
# GEO

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_GEO.csv')

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                         'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_noFiltering.csv')


modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'GEO', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                         'TCGA-AZ-6600_cpd07j08_GEO_noFiltering.csv')

#####
# CMAP

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_CMAP.csv')



##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'CMAP', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_CMAP_noFiltering.csv')

#####
# DrugMatrix

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_DrugMatrix.csv')



##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_DrugMatrix_noFiltering.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_DrugMatrix_noFiltering.csv')


#####
# L1000

##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_L1000.csv')



##modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp' , \
##                                                           'TCGA-AZ-6600_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_L1000_noFiltering.csv')

modulation_data[('TCGA-AZ-6600', 'cpd07j08', 'L1000', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-cp' , \
                                                           'TCGA-AZ-6600_cpd07j08_L1000_noFiltering.csv')


############################
# Module Discovery MD


#####
# GEO

modulation_data[('TCGA-AZ-6600', 'MD', 'GEO', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_GEO.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'GEO', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                         'TCGA-AZ-6600_MD_Colorectal_cancer_GEO_noFiltering.csv')


#####
# CMAP

modulation_data[('TCGA-AZ-6600', 'MD', 'CMAP', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_CMAP.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'CMAP', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_CMAP_noFiltering.csv')


#####
# DrugMatrix

modulation_data[('TCGA-AZ-6600', 'MD', 'DM', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_DrugMatrix.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'DM', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                      'TCGA-AZ-6600_MD_Colorectal_cancer_DrugMatrix_noFiltering.csv')


#####
# L1000
modulation_data[('TCGA-AZ-6600', 'MD', 'L1000', 'filter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_L1000.csv')

modulation_data[('TCGA-AZ-6600', 'MD', 'L1000', 'nofilter')] = ('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/TCGA-AZ-6600_colorectal_cancer/TCGA-AZ-6600_colorectal-cancer-MD' , \
                                                           'TCGA-AZ-6600_MD_Colorectal_cancer_L1000_noFiltering.csv')


