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
drugs_to_check['TCGA-AZ-6600-bis']  = ['Irinotecan','Fluorouracil','Leucovorin'] # record
drugs_to_check_lable['TCGA-AZ-6600-bis'] = 'crc_drug_used_for_TCGA-AZ-6600-bis'


##
drugs_to_check['TCGA-HC-7742'] = ['Leuprolide', # record
                  'Leuprorelin', 'Lupron', 'Eligard', 'Fensolvi', 'Viadur', 'Lupaneta', # aliases 
                  ]
drugs_to_check_lable['TCGA-HC-7742'] = 'prc_drug_used_for_TCGA-HC-7742'



##
##drugs_to_check['TCGA-A6-2671'] = ['Conatumumab', 'Leucovorin', 'Fluorouracil', 'Bevacizumab', 'Oxaliplatin'] # record dor ony first cyecle
##drugs_to_check_lable['TCGA-A6-2671'] = 'crc_drug_used_for_TCGA-A6-2671 in first cycle'
drugs_to_check['TCGA-A6-2671-bis'] = ['Conatumumab', 'Leucovorin', 'Fluorouracil', 'Bevacizumab', 'Oxaliplatin'] # record
drugs_to_check_lable['TCGA-A6-2671-bis'] = 'crc_drug_used_for_TCGA-A6-2671-bis'


drugs_to_check['TCGA-A6-2671'] = ['Conatumumab', 'Leucovorin', 'Fluorouracil', 'Bevacizumab', 'Oxaliplatin', 'Irinotecan', 'Cetuximab'] # record for all cycles
drugs_to_check_lable['TCGA-A6-2671'] = 'crc_drug_used_for_TCGA-A6-2671 in all cycles '


##
drugs_to_check['TCGA-AA-3517'] = ['Leucovorin','Fluorouracil','Irinotecan','Bevacizumab','Oxaliplatin'] # record
drugs_to_check_lable['TCGA-AA-3517'] = 'crc_drug_used_for_TCGA-AA-3517'
drugs_to_check['TCGA-AA-3517-bis'] = ['Leucovorin','Fluorouracil','Irinotecan','Bevacizumab','Oxaliplatin'] # record
drugs_to_check_lable['TCGA-AA-3517-bis'] = 'crc_drug_used_for_TCGA-AA-3517-bis'


##
drugs_to_check['TCGA-AZ-6605'] = []
drugs_to_check_lable['TCGA-AZ-6605'] = 'TCGA-AZ-6605'
drugs_to_check['TCGA-AZ-6605-bis'] = []
drugs_to_check_lable['TCGA-AZ-6605-bis'] = 'crc_drug_used_forTCGA-AZ-6605-bis'


##
drugs_to_check['TCGA-A6-5667'] = []
drugs_to_check_lable['TCGA-A6-5667'] = 'TCGA-A6-5667'
drugs_to_check['TCGA-A6-5667-bis'] = []
drugs_to_check_lable['TCGA-A6-5667-bis'] = 'crc_drug_used_forTCGA-A6-5667-bis'

####################
# Irinotecan series
##
drugs_to_check['TCGA-AA-3679'] = ['Leucovorin', 'Fluorouracil', 'Irinotecan', 'Bevacizumab']
drugs_to_check_lable['TCGA-AA-3679'] = 'crc_drug_used_for_TCGA-A6-3679'


drugs_to_check['TCGA-AA-3869'] =['Capecitabine', 'Leucovorin', 'Fluorouracil', 'Irinotecan', 'Bevacizumab', 'Cetuximab']
drugs_to_check_lable['TCGA-AA-3869'] = 'crc_drug_used_for_TCGA-A6-3869'


drugs_to_check['TCGA-CM-5868'] = ['Fluorouracil', 'Oxaliplatin', 'Irinotecan', 'Dexamethasone', 'Leucovorin']
drugs_to_check_lable['TCGA-CM-5868']     =  'crc_drug_used_for_TCGA-CM-5868'


drugs_to_check['TCGA-NH-A6GA'] = ['Oxaliplatin', 'Bevacizumab', 'Leucovorin', 'Fluorouracil', 'Irinotecan']
drugs_to_check_lable['TCGA-NH-A6GA']     =  'crc_drug_used_for_TCGA-NH-A6GA'


drugs_to_check['TCGA-AY-4070'] = ['Fluorouracil', 'Oxaliplatin', 'Capecitabine', 'Irinotecan']
drugs_to_check_lable['TCGA-AY-4070']     =  'crc_drug_used_for_TCGA-AY-4070'


drugs_to_check['TCGA-AA-A01R'] = ['Fluorouracil', 'Leucovorin', 'Irinotecan']
drugs_to_check_lable['TCGA-AA-A01R']     =  'crc_drug_used_for_TCGA-AA-A01R'


drugs_to_check['TCGA-AA-A02F']  = ['Leucovorin', 'Fluorouracil', 'Irinotecan', 'Bevacizumab', 'Oxaliplatin']
drugs_to_check_lable['TCGA-AA-A02F']     =  'crc_drug_used_for_TCGA-AA-A02F'

drugs_to_check['TCGA-AA-3973'] = ['Fluorouracil', 'Leucovorin', 'Irinotecan']
drugs_to_check_lable['TCGA-AA-3973'] =  'crc_drug_used_for_TCGA-AA-3973'

drugs_to_check['TCGA-AA-3548'] = ['Fluorouracil', 'Leucovorin', 'Irinotecan']
drugs_to_check_lable['TCGA-AA-3548'] =  'crc_drug_used_for_TCGA-AA-3548'

drugs_to_check['TCGA-F4-6806'] = ['Irinotecan', 'Leucovorin', 'Fluorouracil', 'Mitomycin', 'Capecitabine', 'Oxaliplatin', 'Bevacizumab']
drugs_to_check_lable['TCGA-F4-6806'] =  'crc_drug_used_for_TCGA-F4-6806'

drugs_to_check['TCGA-A6-5664'] = ['Oxaliplatin', 'Bevacizumab', 'Fluorouracil', 'Leucovorin', 'Irinotecan', 'Bevacizumab']
drugs_to_check_lable['TCGA-A6-5664'] =  'crc_drug_used_for_TCGA-A6-5664'

drugs_to_check['TCGA-A6-2683'] =  ['Fluorouracil' , 'Leucovorin' , 'Irinotecan' ,'Bevacizumab', 'Oxaliplatin', 'Mitomycin']
drugs_to_check_lable['TCGA-A6-2683'] =  'crc_drug_used_for_TCGA-A6-2683'


drugs_to_check['TCGA-AA-A01D'] =  ['Capecitabine', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab', 'Mitomycin',' Panitumumab']
drugs_to_check_lable['TCGA-AA-A01D'] =  'crc_drug_used_for_TCGA-AA-A01D'

drugs_to_check['TCGA-CM-4748'] =  ['Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab', 'Cetuximab']
drugs_to_check_lable['TCGA-CM-4748'] =  'crc_drug_used_for_TCGA-CM-4748'

drugs_to_check['TCGA-CM-6675'] =  ['Fluorouracil', 'Leucovorin',  'Oxaliplatin', 'Irinotecan', 'Pegfilgrastim', 'Cetuximab']
drugs_to_check_lable['TCGA-CM-6675'] =  'crc_drug_used_for_TCGA-CM-6675'


drugs_to_check['TCGA-CM-6678'] =  ['Fluorouracil', 'Oxaliplatin', 'Irinotecan']
drugs_to_check_lable['TCGA-CM-6678'] =  'crc_drug_used_for_TCGA-CM-6678'


drugs_to_check['TCGA-G4-6297'] =  ['Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab', 'Cetuximab','Capecitabine']
drugs_to_check_lable['TCGA-G4-6297'] =  'crc_drug_used_for_TCGA-G4-6297'


drugs_to_check['TCGA-G4-6303'] =  ['Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab', 'Cetuximab']
drugs_to_check_lable['TCGA-G4-6303'] =  'crc_drug_used_for_TCGA-G4-6303'


drugs_to_check['TTCGA-AA-3494'] =  ['Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Cetuximab','Capecitabine']
drugs_to_check_lable['TCGA-AA-3494'] =  'crc_drug_used_for_TCGA-AA-3494'


drugs_to_check['TCGA-AZ-4682'] = ['Irinoteca','Leucovorin', 'Fluorouracil', 'Bevacizumab']
drugs_to_check_lable['TCGA-AZ-4682'] =  'crc_drug_used_for_TCGA-AZ-4682'


drugs_to_check['TCGA-F5-6571'] = ['Capecitabine', 'Mitomycin', 'Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab', 'Gemcitabine']
drugs_to_check_lable['TCGA-F5-6571'] =  'crc_drug_used_for_TCGA-F5-6571'

drugs_to_check['TCGA-DC-6156'] = [ 'Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab']
drugs_to_check_lable['TCGA-DC-6156'] =  'crc_drug_used_for_TCGA-DC-6156'


drugs_to_check['TCGA-AG-A02G'] = [ 'Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan', 'Bevacizumab']
drugs_to_check_lable['TCGA-AG-A02G'] =  'crc_drug_used_for_TCGA-AG-A02G'


drugs_to_check['TCGA-AG-3887'] = ['Capecitabine', 'Irinotecan', 'Bevacizumab']
drugs_to_check_lable['TCGA-AG-3887'] =  'crc_drug_used_for_TCGA-AG-3887'


drugs_to_check['TCGA-AG-3731'] = [ 'Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan']
drugs_to_check_lable['TCGA-AG-3731'] =  'crc_drug_used_for_TCGA-AG-3731'


drugs_to_check['TCGA-AG-3584'] = [ 'Fluorouracil', 'Leucovorin', 'Oxaliplatin', 'Irinotecan']
drugs_to_check_lable['TCGA-AG-3584'] =  'crc_drug_used_for_TCGA-AG-3584'


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
# patient TCGA-AZ-6600-bis  con il nuovo algoritmo


path_cp['TCGA-AZ-6600-bis'] = '../../New_R_code_2025_04_15//DrugMerge_SP_Final/3-Drugs/TCGA-AZ-6600_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp'


######
# core and peel

#####
# GEO

modulation_data[('TCGA-AZ-6600-bis', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6600-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')

# CMAP
modulation_data[('TCGA-AZ-6600-bis', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6600-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')



#####################
# patient TCGA-A6-2671-bis  con il nuovo algoritmo


path_cp['TCGA-A6-2671-bis'] = '../../New_R_code_2025_04_15//DrugMerge_SP_Final/3-Drugs/TCGA-A6-2671_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp'


######
# core and peel

#####
# GEO

modulation_data[('TCGA-A6-2671-bis', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-A6-2671-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')

# CMAP
modulation_data[('TCGA-A6-2671-bis', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-A6-2671-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')



#####################
# patient TCGA-AA-3517-bis  con il nuovo algoritmo


path_cp['TCGA-AA-3517-bis'] = '../../New_R_code_2025_04_15//DrugMerge_SP_Final/3-Drugs/TCGA-AA-3517_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp'


######
# core and peel

#####
# GEO

modulation_data[('TCGA-AA-3517-bis', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-AA-3517-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')

# CMAP
modulation_data[('TCGA-AA-3517-bis', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-AA-3517-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')


#####################
# patient TCGA-AZ-6605-bis  con il nuovo algoritmo


path_cp['TCGA-AZ-6605-bis'] = '../../New_R_code_2025_04_15//DrugMerge_SP_Final/3-Drugs/TCGA-AZ-6605_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp'


######
# core and peel

#####
# GEO

modulation_data[('TCGA-AZ-6605-bis', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6605-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')
# CMAP
modulation_data[('TCGA-AZ-6605-bis', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-AZ-6605-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')



#####################
# patient TCGA-A6-5667-bis  con il nuovo algoritmo


path_cp['TCGA-A6-5667-bis'] = '../../New_R_code_2025_04_15//DrugMerge_SP_Final/3-Drugs/TCGA-A6-5667_colorectal_cancer/Drug-modulation-data/colorectal-cancer-cp'


######
# core and peel

#####
# GEO

modulation_data[('TCGA-A6-5667-bis', 'cpd07j08', 'GEO', 'nofilter', 'crc')] = (path_cp['TCGA-A6-5667-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_GEO_noFiltering.csv')

# CMAP
modulation_data[('TCGA-A6-5667-bis', 'cpd07j08', 'CMAP', 'nofilter', 'crc')] = (path_cp['TCGA-A6-5667-bis'] , \
                                                          'Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8_Colorectal_cancer_CMAP_noFiltering.csv')



