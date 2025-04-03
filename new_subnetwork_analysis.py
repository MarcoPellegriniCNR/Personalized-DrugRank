# -*- coding: cp1252 -*-
# 2025-03-19
# This file holds the code for computing the ranking of the active subnetwoerks  of data uploaded
# from file


import os

################################################
# main mapping parameter to location

#from new_mapping_modulation_data import modulation_data
from new_mapping_subnetwork_data import  subnetwork_data

################################################
# generic utilities

from generic_util import printmap, printlist, increment_count, domain_restriction

##############################################
# import file handling utils

from file_handling_util import testexistencedirectory, testexistencefile


####################################
# loading function and aprsinf functions

from new_loader_subnetwork_routines import load_data_subnetwork, load_data_subnetwork_from_conf

from new_drug_strings_parsing import extract_drug_name_direction

from new_drug_string_parsing_subnetwork import extract_drug_name_only



##################################
# loaign name filtering routines

from random_drug_scoring_model import filterbyfdaapproveddrugs

from new_directionality_analysis  import test_pval

######################################################
# main function

def subnetwork_analysis(conf, separator, drug_format,  fdaonly) :

    assert drug_format  in ['GEO', 'DM', 'CMAPup', 'L1000up', 'CMAPdown', 'L1000down', 'CMAPupdown', 'L1000updown']

    dataarray , clean_header = load_data_subnetwork_from_conf(conf, separator)   

    print ("Number of records ", len(dataarray))
    if len(dataarray) > 0 :
        print ('First raw : ' , dataarray[0])
    else :
        return []
        
#    assert len(dataarray[0]) == 9 #code is only for the foramt at 9 fields
    clean_header_map = {clean_header[x] : x for x in range(len(clean_header))}
    printmap(clean_header_map, 'clean_header_map') # maps field names to field indices spcific for the filter files.

    # p-value filter probably not so relevant if ddata is already filtered.
    pval_filtered_dataarray = [x for x in dataarray if test_pval(x[clean_header_map['Adjusted.P.value']], 0.05)] # use the raw p_value not the adjusted one 
    print ("Number of records after pval filtering ", len(pval_filtered_dataarray))

# add here the toxic substance filter  ?

    # extraction fo the relevant fields.

    drug_score_map = {}

    for x in pval_filtered_dataarray :
        drug_name= extract_drug_name_only(x[clean_header_map['Term']], drug_format)
        p_val = float(x[clean_header_map['Adjusted.P.value']])
        if drug_name  not in drug_score_map.keys() :
            drug_score_map[drug_name] = p_val
        else  :
            drug_score_map[drug_name] = min(drug_score_map[drug_name], p_val) # we use the smallest p_value
            

    drugscore_items_sorted = sorted(drug_score_map.items(), key=lambda tup: tup[1], reverse=False)
    printlist(drugscore_items_sorted, 'drugscore_items_sorted', 20)

# fda filter
    if fdaonly == True :
        list_of_drug_names = [x[0] for x in drug_score_map.items()]
        fda_filtered_list_of_drug_names = filterbyfdaapproveddrugs(list_of_drug_names, drug_format)
        drugscore_fda_filtered = domain_restriction(fda_filtered_list_of_drug_names, drug_score_map, 'nosubset')
        drugscore_items = drugscore_fda_filtered.items()
        drugscore_items_sorted = sorted(drugscore_items, key=lambda tup: tup[1], reverse=False)
        printlist(drugscore_items_sorted, 'drugscore_items_sorted after fda filtering', 20)


    return drugscore_items_sorted

#####################################################
# tesing the pvalue  score  computations


if __name__ == "__main__":
#    pass
    print('Testing computation of the pvalue score')

    separator = ','
    counter = 1
    counter_test_limit = 10
    
#    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
#    patient_list = ['TCGA-HC-7742']
#    patient_list = ['TCGA-A6-2671']
#    patient_list = ['TCGA-AA-3517']
    patient_list = ['TCGA-AZ-6605',]
    alg_list = ['cpd07j08', 'MD' , 'KPM', 'Degas', 'clex3800' ]
#    data_list = [ 'GEO', 'DM', 'CMAPup' ]
#    data_list = [  'CMAPup' , 'L1000up' ]
#    data_list = [  'CMAPdown' , 'L1000down' ]
    data_list = [  'CMAPupdown' , 'L1000updown' ]


    for conf in subnetwork_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list  and conf[2] in data_list :
            print('--------------------------- Begin --------------------------------')
            drug_format = conf[2]
            positiveonly = True #per ora  non usato
            fdaonly = True  #per ora  usato
#        weights = "deltadeg"
#       weights = "size_x_deltadeg"
            
            counter +=1
            result_accumulator = subnetwork_analysis(conf, separator, drug_format,   fdaonly)
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break
