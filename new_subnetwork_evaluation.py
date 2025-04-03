# -*- coding: cp1252 -*-
# 2025-03-19
# This file holds the code for computing the ranking of the active subnetwoerks  of data uploaded
# from file


import os

################################################
# main mapping parameter to location

#from new_mapping_modulation_data import modulation_data
from new_mapping_subnetwork_data import  subnetwork_data
from new_mapping_modulation_data import drugs_to_check, drugs_to_check_lable

################################################
# generic utilities

from generic_util import printmap, printlist, increment_count, domain_restriction

##############################################
# import file handling utils

from file_handling_util import testexistencedirectory, testexistencefile


####################################
# loading function and aprsinf functions

from new_loader_subnetwork_routines import load_data_subnetwork

from new_drug_strings_parsing import extract_drug_name_direction



##################################
# loaign name filtering routines

from random_drug_scoring_model import filterbyfdaapproveddrugs, filter_two_lists

from new_directionality_analysis  import test_pval

######################################################
# main function

from new_evaluation_analysis import new_ranking_evaluation_via_TTD_data, check_specific_drugs
from new_subnetwork_analysis  import subnetwork_analysis

##def subnetwork_analysis(directory, filename, separator, drug_format,  fdaonly) :
##
##    assert drug_format  in ['GEO', 'DM', 'CMAP', 'L1000']
##
##    dataarray , clean_header = load_data_subnetwork(directory, filename, separator)   
##    print ("Number of records ", len(dataarray))
##    print ('First raw : ' , dataarray[0])
###    assert len(dataarray[0]) == 9 #code is only for the foramt at 9 fields
##    clean_header_map = {clean_header[x] : x for x in range(len(clean_header))}
##    printmap(clean_header_map, 'clean_header_map') # maps field names to field indices spcific for the filter files.
##
##    # p-value filter probably not so relevant if ddata is already filtered.
##    pval_filtered_dataarray = [x for x in dataarray if test_pval(x[clean_header_map['Adjusted.P.value']], 0.05)] # use the raw p_value not the adjusted one 
##    print ("Number of records after pval filtering ", len(pval_filtered_dataarray))
##
##    # extraction fo the relevant fields.
##
##    drug_score_map = {}
##
##    for x in pval_filtered_dataarray :
##        (drug_name, drug_direction) = extract_drug_name_direction(x[clean_header_map['Term']], drug_format)
##        p_val = float(x[clean_header_map['Adjusted.P.value']])
##        if drug_name  not in drug_score_map.keys() :
##            drug_score_map[drug_name] = p_val
##        else  :
##            drug_score_map[drug_name] = min(drug_score_map[drug_name], p_val) # we use the smallest p_value
##            
##
##    drugscore_items_sorted = sorted(drug_score_map.items(), key=lambda tup: tup[1], reverse=False)
##    printlist(drugscore_items_sorted, 'drugscore_items_sorted', 20)
##
### fda filter
##    if fdaonly == True :
##        list_of_drug_names = [x[0] for x in drug_score_map.items()]
##        fda_filtered_list_of_drug_names = filterbyfdaapproveddrugs(list_of_drug_names, drug_format)
##        drugscore_fda_filtered = domain_restriction(fda_filtered_list_of_drug_names, drug_score_map)
##        drugscore_items = drugscore_fda_filtered.items()
##        drugscore_items_sorted = sorted(drugscore_items, key=lambda tup: tup[1], reverse=False)
##        printlist(drugscore_items_sorted, 'drugscore_items_sorted after fda filtering', 20)
##
##
##    return drugscore_items_sorted

#####################################################
# tesing the pvalue  score  computations


if __name__ == "__main__":
#    pass
    print('Testing computation of the pvalue score')

    separator = ','
    counter = 1
    counter_test_limit = 100
    
#    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
#    patient_list = ['TCGA-HC-7742']
#    patient_list = ['TCGA-A6-2671']
#    patient_list = ['TCGA-AA-3517']
#    patient_list = ['TCGA-AZ-6605',]
    patient_list = ['TCGA-A6-5667']
    alg_list = ['cpd07j08', 'MD' , 'KPM', 'Degas', 'clex3800' , 'clex3000']
#    data_list = [ 'GEO', 'DM']
    data_list = [  'CMAPup' , 'L1000up' ,  'CMAPdown' , 'L1000down', 'CMAPupdown' , 'L1000updown',  'GEO', 'DM']

    collect_results_map = {}
    for conf in subnetwork_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list  and conf[2] in data_list :
            print('--------------------------- Begin --------------------------------')
##            print(conf, subnetwork_data[conf])
##            directory = subnetwork_data[conf][0]
##            filename = subnetwork_data[conf][1]
            drug_format = conf[2]
            positiveonly = True #per ora  non usato
            fdaonly = True  #per ora  usato
#        weights = "deltadeg"
#       weights = "size_x_deltadeg"
            
            counter +=1
            result_accumulator = subnetwork_analysis(conf, separator, drug_format,   fdaonly)
            new_ranking_evaluation_via_TTD_data(conf, result_accumulator, fdaonly)
            drug_names_list = [x[0]   for x in result_accumulator ]
            collect_results_map[conf] = result_accumulator

    # check if therapy drigs dta is available form apient and in teh matrix
            if conf[0] not in drugs_to_check.keys() :
                print('No drugs_to_check for {}'.format(conf[0]))
            elif drugs_to_check[conf[0]] == [] :
                print('No drugs_to_check for {}'.format(conf[0]))
            else :
                target_drugs = filter_two_lists(drugs_to_check_lable[conf[0]], drugs_to_check[conf[0]], str(conf), drug_names_list )  
                check_specific_drugs(result_accumulator, target_drugs, drug_format)
            
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break

# here goes the merging code for the individual ranings -  yet to be developed.




