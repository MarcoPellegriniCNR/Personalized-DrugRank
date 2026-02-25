# -*- coding: cp1252 -*-
# 2025-02-27
# This file holds the code forevaluating the direcionality score
# We introduce filtering routines based on pvaluees
# We introduce variants due to the stricure of th record in filter/nofiltering files

import functools

################################################
# main mapping parameter to location

from new_mapping_modulation_data import modulation_data, drugs_to_check, drugs_to_check_lable

################################################
# laoding form fiel routines

from new_loader_routines import load_data

#############################################
# import the directionality score  routine

from new_directionality_analysis import direction_analysis

#####################################
# importing TTD data

from random_drug_scoring_model import  filterbycolorectcandata, filterbyprostcandata, filter_by_brca_data, compute_zscores_pval, filter_two_lists

from new_pval_computations import new_compute_zscores_pval

#############################
# functions for rank merging

from testingrra import multimergeranking

################
# functions for name analysis

from new_drug_strings_parsing import extract_drug_name_direction

#####################
# generic function

from generic_util import printlist, increment_count, printmap

####################
# extract list of drugs from file.

#@functools.lru_cache(maxsize=None)
def extract_drug_names_from_file(conf, verbose = 'noverbose'):

    separator = ','

    if not conf in modulation_data.keys() :
        print('Wrong configuration in extract_drug_names_from_file : {}'.format(conf))
        assert False
    directory = modulation_data[conf][0]
    filename = modulation_data[conf][1]
    drug_format = conf[2]

    result_accumulator, _ = load_data (directory, filename, separator)

    unique_drugs = {}
    for x in result_accumulator :
        result =  extract_drug_name_direction(x[1], drug_format)
    #   print(result)
        increment_count(unique_drugs, result[0], 1)

    drug_names_in_matrix = [y for y in  unique_drugs.keys()]
    if verbose == 'verbose' :
        printmap(unique_drugs, 'unique_drugs')
    return  drug_names_in_matrix


    
    



###########################
# main function

# upload the TTD data

def trial_drugs_TTD_data(conf, initial_list, fdaonly) :

    drugnames = [x[0] for x in initial_list]

    if conf[4] == 'crc' :
        trial_drug_list = filterbycolorectcandata(drugnames, conf)
    elif conf[4] == 'prc' :
        trial_drug_list = filterbyprostcandata(drugnames, conf)
    elif conf[4] == 'brca' :
        trial_drug_list = filter_by_brca_data(drugnames, conf)
    else:
        print('Wrong disease code ', conf[4])
        assert False

    return trial_drug_list

####

def new_ranking_evaluation_via_TTD_data(conf, initial_list, fdaonly, weight) :

    drugnames = [x[0] for x in initial_list]

    if conf[4] == 'crc' :
        trial_drug_list = filterbycolorectcandata(drugnames, conf)
    elif conf[4] == 'prc' :
        trial_drug_list = filterbyprostcandata(drugnames, conf)
    elif conf[4] == 'brca' :
        trial_drug_list = filter_by_brca_data(drugnames, conf)

    else:
        print('Wrong disease code ', conf[4])
        assert False

    compute_performance_simplified(initial_list, trial_drug_list, conf, fdaonly, weight)

    


###########################
# auxilairy functions


def compute_performance_simplified(initiallist, trialdruglist, conf, fdafilterflag, weight):

    print_limit = 20
    print("***********************************************")
    print("Listing with individual drug scores for {} drugs at print limit {}".format(len(initiallist), print_limit))
    suminvrserank = 0
    precisionat20 = 0
    rankedlist = []
    accum = [0,0,0,0]
    hits_size = 0
    for x in range(len(initiallist)) :
        drugname = initiallist[x][0] # drug name in first postion in the record                    
        if drugname in trialdruglist :
            hits_size +=1
            if x < print_limit :
                print (x, drugname, 1.0/float(x+1))
            rankedlist.append(( x, drugname, 1.0/float(x+1)))
            suminvrserank += 1.0/float(x+1)
            if x <= 19 :
                precisionat20 += 1
        else :
            if x < print_limit :
                print (x,  drugname, "no trial/TTD_Data ")
            rankedlist.append((x,  drugname, "no trial "))

    print ("@@ suminvrserank (RHR), prec-at-20 : ", conf,  weight, precisionat20, suminvrserank)
    include_old_pval_code = False
    if include_old_pval_code :
        ((z_prec,p_prec),(z_rhr,p_rhr)) = compute_zscores_pval(conf[2], conf[4], fdafilterflag , precisionat20, suminvrserank)
        print( '@@ old P_value_precision : {}, p_value_rhr : {}'.format(p_prec, p_rhr) )
    ((z_prec2,p_prec2),(z_rhr2,p_rhr2)) = new_compute_zscores_pval(len(initiallist), hits_size , precisionat20, suminvrserank)
    print( '@@ new P_value_precision : {}, p_value_rhr : {}'.format(p_prec2, p_rhr2) )
    
    return (conf , suminvrserank, precisionat20, accum, rankedlist)


##################################################################
# check fro specific drugs in a lsting

def check_specific_drugs(initial_rank_list, target_list, drug_format, trial_drug_list, weight):

    if not(set(target_list).issubset(set( trial_drug_list))) :
        print('Warning: Mismatches in names: ', set(target_list), set( trial_drug_list))

##    if  'Doxorubicin Hydrochloride' in target_list or 'Doxorubicin_Hydrochloride' in target_list :
##        target_list.append('Doxorubicin')
##
##    target_list2 = [x.lower() for x in target_list]
##    trial_drug_list2 = [x.lower() for x in trial_drug_list]

    print('Target list  :'  , target_list)
    if 'doxorubicin_hydrochloride' in target_list :
        target_list.remove('doxorubicin_hydrochloride')
        target_list.append('doxorubicin')
    

    print('Iniltal list size {}, trial drugs size {}, target list size {}'.format(len(initial_rank_list), len(trial_drug_list), len(target_list)))
          
    if len(initial_rank_list) == 0 or  len(target_list) == 0 :
        return []
    rel_rank = 1
    for x in range(len(initial_rank_list)):
        drugname = initial_rank_list[x][0] # drug name in first postion in the record                    
        if drugname in target_list :
            print('@@ For {} Drug {} weight {}, in position {} on total of {}, relative rank is {} of {}'.format(drug_format, drugname,  weight, x, len(initial_rank_list), rel_rank, len(trial_drug_list)))
        if drugname in trial_drug_list :
            rel_rank += 1

    return []
        


#####################################################
# tesing the evaluation of the directionality score  computations


if __name__ == "__main__":
#    pass
    print('Testing computation of the directionality score')

    separator = ','
    counter = 1
    counter_test_limit = 100
#    with_merge_step = True
    with_merge_step = False


#    patient_list = ['pilot', ]
 #   patient_list = ['TCGA-AZ-6600', ]
#    patient_list = ['TCGA-HC-7742',]
#    patient_list = ['TCGA-A6-2671',]
#    patient_list = ['TCGA-AA-3517',]
#    patient_list = ['TCGA-AZ-6605',]
#    patient_list = ['TCGA-A6-5667',]

#    patient_list = ['TCGA-AZ-6600-bis', ]
#    patient_list = ['TCGA-A6-2671-bis',]
#    patient_list = ['TCGA-AA-3517-bis',]
#    patient_list = ['TCGA-AZ-6605-bis',]
    patient_list = ['TCGA-A6-5667-bis',]

#    alg_list = ['cpd07j08', 'MD' ]
    alg_list = ['cpd07j08', ]

# list of datrabases

    db_list = ['CMAP',]

    weights = ["size_x_deltadeg", "deltadeg", "pval_x_deltadeg", "pval_x_log_deltadeg", 'pval_x_drug_hits', 'pval_x_deltadeg_x_drug_hits']

    print('Paremetrs: counter_test_limit {}, with_merge_step {}, patient_list {}, alg_list {}, db_list {}, weights {} '.format(counter_test_limit, with_merge_step, patient_list, alg_list, db_list, weights))
    collect_results_map = {}
    for conf in modulation_data.keys() :
        if conf[0] in patient_list  and  conf[1] in alg_list and conf[2] in db_list  :
            print('--------------------------- Begin --------------------------------')
            print(conf, modulation_data[conf])
            directory = modulation_data[conf][0]
            filename = modulation_data[conf][1]
            drug_format = conf[2]
            filtertype = conf[3]
            positiveonly = True #per ora  usato
            fdaonly = True  #per ora  usato
#            fdaonly = False
            pvalue_filter = True
            p_value_limit = 0.01
#            positiveonly = False #per ora  usato
#            fdaonly = False  #per ora  usato

#            weights = "deltadeg"
#            weights = "size_x_deltadeg"
#            weights = "pval_x_deltadeg"
#            weights = "pval_x_log_deltadeg"
#            weights = 'pval_x_drug_hits'
#            weights = 'pval_x_deltadeg_x_drug_hits'

            print('Paremetrs: positiveonly {}, fdaonly {}, weights {}, pvalue_filter {},  p_value_limit {} '.format(positiveonly, fdaonly, weights, pvalue_filter, p_value_limit ))
            

    #        if not filename.endswith('noFiltering.csv') and not filename.endswith('noF.csv'):
            if conf[3] in ['filter', 'nofilter']  :
                for weight in weights :
                    counter +=1
                    result_accumulator = direction_analysis(directory, filename, separator, drug_format, positiveonly,  fdaonly ,pvalue_filter, p_value_limit, weight)
                    new_ranking_evaluation_via_TTD_data(conf, result_accumulator, fdaonly, weight)

        # check if therapy drigs dta is available form apient and in teh matrix
                    if conf[0] not in drugs_to_check.keys() :
                        print('No drugs_to_check for {}'.format(conf[0]))
                    elif drugs_to_check[conf[0]] == [] :
                         print('No drugs_to_check for {}'.format(conf[0]))
                    else :
                        drug_names_in_matrix = extract_drug_names_from_file(conf, verbose = 'noverbose')
                        target_drugs = filter_two_lists(drugs_to_check_lable[conf[0]], drugs_to_check[conf[0]], filename, drug_names_in_matrix )
                        trial_drugs = trial_drugs_TTD_data(conf, result_accumulator, fdaonly)
                        check_specific_drugs(result_accumulator, target_drugs, drug_format, trial_drugs, weight)
                    
                    collect_results_map[str(conf) + ' ' + weight] = result_accumulator
            else :
                print('No code for this file format yet', filename)
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break

#  merge multiple configurations
#  collect GA drug data lables.

    if with_merge_step : 
        drug_ge_dataset_lables = [conf[2] for conf in collect_results_map.keys()]
        drug_ge_dataset_lables_set = set(drug_ge_dataset_lables)
        for pat in patient_list :
            for  drug_ge in drug_ge_dataset_lables_set :
                list_confs_to_merge = [conf for conf in collect_results_map.keys() if conf[0] == pat and conf[2] == drug_ge]
                nconfs=len(list_confs_to_merge)
                print('--------------------------- Begin --------------------------------')
                print('Processing list of {} confs {}'.format(nconfs, list_confs_to_merge))
                if len(list_confs_to_merge) > 1 :
                    list_of_lists = [dict(collect_results_map[x]) for x in list_confs_to_merge]
                    (totallist, mergedvector, outmatrix2) = multimergeranking(list_of_lists)
                    inputoutputmatrix2 = [outmatrix2[x] + [mergedvector[x], totallist[x]] for x in range(len(totallist))]
                    lexsorting12 = sorted(inputoutputmatrix2, reverse = True, key=lambda tup: sum(tup[:nconfs]) )
                    lexsorting22 = sorted(lexsorting12, key=lambda tup: tup[nconfs])
                    printlist( lexsorting22  , "inputoutputmatrix", 20)
                    result_accumulator = [(x[-1], x[-2]) for x in lexsorting22]
                    new_ranking_evaluation_via_TTD_data(list_confs_to_merge[0], result_accumulator, fdaonly)
                print('--------------------------- End --------------------------------')    

