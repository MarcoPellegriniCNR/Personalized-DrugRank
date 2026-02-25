#a Command to do tdrug-unctional term enrichement analysis
# calls the functions from the code for drugs-only and terms-only analyis

from generic_util import printmap, printlist, minuslog10, domain_restriction

from new_loader_routines_modulation import main_function

from analysis_functional_enrichments  import main_function_functional_analysis, score

import argparse

#verbose = 'verbose'
verbose = 'no-verbose'

def main_function_2():
    parser = argparse.ArgumentParser(description='Command line to analyze drug directions in patient data from enrichR.')

# parameters for the term enrichment  analysis

#obbligatorio
    parser.add_argument('disease_id',
                        action='store',
                        type = str,
                        choices=["CRC", "BRCA", 'STAD'],
                        help='Defines the disease ID. Now only cancer types CRC, BRCA or STAD')

    parser.add_argument('drug_id',
                        action='store',
                        type = str,
                        choices=['fluorouracil', 'irinotecan', 'docetaxel', 'paclitaxel', 'doxorub_hyd' , 'cyclophosphamide'],
                        help='Defines the key drug ID for sub-population identification')

    
    parser.add_argument('patient_id',
                        action='store',
                        type = str,
                        help='Defines the patient id')

#obbligatorio
    parser.add_argument('database',
                        action='store',
                        type = str,
                        choices=['GO_Biological_Process_2025', 'Reactome_Pathways_2024',  'GO_Molecular_Function_2025', 'KEGG_2021_Human'],
                        help='Defines the annotation database to be used. We use the enrichR terminology')
    

#obbligatorio
    parser.add_argument('p_value_limit',
                        action='store',
                        type = float,
                        choices=[0.05, 0.01, 0.001 , 0.0001],
                        help='Defines the pvalue filtering threshold for individual records from enrichR runs')


# paramters for the drug enrichement analysis

#obbligatorio
    parser.add_argument('drugdatatype',
                        action='store',
                        type = str,
                        choices=["GEO", "CMAP", 'L1000', 'GEO_43_fluorouracil',
                                 'GEO_brca_fluorouracil', 'GEO_brca_doxorubicin',
                                 'GEO_brca_paclitaxel', 'GEO_brca_docetaxel', 'GEO_stad_fluorouracil'],
                        help='Defines the drug database')
    

#obbligatorio
    parser.add_argument('drug_p_value_limit',
                        action='store',
                        type = float,
                        choices=[0.05, 0.01, 0.001 , 0.0001],
                        help='Defines the pvalue filering threshold for individual records')


    args = parser.parse_args()

    print(('disease_id ', args.disease_id))
    print(('key drug_id ', args.drug_id ))
    print(('patient_id ', args.patient_id))
    print(("database ", args.database))
    print(('p_value_limit ', args.p_value_limit))
    print(("drugdatatype ", args.drugdatatype))
    print(('drug_p_value_limit ', args.drug_p_value_limit))    

   
    module_Terms_map, terms_score_map_items_sorted_accumulator = main_function_functional_analysis(args.disease_id,  args.drug_id, args.patient_id, args.database, args.p_value_limit)
    module_Drugs_map = main_function(args.patient_id, args.drugdatatype, args.drug_p_value_limit)

    weights = ['log10', 'c_score', 'deltadeg', 'pval_x_deltadeg', 'pval_x_drug_hits',  'pval_x_deltadeg_x_drug_hits',  'odds_ratio_log', 'c_score_x_deltadeg' ]
#    weight = 'pval_x_deltadeg'

    id_to_terms_map = {}
    id_to_drugs_map = {}
    for x in module_Terms_map.keys():
        if x[0] not in id_to_terms_map.keys() :
            id_to_terms_map[x[0]] = [x[1],]
        else :
            id_to_terms_map[x[0]].append(x[1])
    for x in module_Drugs_map.keys():
        if x[0] not in id_to_drugs_map.keys() :
            id_to_drugs_map[x[0]] = [x[1],]
        else :
            id_to_drugs_map[x[0]].append(x[1])

    for weight in weights :
        drug_term_to_values_map = {}
        drug_term_to_values_counter_map = {}

        for module in id_to_drugs_map.keys() :
            if module in id_to_terms_map.keys() :
                list_drugs = id_to_drugs_map[module]
                list_terms = id_to_terms_map[module]
                for drug, term in zip(list_drugs, list_terms) :
                    newscore = score(module_Terms_map[(module, term )], 'deg_in_module', weight ) * score(module_Drugs_map[(module, drug)], 'deg_in_module', weight ) 
                    if (drug, term) not  in drug_term_to_values_map.keys() :
                        drug_term_to_values_map[(drug, term)] = newscore
                        drug_term_to_values_counter_map[(drug, term)] = 1
                    else :
                        drug_term_to_values_map[(drug, term)] += newscore
                        drug_term_to_values_counter_map[(drug, term)] += 1
                    if verbose == 'verbose' : 
                        print('Processing pair {}  - {}, ----- value {}'.format(drug, term, newscore))
                        
        drug_term_to_values_items = drug_term_to_values_map.items()
        drug_term_to_values_items_2 = []
        for y in drug_term_to_values_items :
            drug_term_to_values_items_2.append((y[0], y[1], drug_term_to_values_counter_map[y[0]]))
            
        drug_term_to_values_items_sorted = sorted(drug_term_to_values_items_2, key=lambda tup: tup[1], reverse=True)
        printlist(drug_term_to_values_items_sorted, 'drug_term_to_values_items_sorted', 50)


        for selected_drug in ['fluorouracil', 'irinotecan', 'cyclophosphamide', 'doxorubicin', 'decitabine']:
            print('*****************************************')
            print('Drug specific ranking for {} with weight {}'.format(selected_drug, weight))
            selected_sorted_list = [ x for x in drug_term_to_values_items_sorted if x[0][0] == selected_drug ]
            printlist(selected_sorted_list, selected_drug, 20)
            evaluate_hits(selected_sorted_list, terms_score_map_items_sorted_accumulator[weight], 20, selected_drug, weight)
        


#########################################################
# function to compute precision at 20, and RHR for the selcteed ranking of functions.
from numpy.linalg import norm

def my_L1(in_v) :
    if len(in_v) == 0 :
        return 0.0
    else :
        return sum([abs(x) for x in in_v])

import math
def my_L2(in_v) :
    if len(in_v) == 0 :
        return 0.0
    else :
        return math.sqrt(sum([x*x for x in in_v]))
    
   


def evaluate_hits(in_drug_terms_list, in_terms_list , precision_limit, selected_drug, weight):

    precision_limit2 = min(precision_limit, len(in_terms_list), len(in_drug_terms_list))

    target_set = set([x[0] for x in in_terms_list[: precision_limit2]])
    term_modules_map = {x[0] : x[2] for x in in_terms_list[: precision_limit2] }
    num_hits = 0
    total_modules = 0
    total_score = 0
    set_hits = set()
    count = 1
    rhr = 0.0
    module_vector = []
    for y in in_drug_terms_list :
        term = y[0][1]
        modules = y[2]
        score_val = y[1]
        if term in  target_set :
            if count <=  precision_limit2 :
                num_hits +=1
                set_hits.add(term)
                total_modules += modules
                total_score += score_val
                module_ratio = float(modules)/ float(term_modules_map[term])
                module_vector.append(module_ratio)
            rhr += 1/float(count)
        count +=1

    set_miss = target_set.difference(set_hits)
    l2_norm = norm(module_vector)
    l1_norm = norm(module_vector, 1)
    if l1_norm != my_L1(module_vector) :
        print('L1 gap : ', module_vector, my_L1(module_vector))
    if l2_norm != my_L2(module_vector) :
        print('L2 gap : ', module_vector, my_L2(module_vector))            
    print ('@@@@@@ For drug {}, weight {} Precision at {} is {}, rhr {}, modules {}, total score {}, L2 {} , L1 {}'.format(selected_drug, weight, precision_limit2, \
                                                                                                            num_hits, rhr, total_modules, total_score, l2_norm , l1_norm))
    printlist(list(set_hits), 'Set of hits functions ')
    printlist(list(set_miss), 'Set of missed functions ')  
            
    

    

###############################################
# main

if __name__ == "__main__":
    main_function_2()

