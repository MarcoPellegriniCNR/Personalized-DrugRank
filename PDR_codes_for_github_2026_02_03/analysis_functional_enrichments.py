# This is the main code for the analysis of functional enrichment in  personalized disease modules for a single patient.
# It is a variation of the code in `new_loader_routines_modualtion.py'

# imports
import os
import os.path
import math

# generic utilities
from new_loader_routines import load_data
from file_handling_util import testexistencedirectory, testexistencefile
from generic_util import printmap, printlist, minuslog10, domain_restriction




###########################################
# main function called by the command line level:

def main_function_functional_analysis(disease_id, drug_id, patient_id, database , p_value_limit) :


# recover the disease modules (not personalized) based on disease id.
    cp_active_modules_file = recover_active_diseae_modules(disease_id)
    active_modules_map = load_modules_form_file_as_map('.', cp_active_modules_file)

# recover the single patients filtered DEG lists up/down  and the path to the directory with the erichR files.
    (path_down, path_up,  patient_up_genes_file, patient_down_genes_file) = recover_patient_data(disease_id, drug_id, patient_id, database)

# load DEG personal files    
    patient_up_genes = load_gene_list_as_set('.', patient_up_genes_file)
    patient_down_genes = load_gene_list_as_set('.', patient_down_genes_file)

# compute length of the personalized disease  modules
    active_modules_map_reduced_with_patient_up_genes = reduce_map_lengths(active_modules_map, patient_up_genes)
    active_modules_map_reduced_with_patient_down_genes = reduce_map_lengths(active_modules_map, patient_down_genes)
    printmap(active_modules_map_reduced_with_patient_up_genes, 'Some Personalized up  genes modules :  ', 10)
    printmap(active_modules_map_reduced_with_patient_down_genes, 'Some Personalize down geens modules : ', 10)

# set up of the main loop: make a list of the enrichR files to upload 

    files_up = sorted([f for f in os.listdir(path_up) if os.path.isfile(os.path.join(path_up, f)) and f.startswith('enrichR_module')])
    printlist(files_up, 'Some of the files up', 10)
    pairs_up = [(path_up, x) for x in files_up ]

    files_down = sorted([f for f in os.listdir(path_down) if os.path.isfile(os.path.join(path_down, f))  and f.startswith('enrichR_module')])
    printlist(files_down, 'Some of the  files down', 10)
    pairs_down = [(path_down, x) for x in files_down ]
    
    all_pairs =  pairs_up + pairs_down

# set up the main loop 
#    verbose = 'verbose'
    verbose = 'no_verbose'
    weights = ['log10', 'c_score', 'deltadeg', 'pval_x_deltadeg', 'pval_x_drug_hits',  'pval_x_deltadeg_x_drug_hits',  'odds_ratio_log', 'c_score_x_deltadeg' ]
    counter = 0
    pvalue_filter = True
    Terms_names = set()  # set of unique generic terms
    Terms_module_map = {}
    counter_limit = 1000000 #an escape value for debugging-testing
#    counter_limit = 10  #an escape value for debugging-testing

    # main loop
    print("@@ Starting p_value : {}".format(p_value_limit))
    for x in all_pairs:
        # exacting module id from filename
        module_id = extract_module_number(x[1])
        counter +=1
        print('Counter {} , module n: {} '.format(counter, module_id ))

        # here the definition si just the sum of teh up and down degrees as it is term-neural.
        delta_deg = active_modules_map_reduced_with_patient_up_genes[module_id] + active_modules_map_reduced_with_patient_down_genes[module_id]

        # load the data from erichR file
        acc ,  header = load_data(x[0], x[1], ',', 'noverbose')

        # apply p-value filter
        if pvalue_filter :
            pval_filtered_dataarray = [z for z in acc  if test_pval(z[record_map['P.value']], p_value_limit)]
            acc = pval_filtered_dataarray

        for row in acc :

            (term , p_val, minuslog, c_score,  DE_drug_num, odds_ratio_log) = extract_from_enrichR_record(row)

            if term == None :
                continue
            new_record =  (p_val, minuslog, c_score, delta_deg, DE_drug_num, odds_ratio_log)
            if verbose == 'verbose' :
                print('Term {}, record : {}'.format(term, new_record))
            Terms_names.add(term) # update set of terms

            # store the data in the module-term map, possible that the same term appears more then once in the fle for a module
            if (module_id, term) not in   Terms_module_map.keys() :        
                Terms_module_map[(module_id, term)] = [new_record]
            else :
                Terms_module_map[(module_id, term)].append(new_record)

        # early escape
        if counter >= counter_limit :
            break

# print statistics of the module-term map
    print('Counter of file opened {} over {} files listed '.format(counter, len(all_pairs)))
    print ("Number of  unique terms : " , len(Terms_names))


#  reorganize counters by enriched terms and make summations of the scores and compute averages
# overwites
    
    terms_score_map_items_sorted_accumulator = {}
    for weight in weights :
        terms_score_map = {}
        terms_modules = {}
        terms_module_counter ={}
        for term in Terms_names :
            terms_score_map[term] = 0
            terms_modules[term] = []
            terms_module_counter[term] = 0 


        for x in Terms_module_map.keys() :
            term = x[1]

            terms_score_map[term] +=  score(Terms_module_map[x], 'deg_in_module', weight)
            terms_modules[term].append(x[0]) # record the module id
            terms_module_counter[term] += 1 

        terms_score_map_items = terms_score_map.items()
        terms_score_map_items_sorted = sorted(terms_score_map_items, key=lambda tup: tup[1], reverse=True)

        # add the number of modules as additional info of the record
        terms_score_map_items_sorted_2 = []
        for y in  terms_score_map_items_sorted :
            terms_score_map_items_sorted_2.append( (y[0], y[1], terms_module_counter[y[0]] ) )
            
        print("********************************************************************************")
        printlist(terms_score_map_items_sorted_2, 'Extended terms_score_map_item_sorted for weight : ' + weight, 30)          
        terms_score_map_items_sorted_accumulator[weight] = terms_score_map_items_sorted_2

    return Terms_module_map, terms_score_map_items_sorted_accumulator

# quality evaluation via RHR  scoring (not yet  implemented) 


# end of the main procedure..


####################################
####################################
# support functions


# recover the map of the active disease module (not personalized)
 # includ file  existence checks
    
def recover_active_diseae_modules(disease_id):

    cp_active_modules_file_CRC = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'
    cp_active_modules_file_BRCA = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/active_subnetworks/Breast_cancer/MultipleModules/BRCA_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'
    cp_active_modules_file_STAD = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/active_subnetworks/Stomach_cancer/MultipleModules/STAD_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'

    if  disease_id == 'CRC':
        testexistencefile('.', cp_active_modules_file_CRC)
        return cp_active_modules_file_CRC
    elif disease_id == 'BRCA':
        testexistencefile('.', cp_active_modules_file_BRCA)
        return cp_active_modules_file_BRCA
    elif disease_id == 'STAD':
        testexistencefile('.', cp_active_modules_file_STAD)
        return cp_active_modules_file_STAD
    else :
        print('Wrong disease ID no file for active disease modules')
        assert False



# recover the path to the list of DEG for the


def recover_patient_data(disease_id, drug_id, patient_id, database) :

    root =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/5-MultipleModules_Functional_Enrichment/'

    if (disease_id== 'CRC') and (drug_id == 'irinotecan'):

        path_down =  root +   'Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =       root +   'Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_bis/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_bis/' + patient_id +  '_4000_DOWN_gene_list.txt'


        
    elif (disease_id== 'CRC') and (drug_id == 'fluorouracil'):
        path_down =  root +   'Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =       root +   'Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'



    elif (disease_id== 'BRCA') and (drug_id == 'fluorouracil') :
        path_down =   root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =        root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'



    elif (disease_id== 'BRCA') and (drug_id =='doxorub_hyd') :
        path_down =   root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =        root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'



    elif (disease_id== 'BRCA') and (drug_id =='paclitaxel' ):
        path_down =   root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =        root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'



    elif (disease_id== 'BRCA') and (drug_id =='docetaxel') :
        path_down =   root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =        root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'


    elif (disease_id== 'BRCA') and (drug_id =='cyclophosphamide') :
        path_down =   root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_down' 
        path_up =        root + 'Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_cyclophosphamide_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_cyclophosphamide_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'


    elif (disease_id== 'STAD') and (drug_id =='fluorouracil' ):
        path_down = root + 'Stomach_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/'+ database + '_down' 
        path_up =     root + 'Stomach_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/' + database + '_up' 
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Stomach_cancer_sp2_fluorouracil/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Stomach_cancer_sp2_fluorouracil/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO' 

    else :
        print('Not valid combination of disease_id  {} and drug_id {}'.format(disease_id, drug_id))
        assert False

# testing existence of files and diectories
    if not testexistencedirectory(path_down):
        assert False
    if not testexistencedirectory(path_up):
        assert False

    # new with common background and 4000 FC filter

    if not testexistencefile('.', patient_up_genes_file):
        assert False        
    if not testexistencefile('.', patient_down_genes_file):
        assert False

    return (path_down, path_up,  patient_up_genes_file, patient_down_genes_file)



# functio to load modules form a modules file

def load_modules_form_file_as_map(directory, filename) :

    separator = ' '
    f = open(os.path.join(directory,filename))
    line = f.readline()
    linenumber = 1
    out_map = {}
    while line != '' :
        out_map[linenumber] = set(line.strip().split(separator))
        line = f.readline()
        linenumber +=1
        
    f.close()
    print('Read lines: {} from file {} '.format( linenumber -1, filename))
    printmap(out_map, 'First Modules', 5)

    return out_map

# function to laog genes form a gene list.

def load_gene_list_as_set(directory, filename) :

    separator = ' '
    f = open(os.path.join(directory,filename))
    line = f.readline()
    linenumber = 1
    out_list = []
    while line != '' :
        gene = line.strip()
        if gene != '' :
            out_list.append(gene)
            linenumber +=1
        line = f.readline()
    f.close()     
    print('Read lines: {} from {} '.format( linenumber -1, filename))
    print('First items : ', out_list[:10])

    return set(out_list)


# computes the lenght of  the reduction of  sets as codomains in a an input map
#

def reduce_map_lengths(in_map, list_as_set) :

    out_map = {}
    for x in in_map.keys() :
        out_map[x] = len(in_map[x] & list_as_set)
    return out_map

# extact module id form the filename:

def extract_module_number(in_string) :

    x = in_string.lstrip('enrichR_module')
    s = x.rstrip('.csv')
    return int(s)

# scheme of fields for enrichR records

record_map = {
'Term' : 0 ,
'Overlap' : 1,
'P.value' : 2 ,
'Adjusted.P.value': 3,
'Old.P.value': 4,
'Old.Adjusted.P.value' : 5,
'Odds.Ratio': 6,
'Combined.Score': 7,
'Genes': 8
    }

# function to test a pvalue from file against a threshold. 
def test_pval(s, pval_limit) :

    if s == None :
        return False
    if s in ['', 'NA', 'na'] :
        return False
    try:
        x = float(s)
#        print (x)
        return x <= pval_limit
    except ValueError:
        print(s, " Not a float")
        assert False

# extracing data from a single enrichR  file row (as a tuple)
# if data has to be extracted from terms the code should go here

def extract_from_enrichR_record(row) :

    term = row[record_map['Term']]
    if term in [None,  '' , 'NA']:
        term = None
        
    p_val = float(row[record_map['P.value']])
    minuslog = minuslog10(p_val)    
    c_score = float(row[record_map['Combined.Score']])
    DE_drug_num = len(row[record_map['Genes']].split(';'))
    odds_ratio_log = math.log10(float(row[record_map['Odds.Ratio']]))

    return  (term , p_val, minuslog, c_score,  DE_drug_num, odds_ratio_log)



#######################################
# sum all scores over a list or records and take the average velue 

def score(list_of_records, typescore, weights) :

    assert len(list_of_records) > 0
    sc = 0
    for x in list_of_records :
#        sc += scoreonerecord(x, typescore, weights)
        if weights == 'log10' :
            sc += x[1]  # add the minuslog10 values..
        elif weights == 'c_score' :
            sc += x[2]
        elif weights == 'deltadeg' :
            sc += x[3]
        elif weights == 'pval_x_deltadeg' :
            sc += x[1]*x[3]
        elif weights == 'pval_x_drug_hits' :
            sc  += x[1]*x[4]
        elif weights == 'pval_x_deltadeg_x_drug_hits' :
            sc += x[1]*x[4]*x[3]
        elif weights == 'odds_ratio_log' :
            sc  += x[5]
        elif  weights == 'c_score_x_deltadeg' :
            sc += x[2] *x[3]
        else :
            print('Warning weight not implemente :', weights)
            assert False

    return float(sc)/float(len(list_of_records))

