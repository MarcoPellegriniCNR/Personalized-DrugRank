# new loader routines for modulation enrichR data, third version of pipeline 2

import os
import os.path
import math

from new_mapping_modulation_data import modulation_data

drug_data_from_file = 'BRCA_all_treatments_2.csv'
#drug_data_from_file = 'dati_clinici_pz_fluorouracil_CRC.csv'
#drug_data_from_file = 'new_mapping_modulation_data'

print('Uplaoding patient treatment data form : ', drug_data_from_file)
if drug_data_from_file == 'dati_clinici_pz_fluorouracil_CRC.csv' :
    from clinical_data_analysis2 import gen_drugs_to_check_maps
    (drugs_to_check, drugs_to_check_lable, reduced_response_map) = gen_drugs_to_check_maps()
    
elif drug_data_from_file == 'BRCA_all_treatments_2.csv' :
    from clinical_data_analysis2 import BRCA_gen_drugs_to_check_maps
    (drugs_to_check, drugs_to_check_lable, reduced_response_map) = BRCA_gen_drugs_to_check_maps()
    
elif  drug_data_from_file == 'new_mapping_modulation_data' :
    from new_mapping_modulation_data import drugs_to_check, drugs_to_check_lable
    
else :
    print('Wrong patient treatmetn data file indicator : ', drug_data_from_file )
    assert False

from new_loader_routines import load_data

from file_handling_util import testexistencedirectory, testexistencefile

from generic_util import printmap, printlist, minuslog10, domain_restriction

from  new_drug_string_parsing_subnetwork import extract_drug_name_only
from new_drug_strings_parsing import extract_drug_name_direction

from random_drug_scoring_model import filterbyfdaapproveddrugs, filter_two_lists

from new_directionality_analysis import test_pval

from new_evaluation_analysis import new_ranking_evaluation_via_TTD_data, extract_drug_names_from_file, check_specific_drugs, trial_drugs_TTD_data


##patient_id = 'TCGA-A6-2671'
##path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/TCGA-A6-2671/Core&Peel-d0.7-j0.8/CMAPdown'
##path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/TCGA-A6-2671/Core&Peel-d0.7-j0.8/CMAPup'

# with threapy
#patient_id = 'TCGA-AZ-6600'
#patient_id = 'TCGA-AA-3517'
#patient_id = 'TCGA-A6-2671'

# without therapy
#patient_id = 'TCGA-AZ-6605'
#patient_id =   'TCGA-A6-5667'


################
# nuova seri di irnotecan patients

#patient_id = 'TCGA-AA-3679'

#patient_id = 'TCGA-A6-2671'

##patient_id = 'TCGA-AA-3869'

##patient_id = 'TCGA-CM-5868'

##patient_id = 'TCGA-NH-A6GA'

##patient_id = 'TCGA-AY-4070'

#patient_id = 'TCGA-AA-A01R'

#patient_id = 'TCGA-AA-A02F'



###########################


##patient_id = 'TCGA-AA-3973'

##patient_id = 'TCGA-AA-3548'

##patient_id = 'TCGA-F4-6806'

##patient_id = 'TCGA-A6-5664'

##patient_id =  'TCGA-A6-2683'

##patient_id =  'TCGA-AA-A01D'

##patient_id = 'TCGA-CM-4748'

#patient_id = 'TCGA-CM-6675'

#patient_id = 'TCGA-CM-6678'

#patient_id = 'TCGA-G4-6297'

#patient_id = 'TCGA-G4-6303'

#patient_id = 'TCGA-AA-3494'

#############################

#patient_id = 'TCGA-AZ-4682'

############################
##

##
##
##

##


################


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


def extract_module_number(in_string) :

    x = in_string.lstrip('enrichR_module')
    s = x.rstrip('.csv')
    return int(s)


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


##############################
# loading module and up/down regulated genes of the patient

##patient_up_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_' + patient_id +  '_UP_gene_list.txt'
##patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_' + patient_id +  '_DOWN_gene_list.txt'

# new with common background and 1.5 FC filter
##patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_1.5_UP_gene_list.txt'
##patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_-1.5_DOWN_gene_list.txt'

# new with common background and 1.0 FC filter
##patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_1.0_UP_gene_list.txt'
##patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_-1.0_DOWN_gene_list.txt'

# new with common background and 4000 FC filter
#patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_UP_gene_list.txt'
#patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_DOWN_gene_list.txt'


cp_active_modules_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/active_subnetworks/Colorectal_cancer/MultipleModules/Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'
cp_active_modules_file_brca = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/active_subnetworks/Breast_cancer/MultipleModules/BRCA_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt'

##testexistencefile('.', patient_up_genes_file)
##testexistencefile('.', patient_down_genes_file)
testexistencefile('.', cp_active_modules_file)

def load_modules(directory, filename) :

    separator = ' '
    f = open(os.path.join(directory,filename))
    line = f.readline()
    linenumber = 1
    out_map = {}
    while line != '' :
        out_map[linenumber] = set(line.strip().split(separator))
        line = f.readline()
        linenumber +=1

    print('read lines: {} from {} '.format( linenumber -1, filename))
    f.close()
    return out_map

#printmap( load_modules('.', cp_active_modules_file), 'mappa moduli', 5)



def load_gene_list(directory, filename) :

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
        

    print('read lines: {} from {} '.format( linenumber -1, filename))
    f.close()
    return set(out_list)

##printlist(load_gene_list('.', patient_up_genes_file), 'liste geni up del pazinte', 10 )
##printlist(load_gene_list('.', patient_down_genes_file), 'liste geni down del pazinte', 10 )

def reduce_map_lengths(in_map, lista) :

    out_map = {}
    for x in in_map.keys() :
        out_map[x] = len(in_map[x] & lista)
    return out_map



##############################

##active_modules_map = load_modules('.', cp_active_modules_file)
##patient_up_genes = load_gene_list('.', patient_up_genes_file)
##patient_down_genes = load_gene_list('.', patient_down_genes_file)
##
##
##active_modules_map_reduced_with_patient_up_genes = reduce_map_lengths(active_modules_map, patient_up_genes)
##active_modules_map_reduced_with_patient_down_genes = reduce_map_lengths(active_modules_map, patient_down_genes)
##
##
##printmap(active_modules_map_reduced_with_patient_up_genes, 'up', 10)
##printmap(active_modules_map_reduced_with_patient_down_genes, 'down', 10)

#################################


def main_function(patient_id, drugdatatype , p_value_limit) :



    if drugdatatype == 'GEO' :
        # for GEO data
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO'

    elif drugdatatype == 'CMAP' :
        # for CMAP data
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/CMAPdown'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/CMAPup'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'CMAP'

    elif drugdatatype == 'L1000' :
        #for L1000 data 
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/L1000down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/L1000up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'L1000'
        
    elif drugdatatype == 'GEO_43_fluorouracil' :
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer_fluorouracil/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Colorectal_cancer_fluorouracil/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO'

    elif drugdatatype == 'GEO_brca_fluorouracil' :
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO'

    elif drugdatatype == 'GEO_brca_doxorubicin' :
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO'

    elif drugdatatype =='GEO_brca_paclitaxel' :
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO'

    elif drugdatatype == 'GEO_brca_docetaxel' :
        path_down = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014down'
        path_up =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/4-MultipleModules_output/Breast_cancer/' + patient_id + '/Core&Peel-d0.7-j0.8/Drug_Perturbations_from_GEO_2014up'
        patient_up_genes_file =      '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients/' + patient_id +  '_4000_UP_gene_list.txt'
        patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients/' + patient_id +  '_4000_DOWN_gene_list.txt'
        drugdatatype1 = 'GEO'
      

    else :
        print('Not valid drugdatatype {}'.format(drugdatatype))
        assert False

    if not testexistencedirectory(path_down):
        assert False
    if not testexistencedirectory(path_up):
        assert False

    # new with common background and 4000 FC filter

    if not testexistencefile('.', patient_up_genes_file):
        assert False
        
    if not testexistencefile('.', patient_down_genes_file):
        assert False
        
    if drugdatatype in ['GEO_brca_fluorouracil', 'GEO_brca_doxorubicin', 'GEO_brca_paclitaxel', 'GEO_brca_docetaxel']  :
        active_modules_map = load_modules('.', cp_active_modules_file_brca)
    else :       
        active_modules_map = load_modules('.', cp_active_modules_file)
        
    patient_up_genes = load_gene_list('.', patient_up_genes_file)
    patient_down_genes = load_gene_list('.', patient_down_genes_file)


    active_modules_map_reduced_with_patient_up_genes = reduce_map_lengths(active_modules_map, patient_up_genes)
    active_modules_map_reduced_with_patient_down_genes = reduce_map_lengths(active_modules_map, patient_down_genes)


    printmap(active_modules_map_reduced_with_patient_up_genes, 'up', 10)
    printmap(active_modules_map_reduced_with_patient_down_genes, 'down', 10)


    #print(os.listdir(path_up))
    files_up = sorted([f for f in os.listdir(path_up) if os.path.isfile(os.path.join(path_up, f)) and f.startswith('enrichR_module')])
    printlist(files_up, 'lista files up', 10)

    pairs_up = [(path_up, x) for x in files_up ]


    files_down = sorted([f for f in os.listdir(path_down) if os.path.isfile(os.path.join(path_down, f))  and f.startswith('enrichR_module')])
    printlist(files_down, 'lista files down', 10)

    pairs_down = [(path_down, x) for x in files_down ]

    verbose = 'verbose'
    #weights = 'log10'
    #weights = 'c_score'
    #weights = 'deltadeg'
    #weights = 'pval_x_deltadeg'
    #weights = 'pval_x_drug_hits'
    #weights = 'pval_x_deltadeg_x_drug_hits'
    #weights = 'odds_ratio_log'
    #weights = 'c_score_x_deltadeg'
    weights = ['log10', 'c_score', 'deltadeg', 'pval_x_deltadeg', 'pval_x_drug_hits',  'pval_x_deltadeg_x_drug_hits',  'odds_ratio_log', 'c_score_x_deltadeg' ]

    drugnames = set()
    drug_module_map = {}
    counter = 0
    fdaonly = True
#    drugdatatype = 'GEO'
#    drugdatatype = 'CMAP'
#    drugdatatype = 'L1000'
    pvalue_filter = True
#    p_value_limit = 0.05
#    p_value_limit = 0.01
#    p_value_limit = 0.001
#    p_value_limit = 0.0001
#    p_value_limits = [0.0001, 0.001, 0.01, 0.05 ]
    
    #counter_limit = 10
    counter_limit = 1000000

    pairs =  pairs_up + pairs_down

#    for p_value_limit in p_value_limits :
    print("@@ Starting p_value : {}".format(p_value_limit))
    for x in pairs:
        module_id = extract_module_number(x[1])
        counter +=1
        print('Counter {} , module n: {} '.format(counter, module_id ))

        if x[0].endswith('up') :
            delta_deg = active_modules_map_reduced_with_patient_up_genes[module_id] - active_modules_map_reduced_with_patient_down_genes[module_id]
        elif x[0].endswith('down') :
            delta_deg = active_modules_map_reduced_with_patient_down_genes[module_id] - active_modules_map_reduced_with_patient_up_genes[module_id]
        else :
            print ('Wrong suffix in ', x[0])
            assert False
      

        acc ,  header = load_data(x[0], x[1], ',', 'noverbose')
        if pvalue_filter :
            pval_filtered_dataarray = [z for z in acc  if test_pval(z[record_map['P.value']], p_value_limit)]
            acc = pval_filtered_dataarray
        
        for y in acc :
            term = y[record_map['Term']]
            if term == '' :
                continue
            p_val = float(y[record_map['P.value']])
            c_score = float(y[record_map['Combined.Score']])
            DE_drug_num = len(y[record_map['Genes']].split(';'))
            minuslog = minuslog10(p_val)
            odds_ratio_log = math.log10(float(y[record_map['Odds.Ratio']]))
            if drugdatatype1 != 'GEO' :
                drug_name = extract_drug_name_only(term, drugdatatype1 + 'up')
            elif drugdatatype1 == 'GEO' :
                (drug_name, direction) = extract_drug_name_direction(term , drugdatatype1)
                if x[0].endswith(direction) :  # check if they have he same direction skip this record
#                    print('Skipping {}, {}'.format(x[0], term))
                    continue
                if drug_name == None :
                    continue
                
            if verbose == 'verbose' :
                new_record =  (p_val, minuslog, c_score, delta_deg, DE_drug_num, odds_ratio_log)
    #            print(drug_name, p_val, minuslog)
        
            drugnames.add(drug_name)
            if (module_id, drug_name) not in   drug_module_map.keys() :        
                drug_module_map[(module_id, drug_name)] = [new_record]
            else :
                drug_module_map[(module_id, drug_name)].append(new_record)

        # early escape
        if counter >= counter_limit :
            break

    print('Counter of file opened: ', counter, ' over ', len(pairs))
    print ("Number of drugs " , len(drugnames))

    # initialize mappings to hold the reorganized results by drugname

    for weight in weights :
        drugscore = {}
        drug_modules = {}
        for x in drugnames :
            drugscore[x] = 0
            drug_modules[x] = []


        for x in drug_module_map.keys() :
        #        drugdirection = drug_module_map[x][0][2] #look at the first record only
            drugname = x[1]

            drugscore[drugname] +=  score(drug_module_map[x], 'deg_in_module', weight)
            drug_modules[drugname].append(x[0]) # record the module id


        drugscore_items = drugscore.items()
        drugscore_items_sorted = sorted(drugscore_items, key=lambda tup: tup[1], reverse=True)
        printlist(drugscore_items_sorted, 'drugscore_items_sorted', 20)


        # fda filter
        if fdaonly == True :
            list_of_drug_names = [x for x in drugscore.keys()]
            fda_filtered_list_of_drug_names = filterbyfdaapproveddrugs(list_of_drug_names, drugdatatype1)
            drugscore_fda_filtered = domain_restriction(fda_filtered_list_of_drug_names, drugscore, 'nosubset')
            drugscore_items = drugscore_fda_filtered.items()
            drugscore_items_sorted = sorted(drugscore_items, key=lambda tup: tup[1], reverse=True)
            printlist(drugscore_items_sorted, 'drugscore_items_sorted after fda filtering', 20)

            if drugdatatype in ['GEO_brca_fluorouracil', 'GEO_brca_doxorubicin', 'GEO_brca_paclitaxel', 'GEO_brca_docetaxel'] :
                conf = (patient_id ,  'cpd07j08', drugdatatype1, 'nofilter', 'brca', )
            else : 
                conf = (patient_id ,  'cpd07j08', drugdatatype1, 'nofilter', 'crc', )
                
            new_ranking_evaluation_via_TTD_data(conf, drugscore_items_sorted, fdaonly, weight)

        # check if therapy drigs dta is available form apient and in teh matrix
            if conf[0] not in drugs_to_check.keys() :
                print('No drugs_to_check for {}'.format(conf[0]))
            elif drugs_to_check[conf[0]] == [] :
                 print('No drugs_to_check for {}'.format(conf[0]))
            else :
#                drug_names_in_matrix = extract_drug_names_from_file(conf, verbose = 'noverbose')
                drug_names_in_matrix = list_of_drug_names
#                filename = modulation_data[conf][1]
                drug_format = conf[2]
#                target_drugs = filter_two_lists(drugs_to_check_lable[conf[0]], drugs_to_check[conf[0]], filename, drug_names_in_matrix )
                target_drugs = filter_two_lists(drugs_to_check_lable[conf[0]], drugs_to_check[conf[0]], conf[0], drug_names_in_matrix )
                trial_drugs = trial_drugs_TTD_data(conf, drugscore_items_sorted, fdaonly)
                trial_drugs_2 = list(set(trial_drugs) | set(target_drugs))
                check_specific_drugs(drugscore_items_sorted, target_drugs, drug_format, trial_drugs_2, weight)


if __name__ == "__main__":
    #    pass
    main_function(patient_id, drugdatatype , p_value_limit)
