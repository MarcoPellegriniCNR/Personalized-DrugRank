# Code to collect the support clincial data for the CRC subpopulations

from generic_util import venndiagramsets , printmap, printlist, increment_lists
from new_loader_routines import load_data
import statistics
import scipy.stats



########
# Flags

run_gen_CRC_maps = True

# subpopulation control flags

#sub_population_flag = 'all'
#sub_population_flag = 'fl_only'
#sub_population_flag = 'fl_ox'
#sub_population_flag = 'fl_leu'

try:
    sub_population_flag
except NameError:
    print("sub_population_flage is not defined: please edit the script uncommenting one of the possible values!")
else:
    print("sub_population_flag is  defined.")

print('sub_population_flag : ', sub_population_flag)

######


################
# compute duration of the therapy 
def compute_duration(start_time, end_time ):

#    print ('start, end', start_time, end_time)

    if start_time in ['' , "'--"] or end_time in ['' , "'--"] :
        return None
    start = int(start_time)
    end = int(end_time)
    assert end >= start
    if end - start  <= 0 :
        return 1  #a conventional small postive duration 
    else:
        return end - start



#####################
# analysis of the fiamma 43 florouracil  CRC patients clincial data:

def gen_drugs_to_check_maps() :
    
    verbose = False

    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'

    CLINICAL_DATA_FILE_NAME = 'dati_clinici_pz_fluorouracil_CRC.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, ';', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index('\ufeffcase_id')
    treatment_type_idx =   header.index( 'clinical (2).treatments.treatment_type')
    treatment_idx  =   header.index('clinical (2).treatments.therapeutic_agents')
    response_idx = header.index('clinical (2).treatments.treatment_outcome')

    start_time_idx = header.index( 'clinical (2).treatments.days_to_treatment_start')
    end_time_idex = header.index( 'clinical (2).treatments.days_to_treatment_end')

    treatments_map = {}
    resp_map = {}
    reduced_response_map ={}
    pat_set = set()
    drug_set = set()
    resp_set = set()
    for x in acc :
        pat_id = x[patient_idx]
        if pat_id == '' :
            continue
        treatment = x[treatment_idx]
        treatment_type = x[treatment_type_idx]
        response = x[response_idx]
        start_time = x[start_time_idx]
        end_time = x[end_time_idex]
        duration = compute_duration(start_time, end_time )
        if treatment_type ==  'Chemotherapy'  and treatment != 'Not Reported' :
            if verbose :
                print (pat_id, treatment, treatment_type, response)
            pat_set.add(pat_id)
            drug_set.add(treatment)
            resp_set.add(response)
            increment_lists(treatments_map, pat_id, treatment)
            if response not in ['"-"', '"--"' , "'--",  'Unknown', 'Treatment Ongoing'] :
                increment_lists(resp_map, pat_id, (treatment,response))
                if treatment == 'Fluorouracil' :
                    increment_lists(reduced_response_map, pat_id, (response, duration))

    print('pat set', len(pat_set), pat_set)
    print('drug set' , len(drug_set) , drug_set)
    print('responses', len(resp_set), resp_set)


    f= open('fiamma_list.tex')
    fiamma_list = f.read().split()

#    printlist(fiamma_list, 'fiamma_list')

    venndiagramsets('fiamma_list', set(fiamma_list), 'clinical_fluorouracil_data', pat_set)
      
    printmap(treatments_map, 'treatments_map')
    printmap(resp_map, 'responses_map')
    printmap(reduced_response_map, 'reduced responses_map')

    drugs_to_check_lable = {}
    for key in treatments_map.keys():
        drugs_to_check_lable[key] = 'crc_drug_used_for_'+ key

    return(treatments_map, drugs_to_check_lable, reduced_response_map)

####################
# run

if run_gen_CRC_maps :
    (CRC_treatments_map, CRC_drugs_to_check_lable, CRC_reduced_response_map) = gen_drugs_to_check_maps()


def filter_with_without(inlist, with_drug_names, without_drug_names) :
    match = False
    for x in inlist :
        if x in without_drug_names :
            return False
    for y in with_drug_names :
        absent = True
        for x in inlist :
            if x == y :
                absent = False
        if absent :
            return False
    return True


# returns the list of patients with specific drug patterns.


def patients_by_groups (in_treatment_map ,  with_drug_names, without_drug_names) :
    print('*************************************')
    print('Start lsiting patietns by groups')
    out_pat_list = []
    # finding rousp of patients by therapy
    for pat in in_treatment_map.keys() :
        if filter_with_without( in_treatment_map[pat]  , with_drug_names, without_drug_names) :
            out_pat_list.append((pat, in_treatment_map[pat]))
    printlist( out_pat_list  , 'out_pat_list  with {} , without{}'.format(with_drug_names, without_drug_names))
    print('*************************************')
    return [x[0]  for x in out_pat_list]


#######################
# this function apply a restcion to the drugmerge mapping and produces three objects need i the cox analysis

def new_subpopulation_data (in_main_map, sub_class_list, subclass_string):

    print('Processing', subclass_string)
    new_map = {}
    
    for x in sub_class_list:
        if x in in_main_map.keys():
            print('Key {} present'.format(x))
            new_map[x] = in_main_map[x]
        else :
            print ('Key {} NOT present'.format(x))

    printmap(new_map, 'P_Drug_rank data for: ' + subclass_string)

    good_prediction= [key for key  in new_map.keys()  if float(new_map[key][1])/  float(new_map[key][2]) < 0.10   ]
    printlist(good_prediction, 'good_prediction :' + subclass_string)

    bad_prediction=  [key for key  in new_map.keys() if float(new_map[key][1])/  float(new_map[key][2]) >=  0.10   ]
    printlist(bad_prediction , 'bad_prediction' +  subclass_string)

    #export as  new global variable to be set to use this data.
    return(new_map, good_prediction, bad_prediction)



######################
# main pdrugmerge indeices for CRC fluorouracil patients


map_fluorouracil_geo_results = {
'TCGA-A6-A566' :  (0.175, 39, 92, 4, 6, 'c_score_x_deltadeg',  0.05),
'TCGA-A6-A56B' :  (0.15, 12, 93, 1, 6, 'odds_ratio_log', 0.05),
'TCGA-A6-A5ZU' : (0.12, 36, 95,2, 6, 'c_score_x_deltadeg',  0.05),
'TCGA-AA-3971' :  (0.18,12,96,1,7,'odds_ratio_log', 0.05),
'TCGA-AA-A00J' :  (0.337, 3, 64, 1, 5, 'odds_ratio_log', 0.001),
'TCGA-AA-A00Q' : (0.178,18, 97, 1, 7,  'c_score_x_deltadeg',  0.05),
'TCGA-AA-A010' :  (0,20, 12, 90, 1, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-AA-A01F' : (0.256, 22, 62, 2, 5, 'c_score_x_deltadeg',  0.001),
'TCGA-AA-A01K' : (0.19, 18,93, 1, 6, 'c_score_x_deltadeg',  0.05),
'TCGA-AA-A01T' : (0.227, 16, 96, 2, 7,  'c_score_x_deltadeg',  0.05),
'TCGA-AD-6965' : (0.401, 5, 96, 1, 7, 'deltadeg', 0.05), # irinotecan
'TCGA-AF-6136' : (0.38, 6, 96, 2, 6,'deltadeg', 0.05),
'TCGA-AF-A56L' : (0.40,6, 95,2, 7, 'deltadeg', 0.05),
'TCGA-AG-3593' : (0.364, 25, 95, 3, 7, 'c_score_x_deltadeg',  0.05),
'TCGA-AG-3600' : (0.457, 5, 81, 2,6,  'deltadeg', 0.01),
'TCGA-AG-3609' : (0.236, 29, 93, 3, 6,'c_score_x_deltadeg',  0.05),
'TCGA-AG-4022' : (0.468, 22, 95, 4, 5,  'c_score_x_deltadeg',  0.05),
'TCGA-AG-A01W' : (0.338, 5, 78, 1, 5,  'deltadeg', 0.01),
'TCGA-AG-A01Y' : (0.268, 11, 64, 1, 5,  'c_score_x_deltadeg',  0.001),
'TCGA-AG-A020' : (0.370, 7, 95, 2, 7, 'pval_x_deltadeg ', 0.05),
'TCGA-AH-6544' : (0.417, 7, 95, 2, 7, 'pval_x_deltadeg ', 0.05),
'TCGA-AH-6643' : (0.319, 6, 69, 1, 4,   'deltadeg', 0.001),
'TCGA-AH-6903' : (0.337, 9, 96, 2, 6, 'deltadeg', 0.05),
'TCGA-AY-A8YK' : (0.287, 17, 94, 2, 6, 'c_score_x_deltadeg',  0.05),
'TCGA-AZ-6606' : (0.393, 5, 98, 1, 6, 'pval_x_deltadeg ', 0.05), # irinotecan
'TCGA-CA-5254' : (0.534, 4, 94, 1, 5, 'c_score_x_deltadeg',  0.05),
'TCGA-CA-5255' : (0.288, 9, 94, 1, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-CA-5256' : (0.240, 18, 95, 1, 7,   'c_score_x_deltadeg',  0.05),
'TCGA-CA-5796' : (0.172, 27, 97, 2, 7,  'c_score_x_deltadeg',  0.05),
'TCGA-CA-6715' : (0.52, 7, 97, 3, 7,  'c_score_x_deltadeg',  0.05),
'TCGA-CA-6716' : (1.13, 65, 94, 6, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-CA-6717' : (0.369, 11, 92, 3, 5, 'c_score_x_deltadeg',  0.05),
'TCGA-CK-4950' : (0.400, 7,96, 2, 6, 'deltadeg', 0.05),
'TCGA-F4-6856' : (1.13, 19, 70, 2, 5, 'c_score', 0.001),
'TCGA-F4-6569' : (0.493, 7, 94, 2, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-F5-6702' : (0.64, 53, 94, 6, 6,   'c_score_x_deltadeg',  0.05),
'TCGA-F5-6863' : (0.378, 8,93, 2, 6,   'deltadeg', 0.05),
'TCGA-F5-6864' : (0.401, 34, 94, 4, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-NH-A50U' : (0.370, 23, 94, 3, 5, 'pval_x_deltadeg_x_drug_hits ', 0.05),
'TCGA-NH-A50V' : (0.228, 60, 94, 4, 6, 'c_score', 0.05),
'TCGA-NH-A6GB' : (0.346, 5,66, 1, 5,  'deltadeg', 0.001),
'TCGA-NH-A8F8' : (0.641, 26, 95,3, 5, 'c_score_x_deltadeg',  0.05),
'TCGA-SS-A7HO' : (0.337, 6, 95, 1, 6, 'pval_x_deltadeg', 0.05)
}


####################
# subpopulations by treatment

if  sub_population_flag == 'all' :
# use all the patients
    subclass_string = 'All'
    reduced_list_all = CRC_treatments_map.keys()

    ( map_CRC_fluorouracil_geo_results, list_CRC_fluorouracil_good_geo_prediction_ids, list_CRC_fluorouracil_bad_geo_prediction_ids) = new_subpopulation_data (map_fluorouracil_geo_results, reduced_list_all, subclass_string)

elif sub_population_flag == 'fl_only' :
    # fluorouracil only
    subclass_string = 'Subpopulation fluorouracil only'
    with_list = [ 'Fluorouracil']
    without_list = ['Irinotecan', 'Etoposide', 'Leucovorin Calcium', 'Regorafenib', 'Leucovorin', 'Bevacizumab', 'Oxaliplatin', 'Capecitabine']
    reduced_lsit_op_patientes_fl_only = patients_by_groups (CRC_treatments_map ,  with_list, without_list) 

    ( map_CRC_fluorouracil_geo_results, list_CRC_fluorouracil_good_geo_prediction_ids, list_CRC_fluorouracil_bad_geo_prediction_ids) =  new_subpopulation_data(map_fluorouracil_geo_results,reduced_lsit_op_patientes_fl_only, subclass_string)


elif sub_population_flag == 'fl_ox' :
    # fluorouracil and oxaliplatin plus other..but not irinotecan
    subclass_string = 'Subpopulation fluorouracil and oxaliplatin  without Leucovorin'
    with_list = [ 'Fluorouracil',  'Oxaliplatin']
    without_list = ['Irinotecan', 'Leucovorin Calcium',  'Leucovorin']
    reduced_lsit_op_patientes_fl_ox = patients_by_groups (CRC_treatments_map ,  with_list, without_list)

    ( map_CRC_fluorouracil_geo_results, list_CRC_fluorouracil_good_geo_prediction_ids, list_CRC_fluorouracil_bad_geo_prediction_ids) =  new_subpopulation_data(map_fluorouracil_geo_results, reduced_lsit_op_patientes_fl_ox, subclass_string)


elif sub_population_flag == 'fl_leu' :
    # fluorouracil and Leucovorin Calcium plus other..but not irinotecan
    subclass_string = 'Subpopulation fluorouracil and Leucovorin'
    with_list = [ 'Fluorouracil',  'Leucovorin Calcium' ]
    without_list = ['Irinotecan']
    reduced_lsit_op_patientes_fl_leu1 = patients_by_groups (CRC_treatments_map ,  with_list, without_list)

    # fluorouracil and Leucovorin  plus other..but not irinotecan
    with_list = [ 'Fluorouracil' , 'Leucovorin']
    without_list = ['Irinotecan']
    reduced_lsit_op_patientes_fl_leu2 = patients_by_groups (CRC_treatments_map ,  with_list, without_list)

    reduced_lsit_op_patientes_fl_leu = reduced_lsit_op_patientes_fl_leu1 + reduced_lsit_op_patientes_fl_leu2

    ( map_CRC_fluorouracil_geo_results, list_CRC_fluorouracil_good_geo_prediction_ids, list_CRC_fluorouracil_bad_geo_prediction_ids) =  new_subpopulation_data(map_fluorouracil_geo_results, reduced_lsit_op_patientes_fl_leu, subclass_string)


else :
    print('Wrong subpopulation lable {}'.format(sub_population_flag ))
    assert False
    

##########################





