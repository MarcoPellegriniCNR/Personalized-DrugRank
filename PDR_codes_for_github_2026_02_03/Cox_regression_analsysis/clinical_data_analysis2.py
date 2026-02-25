### code to analyze the clinical data

from generic_util import venndiagramsets , printmap, printlist, increment_lists
from new_loader_routines import load_data


########

    
def get_patient_map(lista, p_index, target_index) :

    outmap ={}
    for x in lista:
        outmap[x[p_index]]=x[target_index]
    return outmap


def generate_maps() :

    CLINICAL_DATA_PATH_NAME = '..'
    #CLINICAL_DATA_FILE_NAME = 'Irinitecan_coadread_tcga_pan_can_atlas_2018_clinical_data.csv'
    CLINICAL_DATA_FILE_NAME = 'Fluorouracil_coadread_tcga_pan_can_atlas_2018_clinical_data.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, '\t', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index( 'Patient ID')

    #pfs analysis
    pfs_months_idx = header.index('Progress Free Survival (Months)')
    pfs_status_idx = header.index('Progression Free Status')
    event_string_pfs =  '1:PROGRESSION'

    #os analysis
    os_months_idx = header.index('Overall Survival (Months)')
    os_status_idx = header.index('Overall Survival Status')
    event_string_os = '1:DECEASED'

    #dss
    dss_months_idx = header.index('Months of disease-specific survival')
    dss_status_idx = header.index('Disease-specific Survival status')
    event_string_dss = '1:DEAD WITH TUMOR'


    print('patient_idx', patient_idx, 'pfs_months_idx', pfs_months_idx, 'pfs_status_idx', pfs_status_idx)



    verbose = False

    map_pfs_months = get_patient_map(acc, patient_idx, pfs_months_idx)
    map_pfs_status = get_patient_map(acc, patient_idx, pfs_status_idx)
    if verbose :
        print('map_pfs_months', map_pfs_months)
        print('map_pfs_status', map_pfs_status)

    map_os_months = get_patient_map(acc, patient_idx, os_months_idx)
    map_os_status = get_patient_map(acc, patient_idx, os_status_idx)
    if verbose :
        print('map_os_months', map_os_months)
        print('map_os_status', map_os_status)

    map_dss_months = get_patient_map(acc, patient_idx, dss_months_idx)
    map_dss_status = get_patient_map(acc, patient_idx, dss_status_idx)
    if verbose :
        print('map_dss_months', map_dss_months)
        print('map_dss_status', map_dss_status)

    return (map_pfs_months, map_pfs_status, map_os_months , map_os_status,  map_dss_months,  map_dss_status)

############ CMAP #####
list_irinotecan_good_cmap_prediction_ids = [
'TCGA-AZ-6600',
'TCGA-AA-3679',
'TCGA-A6-2671',
'TCGA-AA-3869',
'TCGA-CM-5868',
'TCGA-NH-A6GA',
'TCGA-AA-A02F',
'TCGA-AA-3548',
'TCGA-A6-5664',
'TCGA-AA-A01D',
'TCGA-CM-6675',
'TCGA-CM-6678',
'TCGA-G4-6297',
'TCGA-G4-6303',
'TCGA-F5-6571'
, # new data
'TCGA-DC-6156',
'TCGA-AG-3887',
'TCGA-AG-3731',
'TCGA-AG-3584',
]

list_irinotecan_bad_cmap_prediction_ids = [
    'TCGA-AA-3517',
#    'TCGA-AA-3869', # duplicate
    'TCGA-AY-4070',
    'TCGA-AA-3494',
    'TCGA-CM-4748',
    'TCGA-A6-2683',
    'TCGA-F4-6806',
    'TCGA-AA-3973',
    'TCGA-AA-A01R',
    'TCGA-AG-A02G', # new data
    ]

############# GEO #############
list_irinotecan_good_geo_prediction_ids = [
'TCGA-AZ-6600',
'TCGA-AA-3679',

'TCGA-AA-3869',
'TCGA-CM-5868',
'TCGA-NH-A6GA',
'TCGA-AA-A02F',
'TCGA-AA-3548',
'TCGA-A6-5664',
'TCGA-AA-A01D',
'TCGA-CM-6675',
'TCGA-CM-6678',
'TCGA-G4-6297',
'TCGA-G4-6303',
 # new data
'TCGA-AG-3887',
'TCGA-AG-3731',
'TCGA-AG-3584',
# from bad cmap
 'TCGA-AA-3517',
'TCGA-AA-A01R',
]


# mediaum

     
# bad


list_irinotecan_bad_geo_prediction_ids = [

#    'TCGA-AA-3869', # duplicate
    'TCGA-AY-4070',
    'TCGA-AA-3494',
    'TCGA-CM-4748',
    'TCGA-A6-2683',
    'TCGA-F4-6806',
    'TCGA-AA-3973',
    'TCGA-AG-A02G', # new data
'TCGA-A6-2671', # from good cmap 
'TCGA-F5-6571',
'TCGA-DC-6156',
    ]

import statistics
import scipy.stats

# uses only  data for occurring evnet (not the censored data)

def ttest_ks_analysis(good_predictions_id_list, bad_predictions_id_list, time_to_event_map, status_event_map, event_string, logging_string) :

    l_good =[float(time_to_event_map[x]) for x in good_predictions_id_list  if status_event_map[x] ==event_string]
    l_bad= [float(time_to_event_map[x]) for x in bad_predictions_id_list if status_event_map[x] == event_string]

    print( "For {}, good prediction {}, with event {}, bd predictions {}, with event {} ".format( logging_string, len(good_predictions_id_list),   len(l_good), len(bad_predictions_id_list),  len(l_bad)))

    print( logging_string,  'mean event time for good predictions' , statistics.mean(l_good))
    print( logging_string,  'mean event time for bad predictions' , statistics.mean(l_bad))

    res_tt = scipy.stats.ttest_ind(l_good, l_bad)
    res_ks = scipy.stats.ks_2samp(l_good, l_bad)

    print(res_tt)
    print(res_ks)


from lifelines.statistics import logrank_test
import numpy as np

# Example data:

def log_rank_analysis(good_predictions_id_list, bad_predictions_id_list, time_to_event_map, status_event_map, event_string, logging_string):

    et1=[]
    eo1 = []
    for x in good_predictions_id_list :
        et1.append(float(time_to_event_map[x]))
        if status_event_map[x] == event_string :
            eo1.append(True)
        else :
            eo1.append(False)
            
    et2=[]
    eo2 = []
    for x in bad_predictions_id_list :
        et2.append(float(time_to_event_map[x]))
        if status_event_map[x] == event_string :
            eo2.append(True)
        else :
            eo2.append(False)
      


    event_times1 = np.array(et1)
    event_observed1 = np.array(eo1)
    event_times2 = np.array(et2)
    event_observed2 = np.array(eo2)
    print(logging_string, 'good predictions', et1, eo1)
    print(logging_string, 'bad predictions', et2, eo2)

    results = logrank_test(event_times1, event_times2, event_observed1, event_observed2)
#    print(results)
    print("For {} Log-rank p-value: {}".format(logging_string, results.p_value))

## 
# run
##logging_string_1 = 'CMAP irinotecan pfs analysis'
##log_rank_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_pfs_months, map_pfs_status, event_string_pfs, logging_string_1)
##ttest_ks_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_pfs_months, map_pfs_status, event_string_pfs, logging_string_1)
##    
##    
##logging_string_2 = 'CMAP irinotecan os analysis'
##log_rank_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_os_months, map_os_status, event_string_os, logging_string_2)
##ttest_ks_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_os_months, map_os_status, event_string_os, logging_string_2)
##
##logging_string_3 = 'CMAP irinotecan dss analysis'
##log_rank_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_dss_months, map_dss_status, event_string_dss, logging_string_3)
##ttest_ks_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_dss_months, map_dss_status, event_string_dss, logging_string_3)

##
##logging_string_4 = 'GEO irinotecan pfs analysis'
##log_rank_analysis(list_irinotecan_good_geo_prediction_ids, list_irinotecan_bad_geo_prediction_ids, map_pfs_months, map_pfs_status, event_string_pfs, logging_string_4)
##ttest_ks_analysis(list_irinotecan_good_geo_prediction_ids, list_irinotecan_bad_geo_prediction_ids, map_pfs_months, map_pfs_status, event_string_pfs, logging_string_4)
##
##logging_string_5 = 'GEO irinotecan os analysis'
##log_rank_analysis(list_irinotecan_good_geo_prediction_ids, list_irinotecan_bad_geo_prediction_ids, map_os_months, map_os_status, event_string_os, logging_string_5)
##ttest_ks_analysis(list_irinotecan_good_cmap_prediction_ids, list_irinotecan_bad_cmap_prediction_ids, map_os_months, map_os_status, event_string_os, logging_string_5)
##
##logging_string_6 = 'GEO irinotecan dss analysis'
##log_rank_analysis(list_irinotecan_good_geo_prediction_ids, list_irinotecan_bad_geo_prediction_ids, map_dss_months, map_dss_status, event_string_dss, logging_string_6)
##ttest_ks_analysis(list_irinotecan_good_geo_prediction_ids, list_irinotecan_bad_geo_prediction_ids, map_dss_months, map_dss_status, event_string_dss, logging_string_6)


##########
#

def comparing_lists () :
    fluorouracil_pats = list(map_pfs_months.keys())

    irinotecan_pats = list_irinotecan_good_cmap_prediction_ids + list_irinotecan_bad_cmap_prediction_ids

    fluorouracil_not_irinotecan_pats = list(set(fluorouracil_pats) - set(irinotecan_pats))

    #printlist(fluorouracil_pats, 'fluorouracil_pats')
    printlist(fluorouracil_not_irinotecan_pats, 'fluorouracil_not_irinotecan_pats')

    f= open('fiamma_list.tex')
    fiamma_list = f.read().split()

    printlist(fiamma_list, 'fiamma_list')

    venndiagramsets('fiamma_lsit', set(fiamma_list), 'fluorouracil_not_irinotecan_pats', set(fluorouracil_not_irinotecan_pats))


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

#generate_maps()
#(treatments_map, drugs_to_check_lable, reduced_response_map) = gen_drugs_to_check_maps()


######################################################
# uploading data for the BRCA patients.

def BRCA_gen_drugs_to_check_maps(target_drugs) :

# To be changed for each drug target on BRCA data 
#    target_drugs = [ 'Fluorouracil']
#    target_drugs =['Doxorubicin Hydrochloride', 'Doxorubicin']
 #   target_drugs =[ 'Paclitaxel', ] 
    verbose = False
#    verbose = True

    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'

    CLINICAL_DATA_FILE_NAME = 'BRCA_all_treatments_2.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, '\t', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index('cases.submitter_id')
    treatment_type_idx =   header.index( 'treatments.treatment_type')
    treatment_idx  =   header.index('treatments.therapeutic_agents')
    response_idx = header.index('treatments.treatment_outcome')

    start_time_idx = header.index( 'treatments.days_to_treatment_start')
    end_time_idex = header.index( 'treatments.days_to_treatment_end')

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
        if treatment_type in ['Hormone Therapy',   'Chemotherapy' ] and treatment != 'Not Reported' :
            if verbose :
                print (pat_id, treatment, treatment_type, response)
            pat_set.add(pat_id)
            drug_set.add(treatment)
            resp_set.add(response)
            increment_lists(treatments_map, pat_id, treatment)
            if response not in ['"-"', '"--"' , "'--",  'Unknown', 'Treatment Ongoing'] :
                increment_lists(resp_map, pat_id, (treatment,response))
                if treatment in target_drugs :
                    increment_lists(reduced_response_map, pat_id, (response, duration))

    print('pat set', len(pat_set), pat_set)
    print('drug set' , len(drug_set) , drug_set)
    print('responses', len(resp_set), resp_set)


##    f= open('fiamma_list.tex')
##    fiamma_list = f.read().split()
##
###    printlist(fiamma_list, 'fiamma_list')
##
##    venndiagramsets('fiamma_list', set(fiamma_list), 'clinical_fluorouracil_data', pat_set)
      
    printmap(treatments_map, 'treatments_map')
    printmap(resp_map, 'responses_map')
    printmap(reduced_response_map, 'reduced responses_map')

    drugs_to_check_lable = {}
    for key in treatments_map.keys():
        drugs_to_check_lable[key] = 'BRCA_drug_used_for_'+ key

    return(treatments_map, drugs_to_check_lable, reduced_response_map)




####################
# run                           

#generate_maps()
target_drugs = ['Fluorouracil']
#target_drugs = ['Docetaxel']
#target_drugs =  ['Paclitaxel']

try:
    target_drugs
except NameError:
    print("target_drugs is not defined: please edit the script uncommenting one of the possible values!")
else:
    print("target_drugs is  defined.")

print('target_drugs : ', target_drugs)
 
(treatments_map, drugs_to_check_lable, reduced_response_map) = BRCA_gen_drugs_to_check_maps(target_drugs)

# getting the repsonse data:
# vero se tutti i drugs in with_drugnames sono contenuti
# falso se almeno una delle witout_drugnames sono contenuti

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

def patients_by_groups ( with_drug_names, without_drug_names) :
    print('*************************************')
    print('Start lsiting patietns by groups')
    out_pat_list = []
    # finding rousp of patients by therapy
    for pat in treatments_map.keys() :
        if filter_with_without( treatments_map[pat]  , with_drug_names, without_drug_names) :
            out_pat_list.append((pat, treatments_map[pat]))
    printlist( out_pat_list  , 'out_pat_list  with {} , without{}'.format(with_drug_names, without_drug_names))
    print('*************************************')
    return [x[0]  for x in out_pat_list]


               
#



# mappa dei risutlati seprimentali
# max score, rank in all drugs, all drugs, rank in crc durgs , all crc drugs, function name, p-val-filter0.173

map_flourouracil_geo_results = {
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

list_CRC_fluorouracil_good_geo_prediction_ids = [key for key  in map_flourouracil_geo_results.keys() if float(map_flourouracil_geo_results[key][1])/  float(map_flourouracil_geo_results[key][2]) < 0.10   ]
printlist(list_CRC_fluorouracil_good_geo_prediction_ids, 'list_CRC_fluorouracil_good_geo_prediction_ids')

list_CRC_fluorouracil_bad_geo_prediction_ids =  [key for key  in map_flourouracil_geo_results.keys() if float(map_flourouracil_geo_results[key][1])/  float(map_flourouracil_geo_results[key][2]) >=  0.10   ]
printlist(list_CRC_fluorouracil_bad_geo_prediction_ids, 'list_CRC_fluorouracil_bad_geo_prediction_ids')


#############
# similar data for irinotecan patients.

map_irinotecan_geo_results = {
'TCGA-AZ-6600' : ( 0.359,  6, 95,  1,  5, 'pval_x_deltadeg' , 0.05),
'TCGA-AA-3679' : ( 0.381, 3, 96, 1, 7,  'pval_x_deltadeg' , 0.05),																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																				
'TCGA-A6-2671' : (0.176, 29, 93, 2, 6, 'c_score_x_deltadeg',  0.05),
'TCGA-AA-3869' : (0.297, 6, 96, 1, 7,   'pval_x_deltadeg' , 0.05),
'TCGA-CM-5868' : (0.350, 5, 92,1,6,  'pval_x_deltadeg' , 0.05),
'TCGA-NH-A6GA' : (0.368, 5, 81, 1, 5, 'pval_x_deltadeg' , 0.01),
'TCGA-AA-A02F' : (0.323, 5, 95, 1, 6 , 'deltadeg' , 0.05),
'TCGA-AA-3548' :( 0.249, 10, 95, 1, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-A6-5664' : (0.401, 5, 95, 1, 6, 'pval_x_deltadeg' , 0.05),
'TCGA-AA-A01D' : (0.254, 17, 97, 2, 7, 'c_score_x_deltadeg',  0.05),
'TCGA-CM-6675' : (0.420, 5, 96, 1, 7, 'pval_x_deltadeg' , 0.05),
'TCGA-CM-6678' : (0.538, 2, 96, 1, 7, 'c_score_x_deltadeg',  0.05),
'TCGA-G4-6297' : (0.382, 4, 96, 1, 8,  'deltadeg' , 0.05),
'TCGA-G4-6303' : (0.328, 7, 85, 1,6, 'pval_x_deltadeg' , 0.01),
'TCGA-F5-6571' : (1.184, 30, 95, 3, 7,  'c_score_x_deltadeg',  0.05),
'TCGA-DC-6156' : (0.860, 21, 94, 4, 5, 'c_score_x_deltadeg',  0.05),
'TCGA-AG-3887' : (0.345, 13, 97, 2, 7, 'pval_x_deltadeg' , 0.05),
'TCGA-AG-3731' : (0.411, 5, 94, 1, 6,  'pval_x_deltadeg' , 0.05),
'TCGA-AG-3584' :  (0.282, 5, 80, 1, 5, 'pval_x_deltadeg' , 0.01),
#
'TCGA-AA-3517' : (0.157, 14, 94, 1, 6, 'odds_ratio_log',  0.05),
'TCGA-AY-4070' : (0.184, 20, 96, 1, 7,  'c_score_x_deltadeg',  0.05),
#'TCGA-AA-3494' :  # no data
'TCGA-CM-4748' : (0.163, 18, 96, 1,8 , 'odds_ratio_log',  0.05),
'TCGA-A6-2683' : (0.156, 25, 95, 2, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-F4-6806' : (0.264, 24, 95, 2, 6,  'c_score_x_deltadeg',  0.05),
'TCGA-AA-3973' : (0.133, 35, 96, 2, 6, 'c_score_x_deltadeg',  0.05),
'TCGA-AA-A01R': (0.246, 8, 94, 1, 7, 'c_score_x_deltadeg',  0.05),
'TCGA-AG-A02G' : (0.493, 17, 95, 2, 6, 'c_score_x_deltadeg',  0.05),
}


list_CRC_irinotecan_good_geo_prediction_ids_2 = [key for key  in map_irinotecan_geo_results.keys() if float(map_irinotecan_geo_results[key][1])/  float(map_irinotecan_geo_results[key][2]) < 0.10   ]
printlist(list_CRC_irinotecan_good_geo_prediction_ids_2, 'list_CRC_irinotecan_good_geo_prediction_ids_2')
##
list_CRC_irinotecan_bad_geo_prediction_ids_2 =  [key for key  in map_irinotecan_geo_results.keys() if float(map_irinotecan_geo_results[key][1])/  float(map_irinotecan_geo_results[key][2]) >=  0.10   ]
printlist(list_CRC_irinotecan_bad_geo_prediction_ids_2 , 'list_CRC_irinotecan_bad_geo_prediction_ids_2')

# checking if th irinotecan aptietns are also present in the fluorouracil patients clinical file from fiamma

##missing_list = []
##for pat_id in map_irinotecan_geo_results.keys():
##    if pat_id not in  treatments_map.keys():
##        print ('missing irinotecan patient in fluorouracil clincial data :', pat_id)
##        missing_list.append(pat_id)
##printlist(missing_list, 'missing_list')


#(map_pfs_months, map_pfs_status, map_os_months , map_os_status,  map_dss_months,  map_dss_status) =  generate_maps()

##
##logging_string_7 = 'GEO fluorouracil pfs analysis'
##event_string_pfs =  '1:PROGRESSION'
##log_rank_analysis(list_fluorouracil_good_geo_prediction_ids, list_fluorouracil_bad_geo_prediction_ids, map_pfs_months, map_pfs_status, event_string_pfs, logging_string_7)
##ttest_ks_analysis(list_fluorouracil_good_geo_prediction_ids, list_fluorouracil_bad_geo_prediction_ids, map_pfs_months, map_pfs_status, event_string_pfs, logging_string_7)
##
##logging_string_8 = 'GEO fluorouracil os analysis'
##log_rank_analysis(list_fluorouracil_good_geo_prediction_ids, list_fluorouracil_bad_geo_prediction_ids, map_os_months, map_os_status, event_string_os, logging_string_8)
##ttest_ks_analysis(list_fluorouracil_good_cmap_prediction_ids, list_fluorouracil_bad_cmap_prediction_ids, map_os_months, map_os_status, event_string_os, logging_string_8)
##
##logging_string_9 = 'GEO fluorouracil dss analysis'
##log_rank_analysis(list_fluorouracil_good_geo_prediction_ids, list_ifluorouracil_bad_geo_prediction_ids, map_dss_months, map_dss_status, event_string_dss, logging_string_9)
##ttest_ks_analysis(list_fluorouracil_good_geo_prediction_ids, list_fluorouracil_bad_geo_prediction_ids, map_dss_months, map_dss_status, event_string_dss, logging_string_9)


##############################
# dati intabellati dai un sui dai dei pazienti BRCA fluorouracil

map_BRCA_fluorouracil_geo_results = {
'TCGA-GM-A4E0' : (2.118, 5, 95, 3, 17 , 'odds_ratio_log', 0.05),
'TCGA-GM-A3XN' : (2.681, 0, 69, 1, 15, 'pval_x_deltadeg', 0.001),
'TCGA-GM-A3XL' : (2.276, 6, 95, 4, 17, 'pval_x_drug_hits', 0.05),
'TCGA-GM-A3XG' : (1.702, 21,93, 5,16, 'odds_ratio_log', 0.05),
'TCGA-GM-A3NW' : (2.767, 3, 70, 3,15, 'deltadeg', 0.001),
'TCGA-GM-A2DN' : (2.859, 1, 97, 2, 17, 'c_score_x_deltadeg', 0.05),
'TCGA-GM-A2DM' : (2.247, 5, 96, 4,17, 'c_score_x_deltadeg', 0.05),
'TCGA-GM-A2DI' : (2.588, 2, 96, 3, 16, 'c_score_x_deltadeg', 0.05),
'TCGA-GM-A2DH' : (2.678, 4, 89, 3, 18, ' pval_x_deltadeg', 0.01),
'TCGA-GM-A2DF' : (2.732, 2, 70, 3,14, 'deltadeg', 0.001),
'TCGA-GM-A2DB' : (2.741, 2, 96, 3, 17, 'pval_x_deltadeg', 0.05),
'TCGA-GM-A2DA' : (2.706, 2, 67, 1,14, 'odds_ratio_log', 0.001),
'TCGA-E9-A3X8' : (2.671, 0, 95, 1, 16, 'deltadeg', 0.05),
'TCGA-E9-A3Q9' : (1.196, 41, 84, 12, 15, 'pval_x_drug_hits', 0.01),
'TCGA-E9-A2JS' : (2.711, 3, 76, 3, 15, 'pval_x_deltadeg', 0.001),
'TCGA-E9-A1NI' : (2.692, 1, 74, 2, 15, 'pval_x_deltadeg', 0.001),
'TCGA-E9-A1NE' : (2.792, 1, 96, 2, 16, 'deltadeg', 0.05 ),
'TCGA-C8-A8HR' : ( 2.683, 2, 97, 3, 17, 'deltadeg', 0.05 ),
'TCGA-C8-A27B' : (2.786, 1, 97, 2, 17,  'deltadeg', 0.05 ),
'TCGA-C8-A27A' : (2.776, 0, 75, 1, 15, 'deltadeg', 0.001),
'TCGA-C8-A26X' : ( 2.672, 1, 71,2, 15, 'pval_x_deltadeg', 0.001),
'TCGA-AC-A3YI' : (2.461, 7, 65, 5, 12, 'odds_ratio_log', 0.001),
'TCGA-AC-A2QI': (2.428, 1, 88, 2, 16, 'c_score', 0.01),
'TCGA-A8-A0AB': (2.731, 1, 92,2, 17, 'pval_x_deltadeg', 0.01),
'TCGA-A8-A07G' : (2.798, 1, 88, 2, 18, 'c_score_x_deltadeg', 0.01),
'TCGA-A8-A07C' : (2.708, 3, 96, 3, 16, 'pval_x_deltadeg', 0.05),
}


list_brca_fluorouracil_good_geo_prediction_ids_2 = [key for key  in map_BRCA_fluorouracil_geo_results.keys() if float(map_BRCA_fluorouracil_geo_results[key][1])/  float(map_BRCA_fluorouracil_geo_results[key][2]) < 0.10   ]
printlist(list_brca_fluorouracil_good_geo_prediction_ids_2, 'list_brca_fluorouracil_good_geo_prediction_ids_2')

list_brca_fluorouracil_bad_geo_prediction_ids_2 =  [key for key  in map_BRCA_fluorouracil_geo_results.keys() if float(map_BRCA_fluorouracil_geo_results[key][1])/  float(map_BRCA_fluorouracil_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_fluorouracil_bad_geo_prediction_ids_2 , 'list_brca_fluorouracil_bad_geo_prediction_ids_2')


sub_class_0 = [
'TCGA-GM-A3NW',
'TCGA-C8-A26X',
'TCGA-E9-A3X8',
'TCGA-A8-A07G',
'TCGA-A8-A0AB',
'TCGA-AC-A2QI',
'TCGA-AC-A3YI',
'TCGA-E9-A3Q9',
'TCGA-C8-A8HR',
'TCGA-E9-A1NE',
'TCGA-E9-A1NI',
'TCGA-E9-A2JS',
]

new_map_BRCA_fluorouracil_geo_results = {}
for x in sub_class_0:
    print(x in map_BRCA_fluorouracil_geo_results.keys())
    new_map_BRCA_fluorouracil_geo_results[x] = map_BRCA_fluorouracil_geo_results[x]

print('subclass 1 : ',  ['cyclophosphamide', 'fluorouracil', 'dox hyd'])
list_brca_fluorouracil_good_geo_prediction_ids_2 = [key for key  in sub_class_0 if float(map_BRCA_fluorouracil_geo_results[key][1])/  float(map_BRCA_fluorouracil_geo_results[key][2]) < 0.10   ]
printlist(list_brca_fluorouracil_good_geo_prediction_ids_2, 'list_brca_fluorouracil_good_geo_prediction_ids_2 - sub_class_0')

list_brca_fluorouracil_bad_geo_prediction_ids_2 =  [key for key  in sub_class_0 if float(map_BRCA_fluorouracil_geo_results[key][1])/  float(map_BRCA_fluorouracil_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_fluorouracil_bad_geo_prediction_ids_2 , 'list_brca_fluorouracil_bad_geo_prediction_ids_2 - sub_class_0')

#export as  new global variable to be set to use this data.
map_BRCA_fluorouracil_geo_results =  new_map_BRCA_fluorouracil_geo_results




#######################################
#

# dati intabellati dai un sui dai dei pazienti BRCA fluorouracil

map_BRCA_doxorubicin_geo_results = {
'TCGA-BH-A18H' : (2.654, 0, 95, 1,17,'pval_x_deltadeg', 0.05),
'TCGA-BH-A0DI':   (2.818, 0,96,1,18, 'pval_x_deltadeg', 0.05),
'TCGA-AC-A2FK':  (2.153, 2,88,1,18, 'log10', 0.01),
'TCGA-LL-A740' :  (2.868, 1,87,2,18, 'deltdeg', 0.01),
'TCGA-A7-A4SE' : (2.758, 0, 72, 1, 15, 'deltdeg', 0.001) ,
'TCGA-A2-A3XS' : (2.804,0, 89,1,18, 'deltdeg', 0.01),
'TCGA-E9-A6HE': (2.715,0,95,1,18, ' c_score_x_deltadeg', 0.05 ),
'TCGA-BH-A0B5': (2.672, 1,75,2,14, 'deltdeg', 0.001),
'TCGA-A2-A4RY': (1.561,6,92,2,17, 'odds_ratio_log', 0.05),
'TCGA-OL-A66I': (2.678,0,69,1,14, 'deltdeg', 0.001),
'TCGA-LL-A73Z': (2.888, 2,96,3,17,'deltdeg', 0.05),
'TCGA-A2-A3XX': (2.790,0,71,1,15, 'deltdeg', 0.001),
'TCGA-MS-A51U': (2.411,1,86,2,17, 'odds_ratio_log', 0.01),
'TCGA-A2-A3Y0': (2.677,0,98,1,19, 'pval_x_deltadeg', 0.05),
'TCGA-BH-A0BS': (2.733,1,97,2,18, 'deltdeg', 0.05),
'TCGA-E9-A1R2': (2.702, 1,91,2,17, 'pval_x_deltadeg_x_drug_hits ', 0.01),
'TCGA-AC-A2FF': (2.686, 0,69,1,17, 'deltdeg', 0.001),
'TCGA-AC-A3BB': (2.426,1,94,2,17, 'deltdeg', 0.05),
'TCGA-A2-A4S2': (2.529, 1,76,2,16, 'odds_ratio_log', 0.001),
'TCGA-BH-A0DG': (2.773,1,90,2,17, 'pval_x_deltadeg', 0.01),
'TCGA-AC-A2FE': (2.716,1,75,2,15, 'deltdeg', 0.001),
'TCGA-A2-A0EP': (2.691, 1,83,2,18, 'pval_x_deltadeg', 0.01),
'TCGA-BH-A0C3': (2.750,1,87,2,16,'pval_x_deltadeg', 0.01 ),
'TCGA-OL-A5RW': (2.724,0,96,1,17, 'deltdeg', 0.05),
'TCGA-A2-A3XY': (2.700,1,65,2,16, 'pval_x_deltadeg', 0.001),
'TCGA-A2-A3XW': (2.846,0,96,1,18, 'deltdeg', 0.05),
'TCGA-BH-A0BT': (2.700,0,97,1,18,'pval_x_deltadeg', 0.05),
'TCGA-LL-A7SZ': (2.775, 1,88,2,17,'pval_x_deltadeg_x_drug_hits', 0.01),
'TCGA-A2-A0CR': (2.639,0,93,1,17, 'deltdeg', 0.05),
'TCGA-BH-A0DV': (2.722,1,95,2,18, 'deltdeg', 0.05),
'TCGA-S3-A6ZH': (2.705,1,87,2,17, 'pval_x_deltadeg', 0.01),
'TCGA-A2-A0CK': (2.703,0,88,1,17, 'pval_x_deltadeg', 0.01),
'TCGA-BH-A0DT': (2.727,1,97,2,17, 'deltdeg', 0.05),
'TCGA-AC-A2BK': (2.644,1,87,2,17,'pval_x_deltadeg', 0.01)
}


list_brca_doxorubicin_good_geo_prediction_ids_2 = [key for key  in map_BRCA_doxorubicin_geo_results.keys() if float(map_BRCA_doxorubicin_geo_results[key][1])/  float(map_BRCA_doxorubicin_geo_results[key][2]) < 0.10   ]
printlist(list_brca_doxorubicin_good_geo_prediction_ids_2, 'list_brca_doxorubicin_good_geo_prediction_ids_2')

list_brca_doxorubicin_bad_geo_prediction_ids_2 =  [key for key  in map_BRCA_doxorubicin_geo_results.keys() if float(map_BRCA_doxorubicin_geo_results[key][1])/  float(map_BRCA_doxorubicin_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_doxorubicin_bad_geo_prediction_ids_2 , 'list_brca_doxorubicin_bad_geo_prediction_ids_2')

###################################
#
# dati intabellati dai un sui dai dei pazienti BRCA paclitaxel

exclusive_GEO_brca_paclitaxel_map = {
'TCGA-OL-A5S0': (2.786,28,82,12,16,  'pval_x_deltadeg', 0.01), # paclitaxel only
'TCGA-D8-A3Z6': (2.731, 2, 96,3, 16, 'deltdeg', 0.05),
'TCGA-A7-A6VY': (2.690, 1,69,2,14, 'pval_x_deltadeg', 0.001),
'TCGA-OL-A5RV': (2.129, 1,96,2,16, 'pval_x_deltadeg', 0.05),
'TCGA-C8-A1HM': (2.703, 9,76,1,15, 'deltdeg', 0.001),
'TCGA-BH-A0AZ': (2.714, 0,95,1,15, 'deltdeg', 0.05),
'TCGA-EW-A2FR': (2.791, 35,75,13,15, 'pval_x_deltadeg', 0.001), # paclitaxel only
'TCGA-S3-A6ZG': (2.707, 1,70,2,15, 'deltdeg', 0.001),
'TCGA-A7-A6VW': (2.747,1,85,2,14, 'deltdeg', 0.01 ),
'TCGA-EW-A1OY': (2.710, 24,89,12,16, 'pval_x_deltadeg_x_drug_hits', 0.01),  # paclitaxel only
'TCGA-A7-A4SA': (2.692, 1,89,2,15,  'pval_x_deltadeg', 0.01),
    }

shared_GEO_brca_paclitaxel = ['TCGA-C8-A27A', 'TCGA-LL-A740', 'TCGA-A2-A3XS', 'TCGA-C8-A27B', 'TCGA-GM-A2DF', 'TCGA-GM-A3XL', 'TCGA-AC-A2BK', 'TCGA-GM-A2DM',
                                                      'TCGA-A7-A4SE', 'TCGA-AC-A2FF', 'TCGA-GM-A3XG', 'TCGA-OL-A5RW', 'TCGA-A2-A3XW', 'TCGA-GM-A3XN', 'TCGA-OL-A66I', 'TCGA-BH-A0DI',
                                                      'TCGA-GM-A2DH', 'TCGA-A2-A0CK', 'TCGA-A2-A3XY', 'TCGA-BH-A0B5', 'TCGA-GM-A2DB', 'TCGA-A2-A0EP', 'TCGA-E9-A6HE', 'TCGA-A2-A4S2',
                                                      'TCGA-MS-A51U', 'TCGA-LL-A7SZ', 'TCGA-A2-A0CR', 'TCGA-A2-A3XX', 'TCGA-GM-A2DN']


map_BRCA_paclitaxel_geo_results = {}
for pat in exclusive_GEO_brca_paclitaxel_map.keys() :
    map_BRCA_paclitaxel_geo_results[pat] = exclusive_GEO_brca_paclitaxel_map[pat]
for pat in shared_GEO_brca_paclitaxel :
    if pat in  map_BRCA_paclitaxel_geo_results.keys() :
        print('Duplicate ', pat)
        assert False
    elif pat in map_BRCA_fluorouracil_geo_results.keys():
        map_BRCA_paclitaxel_geo_results[pat] = map_BRCA_fluorouracil_geo_results[pat]
    elif pat in map_BRCA_doxorubicin_geo_results.keys() :
        map_BRCA_paclitaxel_geo_results[pat] = map_BRCA_doxorubicin_geo_results[pat]
    else :
        print('Missing data for ', pat)

list_brca_paclitaxel_good_geo_prediction_ids_2 = [key for key  in map_BRCA_paclitaxel_geo_results.keys() if float(map_BRCA_paclitaxel_geo_results[key][1])/  float(map_BRCA_paclitaxel_geo_results[key][2]) < 0.10   ]
printlist(list_brca_paclitaxel_good_geo_prediction_ids_2, 'list_brca_paclitaxel_good_geo_prediction_ids_2')

list_brca_paclitaxel_bad_geo_prediction_ids_2 =  [key for key  in map_BRCA_paclitaxel_geo_results.keys() if float(map_BRCA_paclitaxel_geo_results[key][1])/  float(map_BRCA_paclitaxel_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_paclitaxel_bad_geo_prediction_ids_2 , 'list_brca_paclitaxel_bad_geo_prediction_ids_2')


# subclass x paracitaxel :
#['cyclophosphamide', 'paclitaxel', 'dox hyd']

sub_class_1 =[
'TCGA-E9-A6HE',
'TCGA-BH-A0B5',
'TCGA-BH-A0DI',
'TCGA-AC-A2FF',
'TCGA-A2-A0CR',
'TCGA-LL-A740',
'TCGA-A2-A4S2',
'TCGA-A2-A3XY',
'TCGA-A7-A4SE',
'TCGA-A2-A3XX',
'TCGA-A2-A3XW',
'TCGA-OL-A66I',
'TCGA-OL-A5RW',
'TCGA-A2-A3XS',
'TCGA-MS-A51U',
'TCGA-LL-A7SZ',
'TCGA-A2-A0EP',
'TCGA-AC-A2BK',
]


sub_class_4 = sub_class_1 + ['TCGA-A7-A6VW', 'TCGA-OL-A5RV', 'TCGA-D8-A3Z6', 'TCGA-BH-A0AZ', 'TCGA-A7-A4SA', 'TCGA-A7-A6VY', 'TCGA-LL-A73Z',
#                             'TCGA-B6-A402',
                             'TCGA-S3-A6ZG', 'TCGA-AC-A2FE', 'TCGA-C8-A1HM']
# chek for duplicates
assert len(sub_class_4) == len(set(sub_class_4))



def new_sub_data_paclitaxel (sub_class_list, subclass_string):

    new_map_BRCA_paclitaxel_geo_results = {}
    for x in sub_class_list:
        if x in map_BRCA_paclitaxel_geo_results.keys() :
            print(x in map_BRCA_paclitaxel_geo_results.keys())
            new_map_BRCA_paclitaxel_geo_results[x] = map_BRCA_paclitaxel_geo_results[x]
        elif x in map_BRCA_doxorubicin_geo_results.keys() :
            new_map_BRCA_paclitaxel_geo_results[x] =  map_BRCA_doxorubicin_geo_results[x]

    print(subclass_string)
    good_geo_prediction= [key for key  in sub_class_list if float(new_map_BRCA_paclitaxel_geo_results[key][1])/  float(new_map_BRCA_paclitaxel_geo_results[key][2]) < 0.10   ]
    printlist(good_geo_prediction, 'good_geo_prediction :' + subclass_string)

    bad_geo_prediction=  [key for key  in sub_class_list if float(new_map_BRCA_paclitaxel_geo_results[key][1])/  float(new_map_BRCA_paclitaxel_geo_results[key][2]) >=  0.10   ]
    printlist(bad_geo_prediction , 'bad_geo_prediction :2' +  subclass_string)

    #export as  new global variable to be set to use this data.
    return(new_map_BRCA_paclitaxel_geo_results, good_geo_prediction, bad_geo_prediction)


(map_BRCA_paclitaxel_geo_results, list_brca_paclitaxel_good_geo_prediction_ids_2 , list_brca_paclitaxel_bad_geo_prediction_ids_2 ) = new_sub_data_paclitaxel (sub_class_4, 'paclitaxel sub-population 4: with Cyclophosphamide , without Fluorouracil and docetaxel')




# check inclusion

new_map = {}
for x in sub_class_1:
    print(x, 'present', x in map_BRCA_paclitaxel_geo_results.keys())
    new_map[x] = map_BRCA_paclitaxel_geo_results[x]

print('subclass 1 : ',  ['cyclophosphamide', 'paclitaxel', 'dox hyd'])
list_brca_Doxorubicin_Hydrochloride_E_1_good_geo_prediction_ids_2 = [key for key  in sub_class_1 if float(map_BRCA_paclitaxel_geo_results[key][1])/  float(map_BRCA_paclitaxel_geo_results[key][2]) < 0.10   ]
printlist(list_brca_Doxorubicin_Hydrochloride_E_1_good_geo_prediction_ids_2 , 'llist_brca_Doxorubicin_Hydrochloride_E_1_good_geo_prediction_ids_2')

list_brca_Doxorubicin_Hydrochloride_E_1_bad_geo_prediction_ids_2 =  [key for key  in sub_class_1 if float(map_BRCA_paclitaxel_geo_results[key][1])/  float(map_BRCA_paclitaxel_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_Doxorubicin_Hydrochloride_E_1_bad_geo_prediction_ids_2 , 'list_brca_Doxorubicin_Hydrochloride_E_1_bad_geo_prediction_ids_2')

#export as  new global variable to be set to use this data.
map_BRCA_Doxorubicin_Hydrochloride_E_1_geo_results =  new_map

###################

###################################
#
# dati intabellati dai un sui dai dei pazienti BRCA docetaxel

exclusive_GEO_brca_docetaxel_map = {
 'TCGA-A7-A4SD' : (2.786, 2, 97, 3, 17, 'deltadeg', 0.05) ,
 'TCGA-AC-A3TM' : (2.879, 3, 85, 4, 16, 'c_score_x_deltadeg', 0.01) ,
 'TCGA-A7-A4SF' : (2.710, 7, 90, 6, 16, 'pval_x_deltadeg', 0.01) ,
 'TCGA-OL-A66O' : (2.731, 3, 97, 3, 17, 'pval_x_deltadeg', 0.05) ,
 'TCGA-GM-A2DC' : (2.137, 10,96, 4,17, 'odds_ratio_log', 0.05) ,
 'TCGA-LL-A5YO' : (2.768, 2, 70, 3, 14, 'deltadeg', 0.001) ,
 'TCGA-A7-A6VX' : (2.777, 2, 77, 3, 15, 'deltadeg', 0.001),
 'TCGA-A2-A3XV' : (2.864, 9, 97, 8, 17, 'pval_x_deltadeg', 0.05),
 'TCGA-EW-A424' : (2.685, 2, 73, 3, 15, 'pval_x_deltadeg', 0.001),
 'TCGA-BH-A0DD' : (2.688, 12, 70, 10, 14, 'pval_x_deltadeg', 0.001),
 'TCGA-A8-A08O' : (2.593, 79, 92, 15, 16, 'pval_x_deltadeg ', 0.01),
 'TCGA-A7-A4SB' : (2.670, 3, 89, 3, 16, 'pval_x_deltadeg ', 0.01),
 'TCGA-EW-A3U0' :(2.670, 1, 96, 2 , 17, 'pval_x_deltadeg', 0.05),
 'TCGA-GM-A5PV' : (2.530, 10, 69, 7, 13, 'odds_ratio_log', 0.001),
 'TCGA-LL-A441' : (2.763, 2, 70, 3, 14 ,  'pval_x_deltadeg', 0.001),
 'TCGA-AC-A2FB' : (2.805, 84, 98, 15, 18, 'deltadeg', 0.05),
 'TCGA-BH-A18F' : (2.717, 2, 87, 3, 15, 'deltadeg', 0.01),
 'TCGA-BH-A0HY' : (2.879, 8, 98, 7, 19,  'deltadeg', 0.05),
 'TCGA-A8-A097' : (2.754, 2, 91, 14, 16, 'deltadeg', 0.01),
 'TCGA-BH-A18I' : (2.684, 2, 98, 3, 17, 'deltadeg', 0.05),
 'TCGA-E2-A573' : (2.713,2, 90,3, 15,'odds_ratio_log', 0.01 ),
 'TCGA-A7-A5ZW' : (2.232, 9, 68,5, 11, 'odds_ratio_log', 0.001),
 'TCGA-A7-A6VV' : (2.715,2, 87, 3, 15, 'deltadeg', 0.01 ),
 'TCGA-BH-A8G0' : (1.725, 12, 95, 4, 16,  'c_score_x_deltadeg', 0.05),
 'TCGA-LD-A7W5' : (2.698, 94, 96, 16, 16,  'deltadeg', 0.05),
 'TCGA-XX-A899' : (2.645, 3, 94, 3, 16,  'pval_x_deltadeg', 0.05),
 'TCGA-S3-A6ZF' : (2.786, 2, 74, 3, 15, 'deltadeg', 0.001),
 'TCGA-A2-A3XT' : (2.681, 1, 73, 2, 15,  'pval_x_deltadeg', 0.001),
 'TCGA-A2-A3XZ' : (2.778, 8, 96, 7, 17, 'deltadeg', 0.05),
 'TCGA-OL-A66P' : (2.804, 3, 96, 3, 18, 'pval_x_deltadeg_x_drug_hits', 0.05),
 'TCGA-BH-A202' : (2.799, 6, 78,5, 15, 'deltadeg', 0.001),
 'TCGA-OL-A6VO' : (2.710, 2, 86, 3, 14, 'deltadeg', 0.01),
 }

shared_GEO_brca_docetaxel = ['TCGA-GM-A2DI', 'TCGA-E9-A1R2', 'TCGA-GM-A4E0', 'TCGA-A2-A3XY', 'TCGA-BH-A0DV', 'TCGA-AC-A3BB', 'TCGA-AC-A2FK', 'TCGA-A2-A3Y0', 'TCGA-A2-A3XS', 'TCGA-BH-A18H', 'TCGA-A2-A0CK', 'TCGA-A8-A07C', 'TCGA-A2-A4RY', 'TCGA-GM-A2DA', 'TCGA-BH-A0DG', 'TCGA-BH-A0DT']



map_BRCA_docetaxel_geo_results = {}
for pat in exclusive_GEO_brca_docetaxel_map.keys() :
    map_BRCA_docetaxel_geo_results[pat] = exclusive_GEO_brca_docetaxel_map[pat]
for pat in shared_GEO_brca_docetaxel :
    if pat in  map_BRCA_docetaxel_geo_results.keys() :
        print('Duplicate ', pat)
        assert False
    elif pat in map_BRCA_fluorouracil_geo_results.keys():
        map_BRCA_docetaxel_geo_results[pat] = map_BRCA_fluorouracil_geo_results[pat]
    elif pat in map_BRCA_doxorubicin_geo_results.keys() :
        map_BRCA_docetaxel_geo_results[pat] = map_BRCA_doxorubicin_geo_results[pat]
    elif pat in map_BRCA_paclitaxel_geo_results.keys() :
        map_BRCA_docetaxel_geo_results[pat] = map_BRCA_paclitaxel_geo_results[pat]
    else :
        print('Missing data for ', pat)

list_brca_docetaxel_good_geo_prediction_ids_2 = [key for key  in map_BRCA_docetaxel_geo_results.keys() if float(map_BRCA_docetaxel_geo_results[key][1])/  float(map_BRCA_docetaxel_geo_results[key][2]) < 0.10   ]
printlist(list_brca_docetaxel_good_geo_prediction_ids_2, 'list_brca_docetaxel_good_geo_prediction_ids_2')

list_brca_docetaxel_bad_geo_prediction_ids_2 =  [key for key  in map_BRCA_docetaxel_geo_results.keys() if float(map_BRCA_docetaxel_geo_results[key][1])/  float(map_BRCA_docetaxel_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_docetaxel_bad_geo_prediction_ids_2 , 'list_brca_docetaxel_bad_geo_prediction_ids_2')


sub_class_2 = [
'TCGA-GM-A2DC',
'TCGA-EW-A424',
'TCGA-S3-A6ZF',
'TCGA-LL-A5YO',
'TCGA-OL-A6VO',
'TCGA-OL-A66P',
'TCGA-OL-A66O',
'TCGA-GM-A5PV',
'TCGA-LL-A441',
'TCGA-EW-A3U0',
'TCGA-AC-A3TM',
'TCGA-A7-A4SB',
'TCGA-A7-A4SD',
'TCGA-A7-A5ZW',
'TCGA-A7-A6VV',
'TCGA-A7-A6VX',
'TCGA-A8-A097',
'TCGA-BH-A18F',
'TCGA-BH-A18I',
'TCGA-BH-A8G0',
'TCGA-E2-A573',
'TCGA-XX-A899',
]


sub_class_3 = sub_class_2 + [  'TCGA-BH-A18H', 'TCGA-A2-A0CK', 'TCGA-BH-A0DG', 'TCGA-BH-A0DV', 'TCGA-A2-A3XY', 'TCGA-AC-A2FK', 'TCGA-A2-A3Y0', 'TCGA-A2-A4RY', 'TCGA-A2-A3XT', 'TCGA-A2-A3XS', 'TCGA-AC-A3BB', 'TCGA-BH-A0DT', 'TCGA-E9-A1R2']
# chek for duplicates
assert len(sub_class_3) == len(set(sub_class_3))

def new_sub_data_docetaxel (sub_class_list, subclass_string):

    new_map_BRCA_docetaxel_geo_results = {}
    for x in sub_class_list:
        print(x in map_BRCA_docetaxel_geo_results.keys())
        new_map_BRCA_docetaxel_geo_results[x] = map_BRCA_docetaxel_geo_results[x]

    print(subclass_string)
    good_geo_prediction= [key for key  in sub_class_list if float(map_BRCA_docetaxel_geo_results[key][1])/  float(map_BRCA_docetaxel_geo_results[key][2]) < 0.10   ]
    printlist(good_geo_prediction, 'good_geo_prediction :' + subclass_string)

    bad_geo_prediction=  [key for key  in sub_class_list if float(map_BRCA_docetaxel_geo_results[key][1])/  float(map_BRCA_docetaxel_geo_results[key][2]) >=  0.10   ]
    printlist(bad_geo_prediction , 'bad_geo_prediction :2' +  subclass_string)

    #export as  new global variable to be set to use this data.
    return(new_map_BRCA_docetaxel_geo_results, good_geo_prediction, bad_geo_prediction)


(map_BRCA_docetaxel_geo_results,list_brca_docetaxel_good_geo_prediction_ids_2 , list_brca_docetaxel_bad_geo_prediction_ids_2 ) = new_sub_data_docetaxel (sub_class_3, 'Docetaxel sub-population 3: with Cyclophosphamide , Docetaxel , without Fluorouracil')




new_map= {}
for x in sub_class_2:
    print(x, 'present', x in map_BRCA_docetaxel_geo_results.keys())
    new_map[x] = map_BRCA_docetaxel_geo_results[x]

print('subclass 2 : ',  ['cyclophosphamide', 'docetaxel'])
list_brca_docetaxel_D_2_good_geo_prediction_ids_2 = [key for key  in sub_class_2 if float(map_BRCA_docetaxel_geo_results[key][1])/  float(map_BRCA_docetaxel_geo_results[key][2]) < 0.10   ]
printlist(list_brca_docetaxel_D_2_good_geo_prediction_ids_2, 'list_brca_docetaxel_D_2_good_geo_prediction_ids_2 - sub_class_2')

list_brca_docetaxel_D_2_bad_geo_prediction_ids_2 =  [key for key  in sub_class_2 if float(map_BRCA_docetaxel_geo_results[key][1])/  float(map_BRCA_docetaxel_geo_results[key][2]) >=  0.10   ]
printlist(list_brca_docetaxel_D_2_bad_geo_prediction_ids_2 , 'list_brca_docetaxel_D_2_bad_geo_prediction_ids_2 - sub_class_2')

#export as  new global variable to be set to use this data.
map_BRCA_docetaxel_D_2_geo_results =  new_map




######################
##local_out = patients_by_groups(['Fluorouracil'],['Docetaxel', 'Paclitaxel'])
##venndiagramsets('sub_class_0', set(sub_class_0), 'local search' , set(local_out)) 
##local_out = patients_by_groups(['Fluorouracil'],['Docetaxel', 'Paclitaxel'])
####################
##local_out = patients_by_groups(['Cyclophosphamide', 'Docetaxel'], ['Fluorouracil', 'Paclitaxel'])
##venndiagramsets( 'sub_class_2', set(sub_class_2),  'local search' , set(local_out))
###############
##local_out = patients_by_groups(['Cyclophosphamide', 'Paclitaxel'], ['Fluorouracil', 'Docetaxel'])
##venndiagramsets( 'sub_class_1', set(sub_class_1),  'local search' , set(local_out))




######################################################
# uploading data for the STAD patients.

def STAD_gen_drugs_to_check_maps(target_drugs) :

# To be changed for each drug target on BRCA data 
#    target_drugs = [ 'Fluorouracil']
#    target_drugs =['Doxorubicin Hydrochloride', 'Doxorubicin']
 #   target_drugs =[ 'Paclitaxel', ]    
    verbose = False
 #   verbose = True

    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'

    CLINICAL_DATA_FILE_NAME = 'STAD_fluorouracil_sample_and_clinical_data_4.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, ';', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index('cases.submitter_id')
    treatment_type_idx =   header.index( 'treatments.treatment_type')
    treatment_idx  =   header.index('treatments.therapeutic_agents')
    response_idx = header.index('treatments.treatment_outcome')

    start_time_idx = header.index( 'treatments.days_to_treatment_start')
    end_time_idex = header.index( 'treatments.days_to_treatment_end')

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
        if treatment_type in ['Hormone Therapy',   'Chemotherapy' ] and treatment != 'Not Reported' :
            if verbose :
                print (pat_id, treatment, treatment_type, response)
            pat_set.add(pat_id)
            drug_set.add(treatment)
            resp_set.add(response)
            increment_lists(treatments_map, pat_id, treatment)
            if response not in ['"-"', '"--"' , "'--",  'Unknown', 'Treatment Ongoing'] :
                increment_lists(resp_map, pat_id, (treatment,response))
                if treatment in target_drugs :
                    increment_lists(reduced_response_map, pat_id, (response, duration))

    print('pat set', len(pat_set), pat_set)
    print('drug set' , len(drug_set) , drug_set)
    print('responses', len(resp_set), resp_set)


##    f= open('fiamma_list.tex')
##    fiamma_list = f.read().split()
##
###    printlist(fiamma_list, 'fiamma_list')
##
##    venndiagramsets('fiamma_list', set(fiamma_list), 'clinical_fluorouracil_data', pat_set)
      
    printmap(treatments_map, 'treatments_map')
    printmap(resp_map, 'responses_map')
    printmap(reduced_response_map, 'reduced responses_map')

    drugs_to_check_lable = {}
    for key in treatments_map.keys():
        drugs_to_check_lable[key] = 'BRCA_drug_used_for_'+ key

    return(treatments_map, drugs_to_check_lable, reduced_response_map)



####################
# run       for  STAD

#target_drugs = ['Fluorouracil']
#(treatments_map, drugs_to_check_lable, reduced_response_map) = STAD_gen_drugs_to_check_maps(target_drugs)



###################
# STAD data




map_stad_fluorouracil_geo_results = {
 'TCGA-3M-AB46' : (2.107, 5, 83, 3, 14, 'deltadeg', 0.01),

 'TCGA-BR-6457' :  (1.131,  20, 91, 6, 14, 'odds_ratio_log', 0.05),
 'TCGA-BR-6458'  : (2.017,  11  , 85 , 6  , 14  , 'odds_ratio_log' ,  0.01 ),
 'TCGA-BR-6706'  : ( 1.774 ,   12,  83 ,5,  14 ,  'c_score_x_deltadeg'  ,   0.01),
 'TCGA-BR-6707'  : (  2.000    ,  8 ,  70  ,  5 ,  14 ,  'odds_ratio_log' ,  0.001 ),
 'TCGA-BR-6710'  : ( 0.869   , 16 ,  95,   8,  14 ,  'odds_ratio_log', 0.05 ),
 'TCGA-BR-6802': ( 2.269, 7 , 72 , 5  , 14 ,  'deltadeg',  0.001),
 
 'TCGA-BR-8373': ( 2.259, 62 , 94 , 9 , 14 ,  'c_score_x_deltadeg'  ,   0.05 ),
 'TCGA-BR-8381': (1.995 , 19 , 95 , 4 , 14 ,  'c_score_x_deltadeg'  ,   0.05 ),
 'TCGA-BR-8590': (1.369 , 12 , 91 , 6 ,  14,  'odds_ratio_log', 0.05 ),
 'TCGA-BR-8592': (1.198 ,  21,  86  , 8 , 14  ,   'odds_ratio_log', 0.05),
 'TCGA-BR-8690': ( 0.961  ,  61,  83  ,  11 , 13,  'odds_ratio_log' ,  0.01),
 'TCGA-BR-A44T': ( 1.000,  22  ,  95,  6  , 14,  'odds_ratio_log', 0.05),

 
 'TCGA-D7-6520': ( 1.760 ,  1,  94,  1, 15 ,  'c_score' ,  0.05),
 'TCGA-F1-A448': (1.985 ,  27,  70,  9, 14 ,  'pval_x_deltadeg',  0.001),
 'TCGA-FP-7829': (2.177 ,  11,  73,  7,  14,  'deltadeg',  0.001) ,
 'TCGA-FP-7998': ( 1.947,  10, 65 , 5 , 14 ,  'deltadeg',  0.001 ),
 'TCGA-FP-A4BE': (2.269 ,  7,  68,  5  , 14 ,  'deltadeg',  0.001),
 'TCGA-HJ-7597': (1.298 ,  65 ,  83,  11,  13,  'odds_ratio_log' ,  0.01 ),
 'TCGA-MX-A666': ( 2.100  ,  6,  69, 4 , 14  ,  'deltadeg',  0.001),
 
 'TCGA-VQ-A8E0': ( 2.263,  5 , 70 ,  3,  14,  'deltadeg',  0.001 ),
 'TCGA-VQ-A8PD': ( 1.884,  53, 77 , 13 ,  14,  'pval_x_deltadeg',  0.01),
 'TCGA-VQ-A8PP': ( 1.546,  18,72  ,  9,   13 ,  'deltadeg',  0.001  ),
 'TCGA-VQ-A91D': ( 1.606, 6 , 87 , 4,14  ,  'deltadeg',  0.01),
 'TCGA-VQ-A91E': ( 1.714  ,  7,  70 ,  5,  14 ,  'deltadeg',  0.001),
 'TCGA-VQ-A91K': (2.099 , 6 , 95 , 4 , 14 ,  'pval_x_deltadeg',  0.05 ),
 'TCGA-VQ-A91N': (2.237 ,  10,  63,  6,  14 ,  'deltadeg',  0.001 ),
 'TCGA-VQ-A91Q': ( 1.312, 19 , 95 , 7 , 14 ,'odds_ratio_log', 0.05 ),
 'TCGA-VQ-A91S': (1.656 ,  11, 74 ,  7,  14  ,  'deltadeg',  0.001),
 'TCGA-VQ-A91V': ( 2.014, 6 ,  70,  4 ,  14,    'deltadeg',  0.001 ),
 'TCGA-VQ-A91W': (2.136 , 11 ,  69,  6,  14, ' deltadeg',  0.001  ),
 'TCGA-VQ-A91Y': ( 2.175  ,  48,  79,  12,  14 ,  'c_score' ,  0.01 ),
 'TCGA-VQ-A91Z': ( 2.090,  5,87  , 3 ,14 , 'deltadeg',  0.01  ),
 'TCGA-VQ-A922': (2.307 , 8 ,  71,  5,  14,  'deltadeg',  0.001),
 'TCGA-VQ-A928': (1.688 ,  9  , 94 ,  5,   14,  'c_score_x_deltadeg'  ,   0.05),
 'TCGA-VQ-A94O': (1.652 ,  8 ,  88 ,  6  , 14 ,  'pval_x_deltadeg',  0.01),
 'TCGA-VQ-A94R': ( 1.863, 24 , 64 , 9 , 14 , 'odds_ratio_log' ,  0.001  ),
 'TCGA-VQ-AA64': ( 2.478,  16,  69,  6  , 13  ,  'c_score_x_deltadeg'  ,   0.001),
 'TCGA-VQ-AA68': (2.349 ,  6, 66 , 4 ,  12,   'deltadeg',  0.001 ),
 'TCGA-VQ-AA6B': (1.668 , 9 ,90  , 5 ,  14,'c_score_x_deltadeg'  ,   0.001 ),
 'TCGA-VQ-AA6D': (2.202 ,  9, 72  ,  6  ,  14,  'deltadeg',  0.001),
 'TCGA-VQ-AA6F': (1.595 ,  7, 70 ,  5,  14,  'deltadeg',  0.001),
 'TCGA-VQ-AA6G': (2.071 ,  7,89  ,  5,  14,  'deltadeg', 0.01),
 'TCGA-VQ-AA6K': ( 2.246, 6 ,83  ,  4,  14,  'deltadeg', 0.01),
    }


list_stad_fluorouracil_good_geo_prediction_ids_2 = [key for key  in map_stad_fluorouracil_geo_results.keys() if float(map_stad_fluorouracil_geo_results[key][1])/  float(map_stad_fluorouracil_geo_results[key][2]) < 0.10   ]
printlist(list_stad_fluorouracil_good_geo_prediction_ids_2, 'list_stad_fluorouracil_good_geo_prediction_ids_2')

list_stad_fluorouracil_bad_geo_prediction_ids_2 =  [key for key  in map_stad_fluorouracil_geo_results.keys() if float(map_stad_fluorouracil_geo_results[key][1])/  float(map_stad_fluorouracil_geo_results[key][2]) >=  0.10   ]
printlist(list_stad_fluorouracil_bad_geo_prediction_ids_2 , 'list_stad_fluorouracil_bad_geo_prediction_ids_2')

    
######################################
# STAD subclasses

# STAD fluorouracil only as drug

sub_STAD_fluorouracil_only = [
'TCGA-BR-6710',
'TCGA-BR-6457',
'TCGA-3M-AB46',
'TCGA-BR-A44T',
'TCGA-FP-7829',
'TCGA-FP-A4BE',
'TCGA-MX-A666',
'TCGA-VQ-A8E0',
'TCGA-VQ-A8PD',
'TCGA-VQ-A8PP',
'TCGA-VQ-A91D',
'TCGA-VQ-A91E',
'TCGA-VQ-A91K',
'TCGA-VQ-A91N',
'TCGA-VQ-A91Q',
'TCGA-VQ-A91S',
'TCGA-VQ-A91W',
'TCGA-VQ-A91Y',
'TCGA-VQ-A922',
'TCGA-VQ-A928',
'TCGA-VQ-A94O',
'TCGA-VQ-A94R',
'TCGA-VQ-AA64',
'TCGA-VQ-AA68',
'TCGA-VQ-AA6B',
'TCGA-VQ-AA6D',
'TCGA-VQ-AA6F',
'TCGA-VQ-AA6G',
'TCGA-VQ-AA6K',
    ]

def new_sub_data_STAD (sub_class_list, subclass_string):

    new_map_STAD_geo_results = {}
    for x in sub_class_list:
        print(x in map_stad_fluorouracil_geo_results.keys())
        new_map_STAD_geo_results[x] = map_stad_fluorouracil_geo_results[x]

    print(subclass_string)
    good_geo_prediction= [key for key  in sub_class_list if float(new_map_STAD_geo_results[key][1])/  float(new_map_STAD_geo_results[key][2]) < 0.10   ]
    printlist(good_geo_prediction, 'good_geo_prediction :' + subclass_string)

    bad_geo_prediction=  [key for key  in sub_class_list if float(new_map_STAD_geo_results[key][1])/  float(new_map_STAD_geo_results[key][2]) >=  0.10   ]
    printlist(bad_geo_prediction , 'bad_geo_prediction :2' +  subclass_string)

    #export as  new global variable to be set to use this data.
    return(new_map_STAD_geo_results, good_geo_prediction, bad_geo_prediction)


(map_stad_fluorouracil_only_geo_results,list_stad_fluorouracil_only_good_geo_prediction_ids_2 , list_stad_fluorouracil_only_bad_geo_prediction_ids_2 ) = new_sub_data_STAD (sub_STAD_fluorouracil_only, 'STAD subpopulation with Fluorouracil only')



