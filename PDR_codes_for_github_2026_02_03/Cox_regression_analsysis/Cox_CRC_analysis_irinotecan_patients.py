# Version for cox CRC irinitecan patient data data anysis 
# Started 2025-11-10


import os
import sys

###########################################
# Accessing functions in siblig directories

path = os.path.join("..", "")
print("Relative path exists", os.path.isdir(path))
abspath = os.path.abspath(path)
print("Absolute path exists", os.path.isdir(abspath))
sys.path.append(abspath)
print("Including in sys path to ", abspath)



from file_handling_util import testexistencedirectory, testexistencefile
from new_loader_routines import load_data
from generic_util import venndiagramsets , printmap, printlist

# threshold at 10\% for good/bad
# removed import embedded the data
#from clinical_data_analysis2 import map_irinotecan_geo_results , list_irinotecan_good_geo_prediction_ids_2, list_irinotecan_bad_geo_prediction_ids_2 



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










CLINICAL_DATA_PATH_NAME = '../../Single_patients_crc_therapy'
testexistencedirectory(CLINICAL_DATA_PATH_NAME)

################
# utilities

def map_update(inmap, key, value):
    if key not in inmap.keys() :
        inmap[key] = [value]
    else :
        inmap[key].append(value)


# compute duration of the therapy 
def compute_duration(start_time, end_time ):

#    print ('start, end', start_time, end_time)

    if start_time in ['' , "'--"] or end_time in ['' , "'--"] :
        return None
    start = int(start_time)
    end = int(end_time)
    assert end >= start
    return end - start


# nuova formulazione
prediction_class_irinotecan_dict_geo= {}
for key in list_CRC_irinotecan_good_geo_prediction_ids_2 :
    prediction_class_irinotecan_dict_geo[key] = 0
for key in list_CRC_irinotecan_bad_geo_prediction_ids_2 :
    prediction_class_irinotecan_dict_geo[key] = 1
    

printmap(prediction_class_irinotecan_dict_geo, 'prediction_class_dict_geo')

all_ids = list_CRC_irinotecan_good_geo_prediction_ids_2 +  list_CRC_irinotecan_bad_geo_prediction_ids_2



def get_result_map () :

    all_res = []
    acc_no_file = []
    acc_no_measue_of_response = []
    acc_no_therapy_ongoing = []
    map_response ={}

    for ids in all_ids:
        filename = ids + '_table.csv'
        res = testexistencefile(CLINICAL_DATA_PATH_NAME, filename)
        all_res.append(res)

        if not res :
            acc_no_file.append(filename)
            
        if res :
           
            data = load_data(CLINICAL_DATA_PATH_NAME, filename, '\t', 'verbose')
            acc = data[0]
            header = data[1]
            print(header)
            if 'MEASURE_OF_RESPONSE' not in header :
                print('@@ MEASURE_OF_RESPONSE missing in {}'.format(ids))
                acc_no_measue_of_response.append(ids)
                continue
            if 'THERAPY_ONGOING' not in header :
                print('@@ THERAPY_ONGOING missing in {}'.format(ids))
                acc_no_therapy_ongoing.append(ids)
                continue

            agent_index = header.index('AGENT')
            response_index =   header.index('MEASURE_OF_RESPONSE')
            ongoing_index = header.index('THERAPY_ONGOING')
            start_idx = header.index('START_DATE')
            stop_idx = header.index('STOP_DATE')

            
            for x in acc:
    #            print ('@@ ongoing {}, Agent {}, response{}'.format(x[ongoing_index], x[agent_index], x[response_index ]))
                if x[ongoing_index] == 'No' and x[agent_index] in ['Irinotecan', 'Fluorouracil']:
                    print ('@@@@ Selected : ongoing {}, Agent {}, response{}'.format(x[ongoing_index], x[agent_index], x[response_index ]))
                    if x[response_index ] != '' :
                        print('time data : ' , x[start_idx], x[stop_idx] )
                        duration = compute_duration( x[start_idx], x[stop_idx])
                        map_update(map_response, ids, (x[response_index ], duration))

                    


    print( 'all res', all_res)
    print('Exist : ', len([x for x in all_res if x ==True]))
    print('Not exist :', len([x for x in all_res if x == False]))
    print('Missing measure or response ', len(acc_no_measue_of_response), acc_no_measue_of_response)
    print('Missing therapy ongoing ', len(acc_no_therapy_ongoing), acc_no_therapy_ongoing)
    print('Missing file : ' , len(acc_no_file), acc_no_file)

    printmap(map_response, 'map_response' )

    # sort by duration

    for key in map_response.keys() :
        map_response[key]. sort(key = lambda y : y[1])

    printmap(map_response, 'map_response sorted by duration (increasing)' )

    return    map_response


# generate the repsonse map

map_irinotecan_response = get_result_map ()


#event_coding_lable = 'PD1_all0'
#event_coding_lable = 'PD2_CRPR1_all0'
#event_coding_lable = 'PDSD1_all0'
event_coding_lable = 'CRPR1_all0'

try:
    event_coding_lable
except NameError:
    print("event_coding_lable is not defined: please edit the script uncommenting one of the possible values!")
else:
    print("event_coding_lable is  defined.")

print('Event coding : ', event_coding_lable)

if event_coding_lable == 'PD1_all0' :
    # progressive disease = 1 all other 0
    print('Coding progressive disease=1 all other case 0 ')
    response_remapping =  {'partial response': 0, 'complete response': 0, 'stable disease':0, 'clinical progressive disease':1}

elif event_coding_lable == 'PDSD1_all0' :
    # progressive disease = 1 all other 0
    print('Coding progressive disease + stable disease =1 all other case 0 ')
    response_remapping =  {'partial response': 0, 'complete response': 0, 'stable disease':1, 'clinical progressive disease':1}


elif event_coding_lable == 'CR1_all0' :
    #compete response 1 all otehr 0
    print('Coding complete respons=1 all other case 0 ')
    response_remapping =  {'partial response': 0, 'complete response': 1, 'stable disease':0, 'clinical progressive disease':0}

elif event_coding_lable == 'CRPR1_all0' :
    #compete response + partial reponse 1 all otehr 0
    print('Codining Progressive + stable disease = 0, complete response+partial repsonse =1')
    response_remapping =  {'partial response': 1, 'complete response': 1, 'stable disease':0, 'clinical progressive disease':0}

elif event_coding_lable == 'PD2_CRPR1_all0' :
    #compete response + partial reponse 1 , progressive disease 2 , all otehr 0
    print('Codining  stable disease = 0, progeressive disease 2, complete response+partial repsonse =1')
    response_remapping =  {'partial response': 1, 'complete response': 1, 'stable disease':0, 'clinical progressive disease':2}
else :
    print('Error in event labelling coding', event_coding_lable)
    assert False


m_remapping =  {'MX': 1 , 'M1A': 2 ,  'M1B' : 2,  "'--" : 1, 'M1' :2 , 'M0': 0}
n_remapping =    {'N1': 2, 'N0': 0 , 'N1B': 2,  'N1A' : 2, 'N1C' : 2, 'N2': 3, "'--": 1 , 'N2A': 3, 'N2B': 3}
t_remapping =  {"'--": 0, 'T4':3, 'T4A':3, 'T2':1, 'T3':2}
#stage_remapping = {'Stage IIIB':4 , 'Stage IV': 7, "'--": 0 , 'Stage IVA': 7, 'Stage IIB': 3, 'Stage IIA': 1, 'Stage III': 5, 'Stage II': 2, 'Stage IIIC' :6, 'Stage I': 0}
stage_remapping = {'STAGE IIIB':2 , 'STAGE IV': 3, "'--": 0 , 'STAGE IVA': 3, 'STAGE IIB': 1, 'STAGE IIA': 1, 'STAGE III': 2, 'STAGE II': 1, 'STAGE IIIC' :2, 'STAGE I': 0}


scalar_drug_rank_ratio_map = {key : float(map_irinotecan_geo_results[key][1])/float(map_irinotecan_geo_results[key][2]) for key in map_irinotecan_geo_results.keys()}
scalar_drug_crc_rank_ratio_map = {key : float(map_irinotecan_geo_results[key][3])/float(map_irinotecan_geo_results[key][4]) for key in map_irinotecan_geo_results.keys()}
scrore_drug_map = {key : map_irinotecan_geo_results[key][0]  for key in map_irinotecan_geo_results.keys()}



def gen_all_clinical_data_maps(reduced_response_map, prediction_class_irinotecan_dict_geo) :

    verbose = True
    CLINICAL_DATA_PATH_NAME = '../../Testing_lifelines'
    CLINICAL_DATA_FILE_NAME = 'irinotecan_coadread_tcga_pan_can_atlas_2018_clinical_data.csv'   

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, '\t', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index('Patient ID')
    m_idx = header.index('American Joint Committee on Cancer Metastasis Stage Code')
    n_idx = header.index('Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code')
    t_idx = header.index('American Joint Committee on Cancer Tumor Stage Code')    
    stage_idx = header.index( 'Neoplasm Disease Stage American Joint Committee on Cancer Code')
    
    age_idx = header.index('Diagnosis Age')
    gender_idx = header.index('Sex')
    diagnosis_idx = header.index('Tumor Type')
    organ_idx = header.index('Tumor Disease Anatomic Site')



    pat_set = set()
    m_map = {}
    m_set = set()
    n_map = {}
    n_set = set()
    t_map = {}
    t_set = set()
    stage_map = {}
    stage_set = set()

    response_set = set()
    prediction_set = set()
  
    for x in acc :
        pat_id = x[patient_idx]
        if pat_id == '' or pat_id in  pat_set or pat_id not in reduced_response_map.keys():
            continue
        
        pat_set.add(pat_id)
        m_map[pat_id] = m_remapping[x[m_idx]]
        m_set.add(m_remapping[x[m_idx]])
        n_map[pat_id] = n_remapping[x[n_idx]]
        n_set.add(n_remapping[x[n_idx]])
        t_map[pat_id] = t_remapping[x[t_idx]]
        t_set.add(t_remapping[x[t_idx]])
        stage_map[pat_id] = stage_remapping[x[stage_idx]]
        stage_set.add(stage_remapping[x[stage_idx]])

        response_set.add(reduced_response_map[pat_id][0])
        prediction_set.add(prediction_class_irinotecan_dict_geo[pat_id])
        
    if verbose :        
        print('m_set : ', len(m_set), m_set)
        printmap(m_map, 'm_map')
        print('n_set : ', len(n_set), n_set)
        printmap(n_map, 'n_map')
        print('t_set : ', len(t_set), t_set)
        printmap(t_map, 't_map')
        print('stage_set : ', len(stage_set), stage_set)
        printmap(stage_map, 'stage_map')
    print('response_set', len(response_set), response_set)
    print(prediction_set, 'prediction_set')

    all_data = []
    for key in  pat_set :
        all_data.append({'Patient_ID' : key ,
                                     'American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced' : m_map[key],
                                     'Neoplasm_Disease_Lymph_Node_Stage_American_Joint_Committee_on_Cancer_Code_Encoded' : n_map[key],
                                      'Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code_Encoded' : t_map[key],
                                     'American_Joint_Committee_on_Cancer_Tumor_Stage_Code_Encoded' :  stage_map[key],
                                      'P_DrugMerge_Prediction' :  prediction_class_irinotecan_dict_geo[key]  ,
                                      'Scalar_prediction' : scalar_drug_rank_ratio_map[key],
                                       'Scalar_crc_prediction' : scalar_drug_crc_rank_ratio_map[key],
                                      'Scalar_score' : scrore_drug_map[key], 
                                       'Response_coded' : response_remapping[reduced_response_map[key][-1][0]],
                                       'Response_duration' : reduced_response_map[key][-1][1]
                                      })
        

    return (m_map, n_map, t_map, stage_map , all_data)
    
(m_map, n_map, t_map, stage_map , all_data_irinotecan) = gen_all_clinical_data_maps(map_irinotecan_response, prediction_class_irinotecan_dict_geo)

printlist(all_data_irinotecan, 'all data_irinotecan')



###################
# start of cox analysis part..
# copied for the fluoroutacil code..with minor modifications...



##############################
from lifelines.utils import k_fold_cross_validation
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import re
import matplotlib.pyplot as plt

        
def testing_parametric_cox_fitters(df, str_formula) :

   

    duration_col = 'Response_duration'
    event_col =  'Response_coded'    
    gen_prediction = False
    k_fold_computation = True

    print("\n" + "="*50 + "\n")
#    knot_values = [4,]
#    knot_values = [2,  3 ]
#    knot_values = [2,  3 , 4, 5]
#    knot_values =  [2,  3,4,5,6,7, 8, 9, 10 ]
    knot_values = [2,3,4,5,6,7,8]

    results = []
    print("\n" + "="*50 + "\n")
    for kv in  knot_values :
        try :
            cph = CoxPHFitter(baseline_estimation_method="spline", n_baseline_knots=kv)
            cph.fit(df, duration_col=duration_col, event_col=event_col, formula=str_formula)
            cph.print_summary()
            lrt_result = cph.log_likelihood_ratio_test()
            p_value = lrt_result.p_value
            log2_p_value = -np.log2(p_value)
            print('recomputed pval and  - log 2 pval ', p_value,  log2_p_value)
            results.append([kv, cph.AIC_, p_value, log2_p_value]) # no concordance index?

            if k_fold_computation  :
                target = "concordance_index"
                for k_fold in knot_values :
#                    print ('target, k = ', target, k_fold)
                    mean_ci = np.mean(k_fold_cross_validation(cph,   df, duration_col=duration_col, event_col=event_col, scoring_method=target, k=k_fold, fitter_kwargs={'formula' : str_formula}))
                    print('target {}, k_fold  {} , mean CI = {}  '.format(target, k_fold, mean_ci))
            
            if gen_prediction :
                print('Predictions on trial data')
                df2 = cph.predict_expectation(df)
                print(df2.head()) # Print the first 5 rows to verify
                print(df2.info()) # Print a summary of the DataFrame
                true_times= df[duration_col].tolist()
                events = df[event_col].tolist()
                predicted_times = df2.tolist()
                counter = 0
                for x in range(len(events)) :
                    if events[x] ==  1 :
                        print('True time {}, predicted {}'.format(true_times[x], predicted_times[x]))
                        counter +=1
                print('Number of predictions . ', counter)

            
            print("\n" + "="*50 + "\n")
        except Exception as e:
            print('Exception with ', kv, e)
            pass

    print('Summary of log2_p_value  for ', str_formula)
    results.sort(key =lambda x : x[3], reverse=True)
    for x in range(len(results)) :
        print(x, results[x])


###############################
        
def gen_cox_analysis_irinotecan_paz_2(df,  formula_str):

#    df = pd.DataFrame(all_data_fluorouracil)

    print('*******************************************************')
    print("\n--- Columns ---")
    print(df.columns.tolist())
    print("\nDataFrame head and info  columns:")
    print(df.head())
    print(df.info())

    duration_col = 'Response_duration'
    event_col =  'Response_coded'    

    cph_formula = CoxPHFitter()

    cph_formula.fit(df,
                    duration_col=duration_col,
                    event_col=event_col,
                    formula=formula_str)

    print("--- Cox Model Summary ---")
    print('formula_str: ', formula_str)
    cph_formula.print_summary()
    print("\n--------------------------\n")

    # Check proportional hazards assumption

    check_PH_flag = True
#    check_PH_flag = False
    if  check_PH_flag :
        import matplotlib.pyplot as plt 
        cph_formula.check_assumptions(df, p_value_threshold=0.05, show_plots=True)
#        plt.show() # Make sure to show the plot if show_plots=True

#    gen_coeff_plot_flag = True
    gen_coeff_plot_flag = False
    if gen_coeff_plot_flag :
        import matplotlib.pyplot as plt 
        cph_formula.plot()    
        plt.show()


###############################
        

if __name__ == "__main__":
    #    pass
    #    gen_cox_analysis_fluorouracil_paz()

    print("**********************************************************************************")
    print("Starting main section")

    f0 = 'C(American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced)'
    f1 = 'American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced'

    f2 = 'C(Neoplasm_Disease_Lymph_Node_Stage_American_Joint_Committee_on_Cancer_Code_Encoded)'
    f2bis = 'Neoplasm_Disease_Lymph_Node_Stage_American_Joint_Committee_on_Cancer_Code_Encoded'

    f3 = 'C(Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code_Encoded)'
    f3bis = 'Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code_Encoded'

    f4 = 'C(American_Joint_Committee_on_Cancer_Tumor_Stage_Code_Encoded)'
    f4bis =  'American_Joint_Committee_on_Cancer_Tumor_Stage_Code_Encoded'

    f5 = 'C(P_DrugMerge_Prediction)'                    
    f6 = 'Scalar_prediction'
    f7 = 'Scalar_crc_prediction'
    f8 =  'Scalar_score'

    formulas = [
          f0, # only methstasis as predictive covariat
         ' + '.join([f0, f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:
         ' + '.join([f1, f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:
           ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:   senza f0 
         ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:   senza f1
           ' + '.join([f2bis,  f3bis , f4bis]), # fine formuae senza indici pdrugmerge 
        ' + '.join([f0,    f5, f6, f7, f8]), # best performant so far               
         ' + '.join([f0,    f6, f7, f8]),  # best senza f5
        ' + '.join([f0,    f5,  f7, f8]), # best senza f6
         ' + '.join([f0,    f5,  f6, f8]), # best senza f7
         ' + '.join([f0,    f5, f6, f7]), # best senza f8
        ' + '.join([  f5, f6, f7, f8]), # best performant so far     senza f0           
         ' + '.join([ f6, f7, f8]),  # best senza f5   senza f0 
        ' + '.join([  f5,  f7, f8]), # best senza f6   senza f0 
         ' + '.join([   f5,  f6, f8]), # best senza f7   senza f0 
         ' + '.join([  f5, f6, f7]), # best senza f8   senza f0 
    
                 ]

    
# variazioni su di una buona configurazioe
    formulas2 = [
     ' + '.join([   f5,  f6, f8]),
      ' + '.join([   f6, f8]),
      ' + '.join([   f5,  f8]),
      ' + '.join([   f5,  f6]),
    f5,
    f6,
    f8
        ]

# variazioni su di una buona configurazioe 2
    formulas3 = [
     ' + '.join([ f0,  f5,  f7, f8]),
      ' + '.join([f0,    f7, f8]),
      ' + '.join([f0,   f5,  f8]),
      ' + '.join([ f0,   f5,  f6]),
    ' + '.join([ f0,f5]),
    ' + '.join([ f0,f7]),
    ' + '.join([ f0,f8]),
        ]


# varaizioni sulla miglioere configurazione x paclitaxel

    formulas4 =[
 ' + '.join([f2bis,  f3bis , f4bis, f6, f7, f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis,  f7, f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f6,  f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f6, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f6]),
 ' + '.join([f2bis,  f3bis , f4bis, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f5]),
 ' + '.join([f2bis,  f3bis , f4bis, f5, f6, f7, f8]),
 ' + '.join([f2bis,  f3bis , f4bis,  f5, f7, f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f5, f6,  f8]),
 ' + '.join([f2bis,  f3bis , f4bis, f5, f6, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f5, f6]),
 ' + '.join([f2bis,  f3bis , f4bis, f5, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f5, f8]),
 ' + '.join([f0,  f2bis , f3bis , f4bis ,f5,  f7, f8 ]) ,
    ]


#    formulas5 = [ ' + '.join([ f2bis , f3bis , f4bis , f7, f8, f5 , f0,  ])  ]

    formulas6 = [f5]

    ###################
 #   formuals without f5

 
    formulas_bis = [
              
 ' + '.join([f0,    f6, f7, f8]),  # best senza f5
' + '.join([f0,     f7, f8]), # best senza f6
 ' + '.join([f0,     f6, f8]), # best senza f7
 ' + '.join([f0,    f6, f7]), # best senza f8
  f0, # only methstasis as predictive covariat
 ' + '.join([f0, f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:
 ' + '.join([f1, f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:      
 ' + '.join([ f6, f7, f8]),  # best senza f5   senza f0 
' + '.join([  f7, f8]), # best senza f6   senza f0 
 ' + '.join([   f6, f8]), # best senza f7   senza f0 
 ' + '.join([ f6, f7]), # best senza f8   senza f0 
 ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:   senza f0 
                  ]

# variazioni su di una buona configurazioe
    formulas2_bis = [
     ' + '.join([   f6, f8]),
    f6,
    f8
        ]

# variazioni su di una buona configurazioe 2
    formulas3_bis = [
      ' + '.join([f0,    f7, f8]),
      ' + '.join([f0,    f8]),
      ' + '.join([ f0,    f6]),
    ' + '.join([ f0,f7,]),
    ' + '.join([ f0,f8]),
        ]


# varaizioni sulla miglioere configurazione x paclitaxel

    formulas4_bis =[
 ' + '.join([f2bis,  f3bis , f4bis, f6, f7, f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis,  f7, f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f6,  f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f6, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f6]),
 ' + '.join([f2bis,  f3bis , f4bis, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f8 ]),

 ' + '.join([f2bis,  f3bis , f4bis, f6, f7, f8]),
 ' + '.join([f2bis,  f3bis , f4bis,  f7, f8]),
 ' + '.join([f2bis,  f3bis , f4bis, f6,  f8 ]),
 ' + '.join([f2bis,  f3bis , f4bis, f6, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f6]),
 ' + '.join([f2bis,  f3bis , f4bis, f7]),
 ' + '.join([f2bis,  f3bis , f4bis, f8 ]),
 ' + '.join([f0,  f2bis , f3bis , f4bis , f7, f8 ]) ,
    ]

# filds to be not considered: metastatic state
#    forbidden = [f0, f1]
    forbidden = []

    df = pd.DataFrame(all_data_irinotecan)
    nan_df = df[df.isnull().any(axis=1)]
    print(nan_df)
    # Create a new DataFrame without any rows containing NaNs
    df_cleaned = df.dropna()
    nan_df_2 = df_cleaned[df_cleaned.isnull().any(axis=1)]
    print(nan_df_2)

    duplicate_formulas = set()
    formulae =   [formulas, formulas2, formulas3, formulas4, formulas6, formulas_bis, formulas2_bis, formulas3_bis, formulas4_bis]
    for formulas in formulae :
        for formula in formulas :
            if formula in  duplicate_formulas :
                print('***Formula already processed, skipped : ', formula)
            elif formula in forbidden  :
                print('***Formula forbidden, skipped : ', formula)
            else :
                for x in forbidden :
                    if x in formula :
                        formula1  =formula.replace(x, "")
                        formula1 = formula1.strip(" +")
                        formula = formula1
                        print('***Formula with forbidden field, changed : ', formula)
                        break
                gen_cox_analysis_irinotecan_paz_2(df_cleaned, formula)
                testing_parametric_cox_fitters(df_cleaned, formula)
                duplicate_formulas.add(formula)           
 
    print ('***Analysed {} formulae'.format(len(duplicate_formulas)))  
