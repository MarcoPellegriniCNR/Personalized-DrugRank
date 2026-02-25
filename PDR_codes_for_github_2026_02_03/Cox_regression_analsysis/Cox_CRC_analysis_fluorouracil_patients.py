# cox analysis on fluorouracil 43 patients data set
#  This file has  the time data and uses one-hot endoding instead of nuemrical ranking-encoding.
# versioen modificata con ni piu' signficativi dei parametri

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


# predictions al 10\% of the drugs.

#from clinical_data_analysis2 import list_CRC_fluorouracil_good_geo_prediction_ids, list_CRC_fluorouracil_bad_geo_prediction_ids
#from clinical_data_analysis2 import map_CRC_fluorouracil_geo_results,  gen_drugs_to_check_maps

from clinical_data_analysis_CRC import list_CRC_fluorouracil_good_geo_prediction_ids, list_CRC_fluorouracil_bad_geo_prediction_ids
from clinical_data_analysis_CRC  import map_CRC_fluorouracil_geo_results,  gen_drugs_to_check_maps


def initial_data_generation() :
    
# vecchia formulazione
##    prediction_class_fluorouracil_dict_geo= {}
##    for key in list_fluorouracil_good_geo_prediction_ids :
##        prediction_class_fluorouracil_dict_geo[key] = 1
##    for key in list_fluorouracil_bad_geo_prediction_ids :
##        prediction_class_fluorouracil_dict_geo[key] = 0

# nuova formulazione
    prediction_class_fluorouracil_dict_geo= {}
    for key in list_CRC_fluorouracil_good_geo_prediction_ids :
        prediction_class_fluorouracil_dict_geo[key] = 0
    for key in list_CRC_fluorouracil_bad_geo_prediction_ids :
        prediction_class_fluorouracil_dict_geo[key] = 1


    printmap(prediction_class_fluorouracil_dict_geo, 'prediction_class_fluorouracil_dict_geo')

    # add othe nuemrical functions obtained form the map_fluorouracil_geo_results

    all_ids =  list_CRC_fluorouracil_good_geo_prediction_ids + list_CRC_fluorouracil_bad_geo_prediction_ids

    return prediction_class_fluorouracil_dict_geo

# produce la mappa delle rispsote gia' ridotta.
prediction_class_fluorouracil_dict_geo  = initial_data_generation()
(treatments_map, drugs_to_check_lable, reduced_response_map) = gen_drugs_to_check_maps()

scalar_drug_rank_ratio_map = {key : float(map_CRC_fluorouracil_geo_results[key][1])/float(map_CRC_fluorouracil_geo_results[key][2]) for key in map_CRC_fluorouracil_geo_results.keys()}
scalar_drug_crc_rank_ratio_map = {key : float(map_CRC_fluorouracil_geo_results[key][3])/float(map_CRC_fluorouracil_geo_results[key][4]) for key in map_CRC_fluorouracil_geo_results.keys()}
score_drug_map = {key : map_CRC_fluorouracil_geo_results[key][0]  for key in map_CRC_fluorouracil_geo_results.keys()}


#event_coding_lable = 'PD1_all0'
#event_coding_lable = 'CRPR1_all0'
#event_coding_lable = 'CR1_all0'
#event_coding_lable = 'PDSD1_all0'
#event_coding_lable = 'PD2_CRPR1_all0'

try:
    event_coding_lable
except NameError:
    print("event_coding_lable is not defined: please edit the script uncommenting one of the possible values!")
else:
    print("event_coding_lable is  defined.")

print('Event coding : ', event_coding_lable)

if event_coding_lable == 'PD1_all0' :
    # progressive disease = 1 all other 0
    print('Coding Progressive disease = 1, all other=0')
    response_remapping =  {'Partial Response': 0, 'Complete Response': 0, 'Stable Disease':0, 'Progressive Disease':1}

elif event_coding_lable == 'PDSD1_all0' :
    # progressive disease = 1 all other 0
    print('Coding Progressive disease + stable disease = 1, all other=0')
    response_remapping =  {'Partial Response': 0, 'Complete Response': 0, 'Stable Disease':1, 'Progressive Disease':1}


elif event_coding_lable == 'CR1_all0' :
        #compete response 1 all otehr 0
    print('Coding Progressive + stable disease +  partial  response= 0, complete response=1')
    response_remapping =  {'Partial Response': 0, 'Complete Response': 1, 'Stable Disease':0, 'Progressive Disease':0}

elif event_coding_lable == 'CRPR1_all0' :
    #compete response + partial reponse 1 all otehr 0
    print('Coding Progressive + stable disease = 0, partial + complete response=1')
    response_remapping =  {'Partial Response': 1, 'Complete Response': 1, 'Stable Disease':0, 'Progressive Disease':0}

elif event_coding_lable == 'CRPRSD1_all0' :
    #compete response + partial reponse + stable disease  1 all otehr 0
    print('Coding Progressive disease = 0, all other=1')
    response_remapping =  {'Partial Response': 1, 'Complete Response': 1, 'Stable Disease':1, 'Progressive Disease':0}

elif event_coding_lable == 'PD2_CRPR1_all0' :
    #compete response + partial reponse 1 , progressive disease 2 , all otehr 0
    print('Coding  stable disease = 0, progeressive disease 2, complete response+partial repsonse =1')
    response_remapping =  {'Partial Response': 1, 'Complete Response': 1, 'Stable Disease':0, 'Progressive Disease':2}
else :
    print('Error in event labelling coding', event_coding_lable)
    assert False



# vecchia remapping x metastasis :
m_remapping =  {'MX': 1 , 'M1a': 2 , "'--" : 1, 'M1' :2 , 'M0': 0}
# nuova remapping x metastasis :
#m_remapping =  {'MX': 1 , 'M1a': 1 , "'--" : 1, 'M1' :1 , 'M0': 0}

n_remapping =    {'N1': 2, 'N0': 0 , 'N1b': 2, 'N2': 3, "'--": 1 , 'N2a': 3, 'N2b': 3}
t_remapping =  {"'--": 0, 'T4':3, 'T4a': 3, 'T2':1, 'T3':2}
#stage_remapping = {'Stage IIIB':4 , 'Stage IV': 7, "'--": 0 , 'Stage IVA': 7, 'Stage IIB': 3, 'Stage IIA': 1, 'Stage III': 5, 'Stage II': 2, 'Stage IIIC' :6, 'Stage I': 0}
stage_remapping = {'Stage IIIB':2 , 'Stage IV': 3, "'--": 0 , 'Stage IVA': 3, 'Stage IIB': 1, 'Stage IIA': 1, 'Stage III': 2, 'Stage II': 1, 'Stage IIIC' :2, 'Stage I': 0}

#################


def gen_all_clinical_data_maps(reduced_response_map, prediction_class_fluorouracil_dict_geo) :

    verbose = False
    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'
    CLINICAL_DATA_FILE_NAME = 'dati_clinici_pz_fluorouracil_CRC.csv'   

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, ';', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index('\ufeffcase_id')
    m_idx = header.index('clinical (2).diagnoses.ajcc_pathologic_m')
    n_idx = header.index('clinical (2).diagnoses.ajcc_pathologic_n')
    t_idx = header.index('clinical (2).diagnoses.ajcc_pathologic_t')    
    stage_idx = header.index( 'clinical (2).diagnoses.ajcc_pathologic_stage')
    
    age_idx = header.index('clinical (2).demographic.age_at_index')
    gender_idx = header.index('clinical (2).demographic.gender')
    diagnosis_idx = header.index('clinical (2).diagnoses.primary_diagnosis')
    organ_idx = header.index('clinical (2).diagnoses.tissue_or_organ_of_origin')



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
        if (pat_id == '') or (pat_id in  pat_set) or (pat_id not in  map_CRC_fluorouracil_geo_results.keys()):
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
        prediction_set.add(prediction_class_fluorouracil_dict_geo[pat_id])
        
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
#    nomi vecchi                     
##                                      'P_DrugMerge_Prediction' :  prediction_class_fluorouracil_dict_geo[key]  ,
##                                       'Scalar_prediction' : scalar_drug_rank_ratio_map[key],
##                                       'Scalar_crc_prediction' : scalar_drug_crc_rank_ratio_map[key],
##                                       'Scalar_score' : score_drug_map[key],
# nomi nuovi
                                      'P_DrugMerge_Prediction' :  prediction_class_fluorouracil_dict_geo[key]  ,
                                       'Scalar_prediction' : scalar_drug_rank_ratio_map[key],
                                       'Scalar_crc_prediction' : scalar_drug_crc_rank_ratio_map[key],
                                       'Scalar_score' : score_drug_map[key],                         
#                         
                                       'Response_coded' : response_remapping[reduced_response_map[key][-1][0]],
                                       'Response_duration' : reduced_response_map[key][-1][1]
                                      })

        
        

    return (m_map, n_map, t_map, stage_map , all_data)
    
(m_map, n_map, t_map, stage_map , all_data_fluorouracil) = gen_all_clinical_data_maps(reduced_response_map, prediction_class_fluorouracil_dict_geo)

printlist(all_data_fluorouracil, 'all data')

#######################

def dump_to_file_for_R(directory, filename, all_data_list_of_dictionaries) :

    assert len(all_data_list_of_dictionaries) > 0
    
    fields_list =  all_data_list_of_dictionaries[0].keys()

    f = open(os.path.join(directory,filename), 'w', encoding="utf-8")

    f.write(' , '.join(fields_list)+ '\n')
    
    for dictionary in all_data_list_of_dictionaries :
        data = [str(dictionary[x]) for x in fields_list ]
        f.write(' , '.join(data)+ '\n')

    f.close()

    



#dump_to_file_for_R('.',  event_coding_lable + '_test_dumping_CRC_fluorouracil_for_R.csv',  all_data_fluorouracil)
      


#######################




import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import re
import matplotlib.pyplot as plt

def gen_cox_analysis_fluorouracil_paz() :

    df = pd.DataFrame(all_data_fluorouracil)

    print('*******************************************************')
    print("\n--- Columns ---")
    print(df.columns.tolist())
    print("\nDataFrame head and info  columns:")
    print(df.head())
    print(df.info())

    duration_col = 'Response_duration'
    event_col =  'Response_coded'

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

    #formula_str = ' + '.join([f1, f2,  f3, f4,  f5])
    #formula_str = ' + '.join([f0, f2,  f3, f4,  f5])
    #formula_str = ' + '.join([f0,  f5])
    #formula_str = f0
    #formula_str = f5
    #formula_str = ' + '.join([f0,  f6])
    #formula_str = ' + '.join([f0,  f6, f7, f8 ])
    #formula_str = ' + '.join([ f6, f7, f8 ])


    #formula_str = ' + '.join([f5,  f6, f7, f8 ])
    # all
#    formula_str = ' + '.join([f0, f2,  f3, f4,  f5, f6, f7, f8])
#    formula_str = ' + '.join([f0, f2,   f4,  f5, f6, f7, f8])
# all senza C
#    formula_str = ' + '.join([f1, f2bis,  f3bis, f4bis,  f5, f6, f7, f8])
#    formula_str = ' + '.join([f0, f2bis,  f3bis, f4bis,   f6, f7, f8])
    formula_str = ' + '.join([f0, f2bis,  f3bis, f4bis ])
    
    #formula_str = ' + '.join([f0, f3, f4,  f5, f6, f7, f8])
    #formula_str = ' + '.join([f0,  f4,  f5, f6, f7, f8])

#################################################
# sistematic analysis

# best so far 
#    formula_str = ' + '.join([f0,    f5, f6, f7, f8]) # best performant so far
# best senza f5
#    formula_str = ' + '.join([f0,    f6, f7, f8])
# best senza f6
#    formula_str = ' + '.join([f0,    f5,  f7, f8])
# best senza f7
#    formula_str = ' + '.join([f0,    f5,  f6, f8])
# best senza f8
#    formula_str = ' + '.join([f0,    f5, f6, f7])

# solo metastasi
#    formula_str = f0 # only methstasis as predictive covariate

# solo mnt e stage com m onehot enceoded:
#    formula_str = ' + '.join([f0, f2bis,  f3bis , f4bis])
# solo mnt e stage com m nuemrico:
#    formula_str = ' + '.join([f1, f2bis,  f3bis , f4bis])


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

    check_PH_flag = False
    if  check_PH_flag :
        import matplotlib.pyplot as plt 
        cph_formula.check_assumptions(df, p_value_threshold=0.05, show_plots=True)
        plt.show() # Make sure to show the plot if show_plots=True

    gen_coeff_plot_flag = True
    if gen_coeff_plot_flag :
        import matplotlib.pyplot as plt 
        cph_formula.plot()    
        plt.show()

        #######################

def gen_cox_analysis_fluorouracil_paz_2(df,  formula_str):

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


########################
from lifelines import AalenAdditiveFitter

def testing_AalenAdditiveFitter_fitters(df, str_formula) :


    print("\n" + "="*50 + "\n")
    coef_penalizer_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#    rossi_dataset['T'] = rossi_dataset[duration_col]
#   rossi_dataset['E'] = rossi_dataset[event_col]
    print(df.head()) # Print the first 5 rows to verify
    print(df.info()) # Print a summary of the DataFrame

    duration_col = 'Response_duration'
    event_col =  'Response_coded'    

    results = []
    print("\n" + "="*50 + "\n")
    for kv in  coef_penalizer_values :
        try :
            cph = AalenAdditiveFitter(coef_penalizer=kv, fit_intercept=False)
            cph.fit(df , duration_col,   event_col=event_col, formula=str_formula)
            cph.print_summary()
#            lrt_result = cph.log_likelihood_ratio_test()
#            p_value = lrt_result.p_value
#            log2_p_value = -np.log2(p_value)
#            print('recomputed pval and  - log 2 pval ', p_value,  log2_p_value)
            results.append([kv, cph.concordance_index_]) # no concordance index?
            print("\n" + "="*50 + "\n")
        except Exception as e:
            print('Exception with ', kv)
            print(e)
            pass

    print('Summary of concordance index for ', str_formula)
    results.sort(key =lambda x : x[1], reverse=True)
    for x in range(len(results)) :
        print(x, results[x])

##########################

##from lifelines import GeneralizedGammaRegressionFitter
##
##def testing_penalization_paramteres_GeneralizedGammaRegressionFitter(df , str_formula) :
##
##    duration_col = 'Response_duration'
##    event_col =  'Response_coded'
##
##
##    print("\n" + "="*50 + "\n")
##    penalizer_values = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08 , 0.09 , 0.1]
##    l1_ratio_values = [0.0, 0.1, 0.2, 0.3, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
##
##    results = []
##    print("\n" + "="*50 + "\n")
##    for pv in penalizer_values :
##        for rv in l1_ratio_values :
##            try :
##                cph = GeneralizedGammaRegressionFitter(penalizer=pv, l1_ratio=rv)
##                cph.fit(df, duration_col=duration_col, event_col=event_col, formula=str_formula)
##                cph.print_summary()
##                lrt_result = cph.log_likelihood_ratio_test()
##                p_value = lrt_result.p_value
##                log2_p_value = -np.log2(p_value)
##                print('recomputed pval and  - log 2 pval ', p_value,  log2_p_value)
##                results.append([pv, rv, cph.concordance_index_, cph.AIC_, p_value, log2_p_value])
##                print("\n" + "="*50 + "\n")
##            except Exception as e:
##                print('Exception with ', pv, rv, e)
##                pass
##            
##    print('Summary of concordance index for ', str_formula)
##    results.sort(key =lambda x : x[5], reverse=True)
##    for x in range(len(results)) :
##        print(x, results[x])
        
##############################
from lifelines.utils import k_fold_cross_validation
        
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


###################################

def old_main() :  #  code with teh old main with only 15 formulae.

    #    pass
    #    gen_cox_analysis_fluorouracil_paz()

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
' + '.join([f0,    f5, f6, f7, f8]), # best performant so far               
 ' + '.join([f0,    f6, f7, f8]),  # best senza f5
' + '.join([f0,    f5,  f7, f8]), # best senza f6
 ' + '.join([f0,    f5,  f6, f8]), # best senza f7
 ' + '.join([f0,    f5, f6, f7]), # best senza f8
  f0, # only methstasis as predictive covariat
 ' + '.join([f0, f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:
 ' + '.join([f1, f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:
' + '.join([  f5, f6, f7, f8]), # best performant so far     senza f0           
 ' + '.join([ f6, f7, f8]),  # best senza f5   senza f0 
' + '.join([  f5,  f7, f8]), # best senza f6   senza f0 
 ' + '.join([   f5,  f6, f8]), # best senza f7   senza f0 
 ' + '.join([  f5, f6, f7]), # best senza f8   senza f0 
 ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:   senza f0 
 ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:   senza f1    
                 ]

    df = pd.DataFrame(all_data_fluorouracil)

    for formula in formulas:
        pass
        gen_cox_analysis_fluorouracil_paz_2(df, formula)
#        testing_AalenAdditiveFitter_fitters(df, formula)
#        testing_penalization_paramteres_GeneralizedGammaRegressionFitter(df, formula)
        testing_parametric_cox_fitters(df, formula)
    

###############################
        

if __name__ == "__main__":
    #    pass
    #    gen_cox_analysis_fluorouracil_paz()

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
    forbidden = [f0, f1]

    df = pd.DataFrame(all_data_fluorouracil)
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
                gen_cox_analysis_fluorouracil_paz_2(df_cleaned, formula)
                testing_parametric_cox_fitters(df_cleaned, formula)
                duplicate_formulas.add(formula)           
 
    print ('***Analysed {} formulae'.format(len(duplicate_formulas)))    
