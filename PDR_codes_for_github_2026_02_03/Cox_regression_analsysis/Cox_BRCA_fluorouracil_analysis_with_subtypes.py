# Cox analysis of BRCA fluorouracil patients - 26 patients



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


from clinical_data_analysis2 import list_brca_fluorouracil_good_geo_prediction_ids_2, list_brca_fluorouracil_bad_geo_prediction_ids_2
from clinical_data_analysis2 import map_BRCA_fluorouracil_geo_results, BRCA_gen_drugs_to_check_maps
from subtypes_testing_2025_12_18 import  res_map_sub_BRCA_fluorouracil_0


def initial_data_generation() :

    # nuova formulazione
    prediction_class_fluorouracil_dict_geo= {}
    for key in list_brca_fluorouracil_good_geo_prediction_ids_2 :
        prediction_class_fluorouracil_dict_geo[key] = 0
    for key in list_brca_fluorouracil_bad_geo_prediction_ids_2 :
        prediction_class_fluorouracil_dict_geo[key] = 1


    printmap(prediction_class_fluorouracil_dict_geo, 'prediction_class_BRCA_fluorouracil_dict_geo')

    return prediction_class_fluorouracil_dict_geo


# produce la mappa delle rispsote gia' ridotta.
prediction_class_fluorouracil_dict_geo  = initial_data_generation()
target_drugs = [ 'Fluorouracil']
(treatments_map, drugs_to_check_lable, reduced_response_map) = BRCA_gen_drugs_to_check_maps(target_drugs)

printmap(treatments_map, 'treatments_map')
printmap(drugs_to_check_lable, 'drugs_to_check_lable')
printmap(reduced_response_map, 'reduced_response_map')


scalar_drug_rank_ratio_map = {key : float(map_BRCA_fluorouracil_geo_results[key][1])/float(map_BRCA_fluorouracil_geo_results[key][2]) for key in map_BRCA_fluorouracil_geo_results.keys()}
scalar_drug_crc_rank_ratio_map = {key : float(map_BRCA_fluorouracil_geo_results[key][3])/float(map_BRCA_fluorouracil_geo_results[key][4]) for key in map_BRCA_fluorouracil_geo_results.keys()}
score_drug_map = {key : map_BRCA_fluorouracil_geo_results[key][0]  for key in map_BRCA_fluorouracil_geo_results.keys()}

printmap(scalar_drug_rank_ratio_map, 'scalar_drug_rank_ratio_map')
printmap(scalar_drug_crc_rank_ratio_map, 'scalar_drug_crc_rank_ratio_map')
printmap(score_drug_map, 'score_drug_map')


#event_coding_lable = 'PD1_all0'
#event_coding_lable = 'PD2_CRPR1_all0'
#event_coding_lable ='CRPR1_all0'
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
    print('Codining Progressive disease = 1, all other=0')
    response_remapping =  {'Partial Response': 0, 'Complete Response': 0, 'Stable Disease':0, 'Progressive Disease':1}

elif event_coding_lable == 'CR1_all0' :
        #compete response 1 all otehr 0
    print('Codining Progressive + stable disease +  partial  response= 0, complete response=1')
    response_remapping =  {'Partial Response': 0, 'Complete Response': 1, 'Stable Disease':0, 'Progressive Disease':0}

elif event_coding_lable == 'CRPR1_all0' :
    #compete response + partial reponse 1 all otehr 0
    print('Codining Progressive + stable disease = 0, partial + complete response=1')
    response_remapping =  {'Partial Response': 1, 'Complete Response': 1, 'Stable Disease':0, 'Progressive Disease':0}

elif event_coding_lable == 'CRPRSD1_all0' :
    #compete response + partial reponse + stable disease  1 all otehr 0
    print('Codining Progressive disease = 0, all other=1')
    response_remapping =  {'Partial Response': 1, 'Complete Response': 1, 'Stable Disease':1, 'Progressive Disease':0}

elif event_coding_lable == 'PD2_CRPR1_all0' :
    #compete response + partial reponse 1 , progressive disease 2 , all otehr 0
    print('Codining  stable disease = 0, progeressive disease 2, complete response+partial repsonse =1')
    response_remapping =  {'Partial Response': 1, 'Complete Response': 1, 'Stable Disease':0, 'Progressive Disease':2}
else :
    print('Error in event labelling coding', event_coding_lable)
    assert False


m_remapping =  {'MX': 1 , 'M1a': 2 , "'--" : 1, 'M1' :2 , 'M0': 0, 'cM0 (i+)': 1}

t_remapping =  {"'--": 0, 'T1':0, 'T1a' : 0, 'T1b': 0, 'T1c': 0,  'T2':1, 'T3':2, 'T3a': 2 ,   'T4':3, 'T4a': 3}

stage_remapping = {"'--": 0, 'Stage I': 0, 'Stage IA': 0, 'Stage II': 1, 'Stage IIB': 1, 'Stage IIA': 1,   'Stage III': 2,  'Stage IIIC' :2, 'Stage IIIB':2 , 'Stage IIIA': 2,  'Stage IV': 3, 'Stage IVA': 3 }

n_remapping =    {'N0': 0 , 'NX':1 ,  'N0 (i+)':1, 'N0 (i-)':1,   "'--": 1 , 'N1': 2,  'N1b': 2, 'N1mi':2, 'N1c':2, 'N1a':2,  'N2': 3,   'N2a': 3, 'N2b': 3, 'N3':4 , 'N3a':4}

sub_type_remapping = {'BRCA_Basal' : 1 , 'BRCA_Her2' :2 ,  'NA' :3  ,  'BRCA_LumB':4, 'BRCA_Normal':5, 'BRCA_LumA' : 6}

#def gen_all_clinical_data_maps(reduced_response_map, prediction_class_fluorouracil_dict_geo) :
def gen_all_clinical_data_maps() :

#    verbose = False
    verbose = True

    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'
    CLINICAL_DATA_FILE_NAME = 'BRCA_all_treatments_2.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, '\t', 'verbose')

    acc = res[0]
    header = res[1]

    patient_idx =  header.index('cases.submitter_id')
    m_idx = header.index('diagnoses.ajcc_pathologic_m')
    n_idx = header.index('diagnoses.ajcc_pathologic_n')
    t_idx = header.index('diagnoses.ajcc_pathologic_t')
    stage_idx = header.index('diagnoses.ajcc_pathologic_stage')


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
        if (pat_id == '') or (pat_id in  pat_set) or (pat_id not in  map_BRCA_fluorouracil_geo_results.keys()) :
            continue
        
        pat_set.add(pat_id)
        m_set.add(x[m_idx])
        n_set.add(x[n_idx])
        t_set.add(x[t_idx])
        stage_set.add(x[stage_idx])

        m_map[pat_id] = m_remapping[x[m_idx]]
        n_map[pat_id] = n_remapping[x[n_idx]]
        t_map[pat_id] = t_remapping[x[t_idx]]
        stage_map[pat_id] = stage_remapping[x[stage_idx]]

        response_set.add(reduced_response_map[pat_id][0])
        prediction_set.add(prediction_class_fluorouracil_dict_geo[pat_id])


    if verbose: 
        print(m_set)
        print(n_set)
        print(t_set)
        print(stage_set)
        printmap(m_map, 'm_map')
        printmap(n_map, 'n_map')
        printmap(t_map, 't_map')
        printmap(stage_map, 'stage_map')
        print('missing patients ', set(map_BRCA_fluorouracil_geo_results.keys()).difference(pat_set))
        print('response_set', len(response_set), response_set)
        print( 'prediction_set', prediction_set,)


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
                                       'Response_duration' : reduced_response_map[key][-1][1],
                                        'Subtype_coded' : sub_type_remapping[res_map_sub_BRCA_fluorouracil_0[key] if key in res_map_sub_BRCA_fluorouracil_0.keys() else 'NA']
                                      })

    return (m_map, n_map, t_map, stage_map , all_data)


####
# run

#gen_all_clinical_data_maps(reduced_response_map, prediction_class_fluorouracil_dict_geo)
(m_map, n_map, t_map, stage_map , all_data_BRCA_fluorouracil) = gen_all_clinical_data_maps()

printlist(all_data_BRCA_fluorouracil, 'all data')

##############################################
# eventual dumping for use in R

def dump_to_file_for_R(directory, filename, all_data_list_of_dictionaries) :

    assert len(all_data_list_of_dictionaries) > 0
    print('Begin dumping in : ', directory, filename)
    
    fields_list =  all_data_list_of_dictionaries[0].keys()
    f = open(os.path.join(directory,filename), 'w', encoding="utf-8")
    f.write(' , '.join(fields_list)+ '\n')
    
    for dictionary in all_data_list_of_dictionaries :
        data = [str(dictionary[x]) for x in fields_list ]
        f.write(' , '.join(data)+ '\n')
    f.close()   
    print('End dumping in : ', directory, filename)
####   

# run
dump_to_file_for_R('.',  event_coding_lable + '_test_dumping_brca_fluorouracil_for_R.csv',  all_data_BRCA_fluorouracil)
           

###############################################



import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import re
import matplotlib.pyplot as plt

#from Cox_2_analysis_fluorouracil_patients_bis import gen_cox_analysis_fluorouracil_paz_2

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


##################

##############################
from lifelines.utils import k_fold_cross_validation
        
def testing_parametric_cox_fitters(df, str_formula) :

   

    duration_col = 'Response_duration'
    event_col =  'Response_coded'    
    gen_prediction = False
    k_fold_computation = True

    print("\n" + "="*50 + "\n")
#    knot_values = [4,]
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
#                for k_fold in [2,  3,4,5,6,7, 8, 9, 10 ] :
                for k_fold in [2,  3 , 4 , 5, 6] :
#                    print ('target, k = ', target, k_fold)
                    print('target, k_fold , mean CI = ', target, k_fold, np.mean(k_fold_cross_validation(cph,   df, duration_col=duration_col, event_col=event_col, scoring_method=target, k=k_fold, fitter_kwargs={'formula' : str_formula})))
            
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


################


if __name__ == "__main__":
    #    pass
    #    gen_cox_analysis_paclitaxel_paz()

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
#    f9 = 'Subtype_coded'  #  ordered numerical version
#    f9 = 'C(Subtype_coded)' # categerical version
    try:
        f9
    except NameError:
        print("f9  is not defined: please edit the script uncommenting one of the possible values!")
    else:
        print("f9 is  defined.")

        
    formulas = [
          f0, # only methstasis as predictive covariat
         ' + '.join([f0, f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:
         ' + '.join([f1, f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:
           ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m onehot enceoded:   senza f0 
         ' + '.join([f2bis,  f3bis , f4bis]), # solo mnt e stage com m nuemrico:   senza f1
 #          ' + '.join([f2bis,  f3bis , f4bis]), # fine formuae senza indici pdrugmerge 
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
    #####################

    # filds to be not considered: metastatic state
#    forbidden = [f0, f1]
    forbidden = []

    df = pd.DataFrame(all_data_BRCA_fluorouracil)
    nan_df = df[df.isnull().any(axis=1)]
    print(nan_df)
    # Create a new DataFrame without any rows containing NaNs
    df_cleaned = df.dropna()
    nan_df_2 = df_cleaned[df_cleaned.isnull().any(axis=1)]
    print(nan_df_2)



    duplicate_formulas = set()
    formula_counter = 0
    formulae =   [formulas, formulas2, formulas3, formulas4, formulas6, formulas_bis, formulas2_bis, formulas3_bis, formulas4_bis]
    for formulas in formulae :
        for formula in formulas :
            formula2 = formula + ' + ' + f9
            formula_counter += 1
            print('************************************************************************')
            print('Processing Formula number {} , values : {} '.format(formula_counter , formula2)) 
            if formula2 in  duplicate_formulas :
                print('***Formula already processed, skipped : ', formula2)
            elif formula2 in forbidden  :
                print('***Formula forbidden, skipped : ', formula2)
            else :
                for x in forbidden :
                    if x in formula2 :
                        formula1  =formula2.replace(x, "")
                        formula1 = formula1.strip(" +")
                        formula2 = formula1
                        print('***Formula with forbidden field, changed : ', formula2)
                        break
                try :
                    gen_cox_analysis_fluorouracil_paz_2(df_cleaned, formula2)
                except Exception as e:
                    print(f"Error occurred: {e}")
                try :
                    testing_parametric_cox_fitters(df_cleaned, formula2)
                except Exception as e:
                    print(f"Error occurred: {e}")
                    
                duplicate_formulas.add(formula2)           
 
    print ('***Analysed {} formulae'.format(len(duplicate_formulas)))  
