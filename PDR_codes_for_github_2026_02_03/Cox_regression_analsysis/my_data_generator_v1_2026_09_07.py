#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pickle
import argparse


import numpy as np
import pandas as pd


####################
def print_list(inlist, explanationstring, limit = None) :

    print ("*************************************************************")
    print (explanationstring, "size", len(inlist), 'with limit at', limit)
    counter = 0
    for x in range(len(inlist)) :
        print (x, inlist[x])
        counter += 1
        if limit != None :
            if counter  > limit :
                break


#####################

def putdatatodump(dumpdirectory, dumpfilename, savedobject, overwrite):

    filenamefordump = dumpdirectory + "/" + dumpfilename

    if os.path.isfile(filenamefordump):
        print ("Existing  file," ,  filenamefordump)
        if overwrite :
            print("Overwriting" ,  filenamefordump)
            dump = True
        else :
            print("Not Overwriting" ,  filenamefordump)
            dump = False
    else :
        print ("New targets file" ,  filenamefordump)
        dump = True

    if dump :
        filefordump = open(filenamefordump, "wb")
        pickle.dump(savedobject, filefordump)
        filefordump.close()



##########################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Run synthetic nested Cox tests, persist CoxEvaluationResult objects "
            "with pickle, and export LaTeX/PNG graphics."
        )
    )

# data loading arguments


    parser.add_argument('experiment_labled',
                        action='store',
                        type = str,
                        choices=['STAD_fluorouracil_only',
                                 'STAD_fluorouracil_only_subtype',
                                 'CRC_irinotecan',
                                 'CRC_fluorouracil',
                                 'CRC_oxaliplatin',
                                 'CRC_leucovorin',
                                 'BRCA_A_paclitaxel',
                                 'BRCA_A_paclitaxel_subtype',
                                 'BRCA_B_cyclophosphamide',
                                 'BRCA_B_cyclophosphamide_subtype',
                                 'BRCA_C_fluorouracil',
                                 'BRCA_C_fluorouracil_subtype',
                                 'BRCA_D_docetaxel',
                                 'BRCA_D_docetaxel_subtype',
                                 'BRCA_E_doxorubicin',
                                 'BRCA_E_doxorubicin_subtype',
                                 ],
                        help='Defines the  experiment type: choose amnong 16 main experiments')

    parser.add_argument('local_event_coding_lable',
                        action='store',
                        type = str,
                        choices=[ 'PD1_all0', 'PDSD1_all0', 'CRPR1_all0', 'CR1_all0'],
                        help='Defines the  experiment coding for the mian events')


# nested cox aguments
    parser.add_argument(
        "--mode",
        choices=("pickle", "just_load"),
        default="just_load",
        help=("just load: loads but produces no effect. picckle: produces the pickeld output file"
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Exlplicitly overwrite an existign file. defauls is not overwrite."
        ),
    )

    parser.add_argument(
        "--output",
        default=None,
        help="Output filename.",
    )


    args = parser.parse_args()

    print('Experiment {}'.format(args.experiment_labled))
    print('Event coding {}'.format(args.local_event_coding_lable))
    print('mode', args.mode)
    print('overwrite', args.overwrite)

    if args.output == None :
        filename = '_'.join([args.experiment_labled, args.local_event_coding_lable]) + ".pkl"
        args.output = filename
        
    print('output', args.output)

    return args

##################

f0 = 'C(American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced)'
f1 = 'American_Joint_Committee_on_Cancer_Metastasis_Stage_Code_Reduced'


def remove_forbidden(in_list, forbidden) :
    out_list = []
    for formula in in_list :
        found = False
        for x in forbidden :
            if x in formula :
                found = True
                formula1  =formula.replace(x, "")
                formula1 = formula1.strip(" +")
                if formula1 != "" :
                    out_list.append(formula1)
                    print('substitution', formula, '====' ,  formula1)
                break
        if not found :        
            out_list.append(formula)
    print_list(out_list, 'Reduced List of fingerprints')
    return out_list

def df_fingerprint(df: pd.DataFrame) -> int:
    # Stable hash across runs for the *content* (including index)
    return int(pd.util.hash_pandas_object(df, index=True).sum())





def loading_main_data(args) :

    experiment_lable = args.experiment_labled

    #STAD
    if experiment_lable == 'STAD_fluorouracil_only' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import  df, formulas_clinical_only, formulas_clinical_plus_PDR_indices, event_coding_lable
        
    elif experiment_lable ==  'STAD_fluorouracil_only_subtype' :
        from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import  df, formulas_clinical_only, formulas_clinical_plus_PDR_indices, event_coding_lable
        
    # CRC
    elif experiment_lable ==  'CRC_irinotecan' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_CRC_analysis_irinotecan_patients  import  df_cleaned, event_coding_lable  
        df = df_cleaned
        
    elif  experiment_lable ==  'CRC_fluorouracil' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_CRC_analysis_fluorouracil_patients import df_cleaned, event_coding_lable  # file da duplicare ed editare
        from clinical_data_analysis_CRC import sub_population_flag
        assert  sub_population_flag == 'fl_only',  "Error: Make sure the subpopluation flag  is aligned with the local  experiment lable CRC_fluorouraci and fl_only"
        df = df_cleaned
        formulas_clinical_only = remove_forbidden(formulas_clinical_only, [f0, f1])
        formulas_clinical_plus_PDR_indices = remove_forbidden(formulas_clinical_plus_PDR_indices, [f0, f1])
        
    elif experiment_lable == 'CRC_oxaliplatin' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_CRC_analysis_fluorouracil_patients import df_cleaned,  event_coding_lable  # file da duplicare ed editare
        from clinical_data_analysis_CRC import sub_population_flag
        assert  sub_population_flag == 'fl_ox',  "Error: Make sure the subpopluation flag  is aligned with the local  experiment lable CRC_oxaliplatin and fl_ox "
        df = df_cleaned
        formulas_clinical_only = remove_forbidden(formulas_clinical_only, [f0, f1])
        formulas_clinical_plus_PDR_indices = remove_forbidden(formulas_clinical_plus_PDR_indices, [f0, f1])   

    elif experiment_lable == 'CRC_leucovorin' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_CRC_analysis_fluorouracil_patients import df_cleaned, event_coding_lable  # file da duplicare ed editare
        from clinical_data_analysis_CRC import sub_population_flag
        assert  sub_population_flag == 'fl_leu',  "Error: Make sure the subpopluation flag  is aligned with the local  experiment lable CRC_leucovrin and fl_leu "
        df = df_cleaned
        formulas_clinical_only = remove_forbidden(formulas_clinical_only, [f0, f1])
        formulas_clinical_plus_PDR_indices = remove_forbidden(formulas_clinical_plus_PDR_indices, [f0, f1])   

    # BRCA
    elif experiment_lable ==  'BRCA_A_paclitaxel' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from   Cox_BRCA_paclitaxel_analysis import df_cleaned,  event_coding_lable  
        df = df_cleaned
        
    elif experiment_lable == 'BRCA_A_paclitaxel_subtype' :
        from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from   Cox_BRCA_paclitaxel_analysis_with_subtypes import df_cleaned,  event_coding_lable 
        df = df_cleaned
        
    elif experiment_lable ==  'BRCA_B_cyclophosphamide' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_docetaxel_analysis import   df_cleaned,  event_coding_lable 
        df = df_cleaned
         
    elif experiment_lable ==  'BRCA_B_cyclophosphamide_subtype':
        from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_docetaxel_analysis_with_subtypes import  df_cleaned,  event_coding_lable 
        df = df_cleaned
         
    elif experiment_lable == 'BRCA_C_fluorouracil':
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_fluorouracil_analysis  import   df_cleaned,  event_coding_lable   
        df = df_cleaned
         
    elif experiment_lable ==  'BRCA_C_fluorouracil_subtype' :
        from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_fluorouracil_analysis_with_subtypes import   df_cleaned,  event_coding_lable
        df = df_cleaned

    elif experiment_lable == 'BRCA_D_docetaxel' :
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_docetaxel_D_2_analysis import  df_cleaned,  event_coding_lable 
        df = df_cleaned
         
    elif experiment_lable == 'BRCA_D_docetaxel_subtype' :
        from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_docetaxel_D_2_analysis_with_subtypes import  df_cleaned,  event_coding_lable 
        df = df_cleaned
         
    elif experiment_lable == 'BRCA_E_doxorubicin':
        from Cox_STAD_fluorouracil_only_analysis_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_doxorubicin_analysis  import  df_cleaned,  event_coding_lable 
        df = df_cleaned
         
    elif experiment_lable == 'BRCA_E_doxorubicin_subtype' :
        from Cox_STAD_fluorouracil_only_analysis_subtype_iterated import formulas_clinical_only, formulas_clinical_plus_PDR_indices
        from Cox_BRCA_doxorubicin_analysis_with_subtypes  import df_cleaned,  event_coding_lable 
        df = df_cleaned

    else :
        print("Wrong experiment_lable : ", experiment_lable )
        assert False

    assert args.local_event_coding_lable == event_coding_lable , "Error: Make sure the imported event_coding_lable is aligned with the local value"

    df = df.rename(columns={
    "Neoplasm_Disease_Lymph_Node_Stage_American_Joint_Committee_on_Cancer_Code_Encoded": "Neoplasm_Disease_Lymph_Node_Stage_AJCCC_Encoded"})


    print("fingerprint:", df_fingerprint(df))
    print(df.dtypes)
    df = df.sort_values('Patient_ID').reset_index(drop=True)
    print("fingerprint:", df_fingerprint(df))
    print(df.dtypes)

    return df , formulas_clinical_only, formulas_clinical_plus_PDR_indices



def main():
    args = parse_arguments()

    base_df , clinical_only_FORMULAS, pdr_and_clinical_FORMULAS = loading_main_data(args)

    print_list(clinical_only_FORMULAS, "clinical_only_FORMULAS")
    print_list(pdr_and_clinical_FORMULAS, "pdr_and_clinical_FORMULAS")

# pickle the pandas frame

    if args.mode == "pickle" :
        
        putdatatodump('.', args.output , base_df , args.overwrite)
        

    



if __name__ == "__main__":
    main()
