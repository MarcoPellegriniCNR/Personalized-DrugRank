# command line version of new_loader_routines_modulation.py

from new_loader_routines_modulation import main_function

import argparse



def main_function_2():
    parser = argparse.ArgumentParser(description='Command line to analyze drug directions in patient data from enrichR.')

#obbligatorio
    parser.add_argument('patient_id',
                        action='store',
                        type = str,
 #                       choices=["brc", "prc", 'aml', 'pfm', 'brc2'],
                        help='Defines the patient id')

#obbligatorio
    parser.add_argument('drugdatatype',
                        action='store',
                        type = str,
                        choices=["GEO", "CMAP", 'L1000', 'GEO_43_fluorouracil', 'GEO_brca_fluorouracil', 'GEO_brca_doxorubicin', 'GEO_brca_paclitaxel', 'GEO_brca_docetaxel'],
                        help='Defines the drug database')
    

#obbligatorio
    parser.add_argument('p_value_limit',
                        action='store',
                        type = float,
                        choices=[0.05, 0.01, 0.001 , 0.0001],
                        help='Defines the pvalue filtering threshold for individual records')

    args = parser.parse_args()

    print(('patient_id', args.patient_id))
    print(("drugdatatype ", args.drugdatatype))
    print(('p_value_limit', args.p_value_limit))    
   
    main_function(args.patient_id, args.drugdatatype, args.p_value_limit)

###############################################
# main

if __name__ == "__main__":
    main_function_2()
