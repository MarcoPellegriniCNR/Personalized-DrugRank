# -*- coding: cp1252 -*-
# 2025-03-19
# this file check the parsing of  drugns rings in the subnetwork data


import os

################################################
# main mapping parameter to location

#from new_mapping_modulation_data import modulation_data
from new_mapping_subnetwork_data import  subnetwork_data

################################################
# generic utilities

from generic_util import printmap, printlist, increment_count

##############################################
# import file handling utils

from file_handling_util import testexistencedirectory, testexistencefile


####################################
# loading function and aprsinf functions

from new_loader_subnetwork_routines import load_data_subnetwork, load_data_subnetwork_from_conf

# from new_drug_strings_parsing import extract_drug_name_direction

#######################################################

from file_parsing import subparse2, subparse3, extractdrugstring

def extract_drug_name_only(instring, file_format) :

    if  file_format in ['L1000up', 'L1000down' , 'L1000updown'] :
        drugstring = subparse2(instring)

    elif file_format in ['CMAPup', 'CMAPdown' , 'CMAPupdown'] :
        drugstring = subparse3(instring)

    elif file_format in ['GEO', 'DM'] :
        drugstring = extractdrugstring(instring, file_format)

    else :
        print('Wrong option in extract_drug_name_only', file_format)
        assert False

    return drugstring
        

    


#####################################################
# tesing the loader routines

#####################################################
# tesing the loader routines



if __name__ == "__main__":
#    pass
    print('Testing access to all files')
    separator = ','
    counter = 1
    counter_test_limit = 15
    
#    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
#    patient_list = ['TCGA-HC-7742']
#    patient_list = ['TCGA-A6-2671']
#    patient_list = ['TCGA-AA-3517']
    patient_list = ['TCGA-AZ-6605',]
    alg_list = ['cpd07j08', 'MD' , 'KPM', 'Degas', 'clex3800' ]
#    data_list = [ 'GEO', 'DM']
#    data_list = [  'CMAPup', 'CMAPdown', 'CMAPupdown', 'L1000up', 'L1000down', 'L1000updown']
#    data_list = [  'CMAPup', 'CMAPdown', 'CMAPupdown', 'L1000up', 'L1000down', 'L1000updown',  'GEO', 'DM']
    data_list = [ 'L1000up', 'L1000down', 'L1000updown']
    
    for conf in subnetwork_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list  and conf[2] in data_list :
            print('--------------------------- Begin --------------------------------')
##            print(conf, subnetwork_data[conf])
##            directory = subnetwork_data[conf][0]
##            filename = subnetwork_data[conf][1]
            drug_format = conf[2]
            
            counter +=1
#            result_accumulator, _ = load_data_subnetwork(directory, filename, separator)
            result_accumulator, _ = load_data_subnetwork_from_conf(conf, separator)

                # parsing all the drug sttrings
            unique_drugs = {}
            unique_direction = {}
            for x in result_accumulator :
              result =  extract_drug_name_only(x[0], drug_format)  # for the subnetwork data the drung names are at  index 0  
#              print(x[0], result)
              increment_count(unique_drugs, result, 1)
#              increment_count(unique_direction, result[1], 1)

            drug_names_in_matrix = [y for y in  unique_drugs.keys()] 
            printmap(unique_drugs, 'unique_drugs')
#            printmap(unique_direction, 'unique_direction')

# cehcking presence of specific drugs            
#            filter_two_lists(drugs_to_check_lable, drugs_to_check, drug_format, drug_names_in_matrix )


            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break


