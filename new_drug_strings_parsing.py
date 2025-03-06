# -*- coding: cp1252 -*-
# 2025-02-25
# This file holds code for extracting information form drug strings
# It is based on parsing techniques and must take many formats into condideration

######################################
# data location mapper

from new_mapping_modulation_data import modulation_data

################################################
# laoding form fiel routines

from new_loader_routines import load_data


####################################################
# drug names normalization routines

from file_parsing import stringnormalize2


#####################
# generic utilities

from generic_util import printmap, printlist, increment_count


#########################################
# Routine to extract drug name and drug action direction


############################################
# checking druglsits names

from random_drug_scoring_model import filter_two_lists
    
def extract_drug_name_direction(s, drugdatatype) :

#    print s

    if drugdatatype == 'GEO' :

        ll = s.split(' ')
        direction = ll[-1]
        if  direction not in ['up', 'down']:
            print ("Warning wrong direction lable", ll)

        index = 0
        while ll[index] not in ['homo', 'rattus', 'mus', 'bos', 'sus', 'macaca', 'oncorhynchus', 'myzus']:
            index +=1
            if index > 6 :
                break   
        if index == 7 :
            print (ll)
            drugname = None
        else :
            drugname = '_'.join(ll[:index]).lower()
            drugname = stringnormalize2(drugname)

        return (drugname,direction)

    elif drugdatatype in ('DM', 'DM2') :

        ll = s.split('-')
        direction = ll[-1]
        if  direction not in ['up', 'dn']:
            print ("Warning wrong direction lable", ll)
            assert False
            direction = None
        else :
            if direction == 'dn' : # correction to have uniform lables.
                direction = 'down'
        lll =  s.split(' ')
        index = 0
        while lll[index] not in ['mg/kg', 'uM', 'ng/ml']:
            index +=1
            if index > len(lll) -1 :
                break   
        if index == len(lll) :
            print (lll)
            drugname = None
        else :
            drugname1 = '_'.join(lll[:index])
            drugname2 = drugname1.split('-')
            drugname3 = drugname2[:-1] # remove the dose
            drugname4 = '_'.join(drugname3)
            drugname5 = drugname4.split(' ')
            drugname =  '_'.join(drugname5).lower()
            drugname = stringnormalize2(drugname)

#        print 'checking : ' , s, (drugname,direction)
        return (drugname,direction)

    elif drugdatatype == 'CMAP' :
        ll = s.split('-')
        direction = ll[-1]
        if  direction not in ['up', 'dn']:
            print ("Warning wrong direction lable", ll)
            direction = None
        else :
            if direction == 'dn' : # correction to have uniform lables.
                direction = 'down'
        drugname = '_'.join(ll[:-2]).lower()
        drugname = stringnormalize2(drugname)
        return (drugname,direction)

    elif drugdatatype == 'L1000' :
        s1 = s.split(' ')
        s2 = ' '.join(s1[2:])
        ll = s2.split('-')
        direction = ll[-1]
        if  direction not in ['up', 'dn']:
            print ("Warning wrong direction lable", ll)
            direction = None
        else :
            if direction == 'dn' : # correction to have uniform lables.
                direction = 'down'
        l2 = ll[:-2]
        l3 = l2[1:]
        l4 = '-'.join(l3)
        drugname= stringnormalize2(l4)
        if drugname == '' :
            print ("Warning : empty drug name : ", (drugname,direction), s1,s2, ll, l3)
            assert False
#        print s, (drugname,direction)
        return (drugname,direction)

    else :
        print ("Wrong drugdatatype : ", drugdatatype)
        assert False


############################################################
#
# a function to cehck the presence of a given drug in an input file taking alises into consideration
#




#####################################################
# tesing the directionality score  computations


if __name__ == "__main__":
#    pass
    print('Testing drug string parsing routines')
    separator = ','
    counter = 1
    counter_test_limit = 16

#    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
    patient_list = ['TCGA-HC-7742']
    alg_list = ['cpd07j08', 'MD' ]

##    drugs_to_check = [  'Fluorouracil',
##                        'Leucovorin',
##                        'Oxaliplatin',
##                        'Bevacizumab',
##                        'Capecitabine',
##                        'Irinotecan',
##                        'Cetuximab',
##                        'Pegfilgrastim',
##                        'Mitomycin',]
##    
##    drugs_to_check_lable = 'crc_most_used_drugs'

    drugs_to_check = ['Leuprolide', 'Leuprorelin', 'Lupron', 'Eligard', 'Fensolvi', 'Viadur', 'Lupaneta' ]
    drugs_to_check_lable = 'prc_drug_used_for_TCGA-HC-7742'
    

    for conf in modulation_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list :
            print('--------------------------- Begin --------------------------------')
            print(conf, modulation_data[conf])
            directory = modulation_data[conf][0]
            filename = modulation_data[conf][1]
            drug_format = conf[2]

            counter +=1
            result_accumulator, _ = load_data (directory, filename, separator)

    # parsing all the drug sttrings
            unique_drugs = {}
            unique_direction = {}
            for x in result_accumulator :
              result =  extract_drug_name_direction(x[1], drug_format)
    #          print(result)
              increment_count(unique_drugs, result[0], 1)
              increment_count(unique_direction, result[1], 1)

            drug_names_in_matrix = [y for y in  unique_drugs.keys()] 
            printmap(unique_drugs, 'unique_drugs')
            printmap(unique_direction, 'unique_direction')
            filter_two_lists(drugs_to_check_lable, drugs_to_check, drug_format, drug_names_in_matrix )
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break
