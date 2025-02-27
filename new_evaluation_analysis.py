# -*- coding: cp1252 -*-
# 2025-02-27
# This file holds the code forevaluating the direcionality score
# We introduce filtering routines based on pvaluees
# We introduce variants due to the stricure of th record in filter/nofiltering files

################################################
# main mapping parameter to location

from new_mapping_modulation_data import modulation_data

################################################
# laoding form fiel routines

from new_loader_routines import load_data

#############################################
# import the directionality score  routine

from new_directionality_analysis import direction_analysis

#####################################
# importing TTD data

from random_drug_scoring_model import  filterbycolorectcandata


###########################
# main function

# upload the TTD data


def new_ranking_evaluation_via_TTD_data(conf, initial_list) :

    drugnames = [x[0] for x in initial_list]
    
    trial_drug_list = filterbycolorectcandata(drugnames, conf)

    compute_performance_simplified(initial_list, trial_drug_list, conf)

    


###########################
# auxilairy functions


def compute_performance_simplified(initiallist, trialdruglist, conf):

    print ("***********************************************")

    suminvrserank = 0
    precisionat20 = 0
    rankedlist = []
    accum = [0,0,0,0]
    for x in range(len(initiallist)) :
        drugname = initiallist[x][0] # drug name in first postion in the record                    
        if drugname in trialdruglist :
            print (x, drugname, 1.0/float(x+1))
            rankedlist.append(( x, drugname, 1.0/float(x+1)))
            suminvrserank += 1.0/float(x+1)
            if x <= 19 :
                precisionat20 += 1
        else :
            print (x,  drugname, "no trial/TTD_Data ")
            rankedlist.append((x,  drugname, "no trial "))

    print ("suminvrserank (RHR), prec-at-20 : ", conf, suminvrserank, precisionat20)
    return (conf , suminvrserank, precisionat20, accum, rankedlist)



#####################################################
# tesing the evaluation of the directionality score  computations


if __name__ == "__main__":
#    pass
    print('Testing computation of the directionality score')

    separator = ','
    counter = 1
    counter_test_limit = 16

    #    patient_list = ['pilot', ]
    patient_list = ['TCGA-AZ-6600', ]

    for conf in modulation_data.keys() :
        if conf[0] in patient_list  :
            print('--------------------------- Begin --------------------------------')
            print(conf, modulation_data[conf])
            directory = modulation_data[conf][0]
            filename = modulation_data[conf][1]
            drug_format = conf[2]
            positiveonly = True #per ora non usato
            fdaonly = True  #per ora non usato
    #        weights = "deltadeg"
            weights = "size_x_deltadeg"
            
            

    #        if not filename.endswith('noFiltering.csv') and not filename.endswith('noF.csv'):
            if conf[3] in ['filter', 'nofilter']  :
                counter +=1
                result_accumulator = direction_analysis(directory, filename, separator, drug_format, positiveonly,  fdaonly ,weights)
                new_ranking_evaluation_via_TTD_data(conf, result_accumulator) 
            else :
                print('No code for this file format yet', filename)
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break
