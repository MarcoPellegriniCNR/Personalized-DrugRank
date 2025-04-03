# -*- coding: cp1252 -*-
# 2025-02-25
# This file holds the code for computing the direcionality score  of data uploaded
# from file.
# We introduce filtering routines based on pvaluees
# We introduce variants due to the strucure of th record in filter/nofiltering files as input

################################################
# main mapping parameter to location

from new_mapping_modulation_data import modulation_data

################################################
# laoding form fiel routines

from new_loader_routines import load_data


################################################
# parse stingfro drugnam and direction of action
#

from new_drug_strings_parsing import extract_drug_name_direction

###################################
# genric utilities

from generic_util import printlist, printmap, domain_restriction, minuslog10, signedlog10

##################################
# loaign name filtering routines

from random_drug_scoring_model import filterbyfdaapproveddrugs

#   newrecord = (module_id, res[0], res[1], module_size, disease_up, disease_down, drug_deg_in_module, drug_deg_hits)

record_field_map = {
    'modules' : 0 ,
    'drug_name' : 1 ,
    'drug_direction' : 2,
    'module_size' : 3,
    'DEG_up' : 4,
    'DEG_down' : 5,
    'drug_module_genes' : 6,
    'drug_deg_hits' : 7,
    'p_val_enrichr' : 8,
    }

#######################
# main function exported

##############################################
# postiveonly = True if we report only drugs with postive correaaltion
# fdaonly = True  if we filter the drugs retainign only theos eapproved by FDA
# drugdatatyep in: GEO, DM, CMAP and L1000
# weights in : "size_x_deltadeg", "deltadeg"

def direction_analysis(directory, filename, separator, drugdatatype, positiveonly, fdaonly ,pvalue_filter,p_value_limit , weights):

    assert weights in ["size_x_deltadeg", "deltadeg", "pval_x_deltadeg", "pval_x_log_deltadeg", 'pval_x_drug_hits', 'pval_x_deltadeg_x_drug_hits']
    assert drugdatatype in ['GEO', 'DM', 'CMAP', 'L1000']
#    assert not filename.endswith('noFiltering.csv') #code is only for the format at 9 fields in the filtered data 
        
    dataarray , clean_header = load_data(directory, filename, separator)              
    print ("Number of records ", len(dataarray))
    print ('First raw : ' , dataarray[0])
#    assert len(dataarray[0]) == 9 #code is only for the foramt at 9 fields
    clean_header_map = {clean_header[x] : x for x in range(len(clean_header))}
    printmap(clean_header_map, 'clean_header_map') # maps field names to field indices spcific for the filter files. 

# p-value filter

#    pval_filtered_dataarray = [x for x in dataarray if test_pval(x[7], 0.05)] # use the raw pval not the adjusted one
    if pvalue_filter :
        pval_filtered_dataarray = [x for x in dataarray if test_pval(x[clean_header_map['p_value']], p_value_limit)] # use the raw p_value not the adjusted one 
        print ("Number of records after p_value filtering ", len(pval_filtered_dataarray))
    else:
       pval_filtered_dataarray =   dataarray
       print('No p_value filtering')


# remapping tupels after extraction of drung names and direction

    drug_module_map = {}
    drugnames = set()

    count_uniform = 0
    count_mixed = 0

    for x in pval_filtered_dataarray :

##        res = extract_drug_name_direction(x[1], drugdatatype)
##        module_id = int(x[0])
##        module_size = int(x[2])
##        disease_up = int(x[3])
##        disease_down = int(x[4])
##        drug_deg_in_module = int(x[5])
##        drug_deg_hits = int(x[6])
        
        (drug_name, drug_direction) = extract_drug_name_direction(x[clean_header_map['drugs']], drugdatatype)
        module_id = int(x[clean_header_map['modules']])
        module_size = int(x[clean_header_map['module_size']])
        disease_up = int(x[clean_header_map['DEG_up']])
        disease_down = int(x[clean_header_map['DEG_down']])
        drug_deg_in_module = int(x[clean_header_map['drug_module_genes']])
        p_val_enrichr = float(x[clean_header_map['p_value']])

# handles difference veteeen filter/nofilter intu files
        if 'DE_drug' in clean_header_map.keys() :
            drug_deg_hits = int(x[clean_header_map['DE_drug']])
        elif 'DE_up_drug' in clean_header_map.keys() and 'DE_down_drug' in clean_header_map.keys() :
 #           drug_deg_hits = int(x[clean_header_map['DE_up_drug']]) +  int(x[clean_header_map['DE_down_drug']])
            if  drug_direction == 'up' :
                drug_deg_hits = int(x[clean_header_map['DE_down_drug']]) - int(x[clean_header_map['DE_up_drug']])
            elif drug_direction == 'down' :
                drug_deg_hits = int(x[clean_header_map['DE_up_drug']]) - int(x[clean_header_map['DE_down_drug']])
            else :
                print('Wrong drug_direction lable', drug_direction)
                assert  False
        else:
            print("Wrong labels in mapping ")
            assert False

        newrecord = (module_id, drug_name, drug_direction, module_size, disease_up, disease_down, drug_deg_in_module, drug_deg_hits, p_val_enrichr)
        # new represetnation with 7 fields remove the pvalues and add the drug name drug direction, use with the record_field_map to index

    # check if any module hd both up and down DEG for the disease.relevant for the nofilter case

        if  disease_up != 0 and disease_down != 0 :
    #        print newrecord
            count_mixed += 1
        else :
            count_uniform += 1

# collect all recods for the same drug on teh same module - intended to average later across the dosages

        if (module_id, drug_name) not in   drug_module_map.keys() :        
            drug_module_map[(module_id, drug_name)] = [newrecord]
        else :
            drug_module_map[(module_id, drug_name)].append(newrecord)

# collect unique drug names

        drugnames.add(drug_name)

    print ("Drug-module action mixed : ", count_mixed, "uniform : " , count_uniform)
    print ("Number of drugs " , len(drugnames))

    # initialize mappings to hold the reorganized results by drugname
    drugscore = {}
    drug_modules = {}
    for x in drugnames :
        drugscore[x] = 0
        drug_modules[x] = []


    for x in drug_module_map.keys() :
#        drugdirection = drug_module_map[x][0][2] #look at the first record only
        drugname = x[1]

        drugscore[drugname] +=  score(drug_module_map[x], 'deg_in_module', weights)
        drug_modules[drugname].append(x[0]) # record the module id

    drugscore_items = drugscore.items()
    drugscore_items_sorted = sorted(drugscore_items, key=lambda tup: tup[1], reverse=True)
    printlist(drugscore_items_sorted, 'drugscore_items_sorted', 20)

# fda filter
    if fdaonly == True :
        list_of_drug_names = [x for x in drugscore.keys()]
        fda_filtered_list_of_drug_names = filterbyfdaapproveddrugs(list_of_drug_names, drugdatatype)
        drugscore_fda_filtered = domain_restriction(fda_filtered_list_of_drug_names, drugscore)
        drugscore_items = drugscore_fda_filtered.items()
        drugscore_items_sorted = sorted(drugscore_items, key=lambda tup: tup[1], reverse=True)
        printlist(drugscore_items_sorted, 'drugscore_items_sorted after fda filtering', 20)

# positivity filter

    if  positiveonly == True :
        drugscore_items_sorted_tmp = [x for x in drugscore_items_sorted if x[1] >= 0.0 ]
        drugscore_items_sorted = drugscore_items_sorted_tmp
        printlist(drugscore_items_sorted, 'drugscore_items_sorted after positivity filtering', 20)


    return drugscore_items_sorted


##############################################################
# auxiliary functions


#######################
def test_pval(s, pval_limit) :

    if s == None :
        return False
    if s in ['', 'NA', 'na'] :
        return False
    try:
        x = float(s)
#        print (x)
        return x <= pval_limit
    except ValueError:
        print(s, "Not a float")


##########################

###############################################
# score a records for counteraction drug-disease


def scoreonerecord(record, typescore, weights) :

##    drugdirection = record[2]
##    upcount = record[4]
##    downcount = record[5]
    drugdirection = record[record_field_map['drug_direction']]
    upcount = record[record_field_map['DEG_up']]
    downcount = record[record_field_map['DEG_down']]
    p_val_enrichr = record[record_field_map['p_val_enrichr']]
    drug_deg_hits = record[record_field_map['drug_deg_hits']]

    if typescore == 'deg_in_module' :
#        drug_in_modue = record[6]
        drug_in_modue = record[record_field_map['drug_module_genes']]
    elif typescore == 'hits_in_module':
#        drug_in_modue = record[7]
        drug_in_modue = record[record_field_map['drug_deg_hits']]
    else :
        print ("wrong typescore ", typescore)
        assert False

    if weights not in ["size_x_deltadeg", "deltadeg", "pval_x_deltadeg", "pval_x_log_deltadeg", 'pval_x_drug_hits', 'pval_x_deltadeg_x_drug_hits'] :
        print ("Wrong weights ", weights)
        assert False

    sc = 0 

    if     drugdirection == 'down':
        if weights== "size_x_deltadeg":
            sc =  drug_in_modue * (upcount - downcount)
        elif weights== "pval_x_deltadeg" :
            sc= minuslog10(p_val_enrichr) * (upcount - downcount)
        elif weights==  "pval_x_log_deltadeg" :
            sc = minuslog10(p_val_enrichr) *signedlog10(upcount - downcount)
        elif weights==  "pval_x_drug_hits" :
            sc = minuslog10(p_val_enrichr) * drug_deg_hits
        elif weights == 'pval_x_deltadeg_x_drug_hits' :
            if  upcount - downcount > 0 and drug_deg_hits > 0 :
                sc = minuslog10(p_val_enrichr) * (upcount - downcount)*drug_deg_hits
            elif upcount - downcount < 0 and drug_deg_hits < 0 :
                sc = minuslog10(p_val_enrichr) * (upcount - downcount)*abs(drug_deg_hits)
            else:
                sc = 0        
        else :
            sc = upcount - downcount

    elif   drugdirection == 'up':
        if weights== "size_x_deltadeg":
            sc = drug_in_modue * (downcount -  upcount)
        elif weights== "pval_x_deltadeg" :
             sc =  minuslog10(p_val_enrichr) *(downcount -  upcount)
        elif weights==  "pval_x_log_deltadeg" :
             sc = minuslog10(p_val_enrichr) *signedlog10(downcount -  upcount)
        elif weights==  "pval_x_drug_hits" :
            sc = minuslog10(p_val_enrichr) * drug_deg_hits
        elif weights == 'pval_x_deltadeg_x_drug_hits' :
            if  downcount -  upcount > 0 and drug_deg_hits > 0 :
                sc = minuslog10(p_val_enrichr) * (downcount -  upcount)*drug_deg_hits
            elif downcount -  upcount < 0 and drug_deg_hits < 0 :
                sc = minuslog10(p_val_enrichr) * (downcount -  upcount)*abs(drug_deg_hits)
            else:
                sc = 0       
        else :
            sc = downcount -  upcount

    else :
        print ("format error ", record)
        assert False

    return sc

#######################################
# sum all scores over a list or records and take the average velue 

def score(list_of_records, typescore, weights) :

    assert len(list_of_records) > 0
    sc = 0
    for x in list_of_records :
        sc += scoreonerecord(x, typescore, weights)

    return float(sc)/float(len(list_of_records))
    


##############
# scitpting

#####################################################
# tesing the directionality score  computations


if __name__ == "__main__":
#    pass
    print('Testing computation of the directionality score')

    separator = ','
    counter = 1
    counter_test_limit = 8
    
    #    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
    patient_list = ['TCGA-HC-7742']
    alg_list = ['cpd07j08', 'MD' ]


    for conf in modulation_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list :
            print('--------------------------- Begin --------------------------------')
            print(conf, modulation_data[conf])
            directory = modulation_data[conf][0]
            filename = modulation_data[conf][1]
            drug_format = conf[2]
            positiveonly = True #per ora  usato
            fdaonly = True  #per ora  usato
            pvalue_filter =True
            p_value_limit = 0.01
    #        weights = "deltadeg"
     #       weights = "size_x_deltadeg"
            weights = "pval_x_deltadeg"
            
            

    #        if not filename.endswith('noFiltering.csv') and not filename.endswith('noF.csv'):
            if conf[3] in ['filter', 'nofilter']  :
                counter +=1
                result_accumulator = direction_analysis(directory, filename, separator, drug_format, positiveonly,  fdaonly , pvalue_filter, p_value_limit , weights)
            else :
                print('No code for this file format yet', filename)
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break
