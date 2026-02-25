# -*- coding: cp1252 -*-
# 2025-02-24
# This file new_loading_procedures holds the basic data laoding code from file.
# It is designed to be robust aginst possible non-uniformity of the data
# in input due tothe possible use of the separator character in the drug identifiers.

import os
import sys

###########################################
# Accessing functions in siblig directories

path = os.path.join("..")
print("Relative path exists", os.path.isdir(path))
abspath = os.path.abspath(path)
print("Absolute path exists", os.path.isdir(abspath))
sys.path.append(abspath)
print("Including in sys path to ", abspath)


################################################
# main mapping parameter to location

from new_mapping_modulation_data import modulation_data


################################################
# generic utilities

from generic_util import printmap, printlist, increment_count

##############################################
# import file handling utils

from file_handling_util import testexistencedirectory, testexistencefile


###############################################################
#
# new version of load_data


######################################
        
# merg the k eccess filds in position from 1 to 1+k and resize the tuple:
# change the separator in the process replacing the standard separator with undescrore in names of drugs.
#

def resize(inlist, k , new_separator):

    assert k >= 0
    assert k <= len(inlist) -1
    new_field_one = new_separator.join(inlist[1:2+k])
    outlist = [inlist[0],] + [new_field_one,] + inlist[2+k:]
    return outlist

####################
#
# If the separator is used withn the drug names it recovers the situation by using underscrore in teh drug name

def load_data (directory, filename, separator, verbose='verbose') :

    existsdir = testexistencedirectory(directory, verbose)
    existsfile = testexistencefile(directory,filename, verbose)

    print( "Opening", directory, filename)

    if not existsdir or not existsfile :
        return [], []

    

    f = open(os.path.join(directory,filename))
    accumulator = []
    fields_counter = {}
    recovered = 0
    normalsize = 0
    linenumber = 0
    testprintnumber = 0


    line = f.readline()
    header = line.split(separator)
    if verbose=='verbose' :
        print ("Header : ", len(header), header)
    clean_header = [x.strip('\n').strip('"') for x in header ]
    if verbose=='verbose' :
        print ('Clean header', len(clean_header), clean_header)
    numfields = len(header)
    line = f.readline()
    while line != '' :
        linenumber +=1
#        print line
        record =  line.split(separator)
        increment_count(fields_counter, len(record), 1)
        delta = len(record) - numfields

        if delta == 0 :
            normalsize +=1
            record1 = [x.strip('\n').strip('"') for x in record ]
            if verbose=='veryverbose' :
                print (record1)
            accumulator.append(record1)
        elif delta > 0 :
#            print "Warning: separator problem: ", len(record), record
            newrecord = resize(record, delta , '_')
            assert len(newrecord) == numfields
            recovered += 1
            if testprintnumber < 10 :
#                print ("sample ", testprintnumber, newrecord)
                testprintnumber += 1
            record1 = [x.strip('\n').strip('"') for x in newrecord ]
            accumulator.append(record1)
        else :
            print ("Warning Too few seprators ", record)

            
        line = f.readline()
        
    f.close()
    if verbose=='verbose' :
        printmap(fields_counter, "distribution of filds number by separator :"+separator)
    print ("normalsize ", normalsize, "recovered", recovered, "matching total", recovered+normalsize , linenumber)

    return accumulator, clean_header


################################
# counting NA fileds in pvalue and adjusted pvalue

def counting_NA_values(conf, result_accumulator, header) :

    pval_index = header.index('p_value')
    apval_index  = header.index('adjusted_pvalue')

    total = len(result_accumulator)
    pval_na_counter = 0
    apval_na_counter= 0
    for x in  result_accumulator :
        if x[pval_index] == 'NA' :
            pval_na_counter += 1
        if x[apval_index] == 'NA' :
            apval_na_counter += 1

    print ('@@ For conf {} over {} lines , pval NA {}, adj pval NA {},  differences {} {}'.format(conf, total, pval_na_counter, apval_na_counter,  total - pval_na_counter,  total -  apval_na_counter))
        
            
import statistics

def compare_two_listing(string1, listing1, string2, listing2) :

    headers =     [
                        'modules', 
                        'drugs', 
                        'module_size', 
                        'DEG_up', 
                        'DEG_down', 
                        'drug_module_genes', 
                        'DE_up_drug', 
                        'DE_down_drug', 
                        'p_value', 
                        'adjusted_pvalue']


    deg_up_index = headers.index('DEG_up')
    deg_down_index = headers.index('DEG_down')
    
    l1 = len(listing1)
    l2 = len(listing2)
    if  l1  != l2 :
        print('Different lengths', l1,  l2)
        return None
    counter = {}
    counter_first_NA = {}
    counter_second_NA = {}
    counter_both_not_NA = {}
    pval_delta_filtered = []
    pval_delta = []
    counter_deg_up_diff = 0
    counter_deg_up_equal = 0
    counter_deg_down_diff = 0
    counter_deg_down_equal = 0
    counter_deg_both_equal = 0
    counter_deg_one_different = 0
    
    
    for field in headers :
        counter[field] = 0
        counter_first_NA[field] = 0
        counter_second_NA[field] = 0
        counter_both_not_NA[field] = 0
        

    for field in headers :
        index = headers.index(field)
        for x in range(l1) :
            if listing1[x][index] ==  listing2[x][index] :
                counter[field]  +=1
            elif listing1[x][index] != 'NA' and listing2[x][index] == 'NA' :
                counter_second_NA[field] += 1
            elif listing1[x][index] == 'NA' and listing2[x][index] != 'NA' :
                counter_first_NA[field] += 1
            else :
                counter_both_not_NA[field] += 1
                if field == 'p_value' :
                    pval_delta.append(abs(float(listing1[x][index]) - float(listing2[x][index])))
                    if float(listing1[x][index]) <= 0.05 and float(listing2[x][index]) <=  0.005 :
                        pval_delta_filtered.append(abs(float(listing1[x][index]) - float(listing2[x][index])))
                    if listing1[x][deg_up_index] != listing2[x][deg_up_index] :
                        counter_deg_up_diff  +=1
                    else :
                        counter_deg_up_equal +=1
                        
                    if listing1[x][deg_down_index] != listing2[x][deg_down_index] :
                        counter_deg_down_diff  +=1
                    else :
                        counter_deg_down_equal +=1
                        
                    if (listing1[x][deg_up_index] == listing2[x][deg_up_index]) and (listing1[x][deg_down_index] == listing2[x][deg_down_index]) :
                        counter_deg_both_equal +=1
                    else :
                        counter_deg_one_different +=1 
                        
                        
                    
                
                
                

    for field in headers :            
        print ("For {},{} for field {},  equal slots {} different slots {}, first_NA {}, second NA {}, both not NA {}".format(string1, string2, field, counter[field] , l1- counter[field] , counter_first_NA[field], counter_second_NA[field], counter_both_not_NA[field] ))

    print('@@ For  {}, {},  size {}, Min pval delta {}, max pval_delta {}, median {}, mean {} '.format(string1, string2,  len(pval_delta), my_min(pval_delta), my_max(pval_delta), my_median(pval_delta), my_mean(pval_delta)))
    print('&& For Filtered  {}, {}, size {}, Min pval delta {}, max pval_delta {}, median {}, mean {} '.format(string1, string2, len(pval_delta_filtered) , my_min(pval_delta_filtered), my_max(pval_delta_filtered), my_median(pval_delta_filtered), my_mean(pval_delta_filtered)))
    print('%% For {},{} with unequal pval , equal up {}  diff up {}, equal down {} diff down {}, both equal {}, one different {} '.format(string1, string2, counter_deg_up_equal, counter_deg_up_diff , \
                                                                                                                                                                                                                counter_deg_down_equal, counter_deg_down_diff, counter_deg_both_equal, counter_deg_one_different))
    return None

###################

def my_min(lista) :
    if len(lista) == 0 :
        return None
    else :
        return min(lista)

def my_max(lista) :
    if len(lista) == 0 :
        return None
    else :
        return max(lista)

def my_mean(lista) :
    if len(lista) == 0 :
        return None
    else :
        return statistics.mean(lista)

def my_median(lista) :
    if len(lista) == 0 :
        return None
    else :
        return statistics.median(lista)
    
  

###################
def testing_identity_of_fields(pat1, patient_list):
    
    separator = ','

#    patient_list = ['TCGA-AZ-6600',  'TCGA-A6-5667', 'TCGA-AZ-6605', 'TCGA-AA-3517',   'TCGA-A6-2671', 'pilot' ]

    
#    pat1 = 'TCGA-A6-2671'
#  alg = 'MD'
    alg = 'cpd07j08'

    drug_data = 'GEO'
    
    conf1 = (pat1, alg , drug_data, 'nofilter', 'crc')
    index1 = patient_list.index(pat1)

    directory1 = modulation_data[conf1][0]
    filename1 = modulation_data[conf1][1]
    result_accumulator1, header1  = load_data (directory1, filename1, separator)


    for  x in patient_list :
        index2  = patient_list.index(x) 
        if index2 > index1  :                
            conf2 = (x, alg , drug_data, 'nofilter', 'crc') 
            directory2 = modulation_data[conf2][0]
            filename2 = modulation_data[conf2][1]
            result_accumulator2, header2  = load_data (directory2, filename2, separator)    
            compare_two_listing(str(conf1), result_accumulator1, str(conf2), result_accumulator2)
        

##################

def testing_access_to_all_files(patient_list, alg_list):
    
    separator = ','
    counter = 1
    counter_test_limit = 48
    
#    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
#    patient_list = ['TCGA-HC-7742']
#    patient_list = ['TCGA-A6-2671']
#    patient_list = ['TCGA-AA-3517']
#    patient_list = ['TCGA-AZ-6605',]
#    patient_list = ['TCGA-A6-5667',]

#    patient_list = ['TCGA-AZ-6600',  'TCGA-A6-5667', 'TCGA-AZ-6605', 'TCGA-AA-3517',   'TCGA-A6-2671', 'pilot' ]
    
#    alg_list = ['cpd07j08', 'MD' ]

    print(patient_list, alg_list)
    
    for conf in modulation_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list  and conf[3] == 'nofilter' :
            print('--------------------------- Begin --------------------------------')
            print(conf, modulation_data[conf])
            directory = modulation_data[conf][0]
            filename = modulation_data[conf][1]
            counter +=1
            result_accumulator, header  = load_data (directory, filename, separator)
            counting_NA_values(conf, result_accumulator, header)
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break



#####################################################
# tesing the loader routines



if __name__ == "__main__":
#    pass

#    patient_list = ['TCGA-AZ-6600',  'TCGA-A6-5667', 'TCGA-AZ-6605', 'TCGA-AA-3517',   'TCGA-A6-2671', 'pilot' ]
#    patient_list = ['TCGA-AZ-6600',  'TCGA-A6-2671']
    patient_list = ['TCGA-AZ-6600-bis','TCGA-AZ-6600']
    alg_list = ['cpd07j08',]

    print('**********************************************************************************')
    print('Testing access to all files')
    testing_access_to_all_files(patient_list , alg_list)

##    print('**********************************************************************************')
##    print('Testing identity of fields in pairs of files')
##    for pat in patient_list :
##        testing_identity_of_fields(pat, patient_list )
    
    
    

    


##
