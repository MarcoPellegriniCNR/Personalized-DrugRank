# -*- coding: cp1252 -*-
# 2025-02-24
# This file new_loading_procedures holds the basic data laoding code from file.
# It is designed to be robust aginst possible non-uniformity of the data
# in input due tothe possible use of the separator character in the drug identifiers.

import os

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

def load_data (directory, filename, separator) :

    existsdir = testexistencedirectory(directory)
    existsfile = testexistencefile(directory,filename)

    print( "Opening", directory, filename)

    if not existsdir or not existsfile :
        return []

    

    f = open(os.path.join(directory,filename))
    accumulator = []
    fields_counter = {}
    recovered = 0
    normalsize = 0
    linenumber = 0
    testprintnumber = 0


    line = f.readline()
    header = line.split(separator)
    print ("Header : ", len(header), header)
    clean_header = [x.strip('\n').strip('"') for x in header ]
    print ('Clean header', len(clean_header), clean_header)
    numfields = len(header)
    line = f.readline()
    while line != '' :
        linenumber +=1
#        print line
        record =  line.split(separator)
        increment_count(fields_counter, len(record))
        delta = len(record) - numfields

        if delta == 0 :
            normalsize +=1
            record1 = [x.strip('\n').strip('"') for x in record ]
#        print record1
            accumulator.append(record1)
        elif delta > 0 :
#            print "Warning: separator problem: ", len(record), record
            newrecord = resize(record, delta , '_')
            assert len(newrecord) == numfields
            recovered += 1
            if testprintnumber < 10 :
                print ("sample ", testprintnumber, newrecord)
                testprintnumber += 1
            record1 = [x.strip('\n').strip('"') for x in newrecord ]
            accumulator.append(record1)
        else :
            print ("Warning Too few seprators ", record)

            
        line = f.readline()
        
    f.close()
    printmap(fields_counter, "distribution of filds number by separator :"+separator)
    print ("normalsize ", normalsize, "recovered", recovered, "matching total", recovered+normalsize , linenumber)

    return accumulator, clean_header


#####################################################
# tesing the loader routines



if __name__ == "__main__":
#    pass
    print('Testing access to all files')
    separator = ','
    counter = 1
    counter_test_limit = 16
    
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
            counter +=1
            result_accumulator, _ = load_data (directory, filename, separator)
            print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break

    


##
