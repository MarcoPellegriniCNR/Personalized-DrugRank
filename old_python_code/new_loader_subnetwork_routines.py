# -*- coding: cp1252 -*-
# 2025-03-19
#This file new_loading_subnetwork_routines  holds the basic data laoding code from file of subnetwrok data.

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

from new_loader_routines import resize

def load_data_subnetwork(directory, filename, separator, file_format) :

    existsdir = testexistencedirectory(directory)
    existsfile = testexistencefile(directory,filename)

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
# laod data from the configuration s input


def load_data_subnetwork_from_conf(conf, separator) :

    file_format = conf[2]
    
    print('--------------------------- Begin --------------------------------')
    if conf[2] not in ['CMAPupdown',  'L1000updown']:
        print(conf, subnetwork_data[conf])
        directory = subnetwork_data[conf][0]
        filename = subnetwork_data[conf][1]
#               counter +=1
        result_accumulator, clean_header  = load_data_subnetwork(directory, filename, separator, file_format)
    elif conf[2] == 'CMAPupdown' :
        conf1 = (conf[0], conf[1], 'CMAPup', conf[3], conf[4]   )
        
        print(conf1, subnetwork_data[conf1])
        directory = subnetwork_data[conf1][0]
        filename = subnetwork_data[conf1][1]
#               counter +=1
        result_accumulator1, clean_header = load_data_subnetwork(directory, filename, separator, file_format)

        conf2 = (conf[0], conf[1], 'CMAPdown', conf[3], conf[4] )
        print(conf1, subnetwork_data[conf2])
        directory = subnetwork_data[conf2][0]
        filename = subnetwork_data[conf2][1]
#               counter +=1
        result_accumulator2, _ = load_data_subnetwork(directory, filename, separator, file_format)
       
        result_accumulator = result_accumulator1 + result_accumulator2

    elif conf[2] == 'L1000updown' :
        conf1 = (conf[0], conf[1], 'L1000up', conf[3], conf[4]   )
        
        print(conf1, subnetwork_data[conf1])
        directory = subnetwork_data[conf1][0]
        filename = subnetwork_data[conf1][1]
#               counter +=1
        result_accumulator1, clean_header = load_data_subnetwork(directory, filename, separator, file_format)

        conf2 = (conf[0], conf[1], 'L1000down', conf[3], conf[4] )
        print(conf1, subnetwork_data[conf2])
        directory = subnetwork_data[conf2][0]
        filename = subnetwork_data[conf2][1]
#               counter +=1
        result_accumulator2, _ = load_data_subnetwork(directory, filename, separator, file_format)
       
        result_accumulator = result_accumulator1 + result_accumulator2
     
        print('--------------------------- End --------------------------------')
        
    return result_accumulator, clean_header

#####################################################
# tesing the loader routines



if __name__ == "__main__":
#    pass
    print('Testing access to all files')
    separator = ','
    counter = 1
    counter_test_limit = 100
    
#    patient_list = ['pilot', ]
#    patient_list = ['TCGA-AZ-6600', ]
#    patient_list = ['TCGA-HC-7742']
#    patient_list = ['TCGA-A6-2671']
 #   patient_list = ['TCGA-AA-3517']
#    patient_list = ['TCGA-AZ-6605',]
    patient_list = ['TCGA-A6-5667']
    alg_list = ['cpd07j08', 'MD' , 'KPM', 'Degas', 'clex3800' ,  'clex3000', ]
    data_list = [ 'GEO', 'DM' , 'CMAPup', 'CMAPdown', 'CMAPupdown', 'L1000up', 'L1000down', 'L1000updown']
#    data_list = [  'CMAPup', 'CMAPdown', 'CMAPupdown', 'L1000up', 'L1000down', 'L1000updown']
#    data_list = [  'L1000up', 'L1000down', 'L1000updown']

    
    for conf in subnetwork_data.keys() :
        if conf[0] in patient_list and  conf[1] in alg_list  and conf[2] in data_list :

            result_accumulator, clean_header = load_data_subnetwork_from_conf(conf, separator)
            if result_accumulator != [] :
                counter +=1
                
##            print('--------------------------- Begin --------------------------------')
##            if conf[2] not in ['CMAPupdown',  'L1000updown']:
##                print(conf, subnetwork_data[conf])
##                directory = subnetwork_data[conf][0]
##                filename = subnetwork_data[conf][1]
##                counter +=1
##                result_accumulator, _ = load_data_subnetwork(directory, filename, separator)
##            elif conf[2] == 'CMAPupdown' :
##                conf1 = (conf[0], conf[1], 'CMAPup', conf[3], conf[4]   )
##                
##                print(conf1, subnetwork_data[conf1])
##                directory = subnetwork_data[conf1][0]
##                filename = subnetwork_data[conf1][1]
##                counter +=1
##                result_accumulator1, _ = load_data_subnetwork(directory, filename, separator)
##
##                conf2 = (conf[0], conf[1], 'CMAPdown', conf[3], conf[4] )
##                print(conf1, subnetwork_data[conf2])
##                directory = subnetwork_data[conf2][0]
##                filename = subnetwork_data[conf2][1]
##                counter +=1
##                result_accumulator2, _ = load_data_subnetwork(directory, filename, separator)
##               
##                result_accumulator = result_accumulator1 + result_accumulator2
##
##            elif conf[2] == 'L1000updown' :
##                conf1 = (conf[0], conf[1], 'L1000up', conf[3], conf[4]   )
##                
##                print(conf1, subnetwork_data[conf1])
##                directory = subnetwork_data[conf1][0]
##                filename = subnetwork_data[conf1][1]
##                counter +=1
##                result_accumulator1, _ = load_data_subnetwork(directory, filename, separator)
##
##                conf2 = (conf[0], conf[1], 'L1000down', conf[3], conf[4] )
##                print(conf1, subnetwork_data[conf2])
##                directory = subnetwork_data[conf2][0]
##                filename = subnetwork_data[conf2][1]
##                counter +=1
##                result_accumulator2, _ = load_data_subnetwork(directory, filename, separator)
##               
##                result_accumulator = result_accumulator1 + result_accumulator2
##             
##                print('--------------------------- End --------------------------------')
            if counter > counter_test_limit :
                break

    
