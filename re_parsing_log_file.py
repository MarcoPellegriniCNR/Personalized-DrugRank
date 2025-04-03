#
# this command takes in a a log file a modality and do a reparsing  of the file to extract useul info
#

import os

##############################################
# import file handling utils

from file_handling_util import testexistencedirectory, testexistencefile

################################################
# generic utilities

from generic_util import printmap, printlist, increment_count

#################################################

def directoytesting_with_paths(path, prefix , suffix):
    

    deeptest = False
#    deeptest = True
    
    cwd = os.getcwd()
    print(cwd)

    if os.path.exists(path):
        print(path + ' : exists')
        if os.path.isdir(path):
            print(path + ' : is a directory')
        else :
            print('Warning '  + path + ' : is not a directory')
    else :
        print('Warning ' +path + ' : does not exists')
        

    rootdir = path
    allfiles_with_paths = []
    for root, dirs, files in os.walk(rootdir):
        print("{0} has {1} files".format(root, len(files)))
        for x in files :
#            print (root ,  x)
            if x.startswith(prefix) and x.endswith(suffix) :
                allfiles_with_paths.append((root, x))
            if deeptest :
                assert testexistencedirectory(root)
                assert testexistencefile(root, x)

    return allfiles_with_paths


#####################################################


def load_lines_two_prefix(directory, filename, prefix1, prefix2) :

    existsdir = testexistencedirectory(directory)
    existsfile = testexistencefile(directory,filename)

    print( "Opening", directory, filename)

    if not existsdir or not existsfile :
        return []


    f = open(os.path.join(directory,filename), 'r')

    linenumber = 0
    useful_lines = 0
    accum = []
    line = f.readline()

    while line != '' :
        linenumber +=1
        if line.startswith(prefix1) or  line.startswith(prefix2):
                accum.append(line)
                useful_lines +=1
        line = f.readline()

    f.close()

    print ('In parserd file {}, scanned {} lines,  useful {}'.format(filename, linenumber, useful_lines))
    return accum


#####################################################


def load_lines_one_prefix(directory, filename, prefix1) :

    existsdir = testexistencedirectory(directory)
    existsfile = testexistencefile(directory,filename)

    print( "Opening", directory, filename)

    if not existsdir or not existsfile :
        return []


    f = open(os.path.join(directory,filename), 'r')

    linenumber = 0
    useful_lines = 0
    accum = []
    line = f.readline()

    while line != '' :
        linenumber +=1
        if line.startswith(prefix1) :
                accum.append(line)
                useful_lines +=1
        line = f.readline()

    f.close()

    print ('In parserd file {}, scanned {} lines,  useful {}'.format(filename, linenumber, useful_lines))
    return accum




# asumes two lines hod the data

def parse_to_extract_pval_values(in_list, verbose= 'noverbose') :

    i = 0
    accum = []
    while 2*i  + 1 < len(in_list) :

        drug_db = in_list[2*i].split(' ')
        if verbose == 'verbose' :
            print(drug_db[8])
        pvals = in_list[2*i+1].split(' ')
        if verbose == 'verbose' :
            print(pvals[3], pvals[6])
        accum.append((drug_db[8], pvals[3].strip(','), pvals[6].strip(',')))
        i+=1
    return(accum)   
        
        

def summarize_triplets(in_list, filename) :

    summary = {}

    for x in in_list :
        if float(x[1]) <= 0.05 :
            increment_count(summary, x[0], 1)
        if float(x[2]) <= 0.05 :
            increment_count(summary, x[0], 1)

    printmap(summary, 'summary for {} '.format(filename))         
             

        


    
                

############

if __name__ == "__main__":
#    pass




#    directory = '../log_files'
    directory = '.'
    
#    filename = 'log_extended_eval_TCGA_AZ_6605_2025_03_14.txt'
    prefix1 = '@@ P_value'
    prefix2 = '@@ suminv'
##    file_prefix = 'log_extended_eval_'
##    file_suffix =  '.txt'

##    file_prefix = 'log_new_subnetwork_eval'
##    file_suffix =  '21.txt'

#    file_prefix = 'log_extended_eval_TCGA_A6_2671_202#5_03_31'
#    file_prefix = 'log_extended_eval_TCGA_A6_2671_2025_04_01'

#    file_prefix = 'log_extended_eval_TCGA_AA_3517_2025_03_31'
#    file_prefix = 'log_extended_eval_TCGA_AA_3517_2025_04_01'

#    file_prefix = 'log_extended_eval_TCGA_AZ_6605_2025_03_31'

#    file_prefix = 'log_extended_eval_TCGA_AZ_6600_2025_03_31'
    file_prefix = 'log_extended_eval_TCGA_AZ_6600_2025_04_01'
    
#    file_prefix = 'log_extended_eval_TCGA_HC_7742_2025_03_31'
    
    file_suffix = 'txt'

    

    file_list = directoytesting_with_paths(directory, file_prefix, file_suffix)
    printlist(file_list, 'file_list')

    exception = ['log_extended_eval_TCGA_AZ_6600_2025_03_14.txt', 'log_extended_eval_TCGA_AZ_6600_2025_03_14_bis.txt']

    for x in file_list  :
        filename = x[1]
        if filename not in exception :
            useful_lines = load_lines_two_prefix(directory, filename, prefix1, prefix2)
            tripets = parse_to_extract_pval_values(useful_lines)
            summarize_triplets(tripets, filename)

    for x in file_list  :
        filename = x[1]
        if filename not in exception :
            useful_lines = load_lines_one_prefix(directory, filename, '@@ For ')
            printlist(useful_lines, 'usefule lines for position of therapy drugs')   
