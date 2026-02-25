#
# this is code to load and analysze many patients datawith a a common generic background
#i it performs balanced filtering of up and down regualted DEG

import os
import os.path
import math

from generic_util import printmap, printlist, minuslog10, domain_restriction, venndiagramsets

from file_handling_util import testexistencedirectory, testexistencefile, testexistenceanscreatedirectory

from new_loader_routines_modulation import load_gene_list

from new_loader_routines import load_data

# function for gene id mappings

import biorosetta as br
idmap = br.IDMapper('all')

###############
def  get_patient_name(instring):

    s0 =  instring.split('_')
    s = s0[-1].split('-')
    assert len(s) == 4 
    return '-'.join(s[:3])

###############
# setting the paths for input and output files. 

# handling CRC irinotecan aptients
#path_in =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp_fluorouracil_patients'
#path_out = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2_fluorouracil_patients'

# handling CRC florouracil _no_irinotecan patients
##path_in =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp'
##path_out = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp2'

#######
# handling BRCA florouracil patients
path_in_BRCA_florouracil =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp_fluorouracil_patients'
path_out_BRCA_florouracil = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_fluorouracil_patients'

# handling BRCA patients Doxorubicin Hydrochloride
path_in_BRCA_Doxorubicin  =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_doxorub_hyd'
path_out_BRCA_Doxorubicin = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_doxorub_hyd_patients'

# handling BRCA patients Paclitaxel
path_in_BRCA_Paclitaxel  =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_paclitaxel'
path_out_BRCA_Paclitaxel = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_paclitaxel_patients'


# handling BRCA patients Docetaxel
path_in_BRCA_Docetaxel =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_docetaxel'
path_out_BRCA_Docetaxel = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_docetaxel_patients'

# handling BRCA patients Cyclophosphamide
path_in_BRCA_Cyclophosphamide =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_cycloph'
path_out_BRCA_Cyclophosphamide = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Breast_cancer_sp2_cyclophosphamide_patients'


# cehcking intersection of lsits of patients

##directories_fl = sorted([d for d in os.listdir(path_in_BRCA_florouracil) if os.path.isdir(os.path.join(path_in_BRCA_florouracil, d)) and d.startswith('TCGA')])
##directories_dox = sorted([d for d in os.listdir(path_in_BRCA_Doxorubicin) if os.path.isdir(os.path.join(path_in_BRCA_Doxorubicin, d)) and d.startswith('TCGA')])
##directories_pac = sorted([d for d in os.listdir(path_in_BRCA_Paclitaxel) if os.path.isdir(os.path.join(path_in_BRCA_Paclitaxel, d)) and d.startswith('TCGA')])
##directories_doc = sorted([d for d in os.listdir(path_in_BRCA_Docetaxel) if os.path.isdir(os.path.join(path_in_BRCA_Docetaxel, d)) and d.startswith('TCGA')])
##directories_cyc = sorted([d for d in os.listdir(path_in_BRCA_Cyclophosphamide) if os.path.isdir(os.path.join(path_in_BRCA_Cyclophosphamide, d)) and d.startswith('TCGA')])


#(_, dox_patients_to_process , _ ) = venndiagramsets('fluorouracil', set(directories_fl), 'Doxorubicin', set( directories_dox) )
#(_, pac_patients_to_process , _ ) = venndiagramsets('fluorouracil + Doxorubicin', set(directories_fl) | set( directories_dox), 'Paclitaxel', set( directories_pac) )
#(_, doc_patients_to_process , _ ) = venndiagramsets('fluorouracil + Doxorubicin+ Paclitaxel', set(directories_fl) | set( directories_dox) |  set( directories_pac) , 'Docetaxel', set( directories_doc) )
#(_, cyc_patients_to_process , _ ) = venndiagramsets('fluorouracil + Doxorubicin+ Paclitaxel + Docetaxel', set(directories_fl) | set( directories_dox) |  set( directories_pac) | set( directories_doc) , 'Cyclophosphamide', set( directories_cyc) )

# setting global path values for BRCA data 

#path_in = path_in_BRCA_Doxorubici
#path_out = path_out_BRCA_Doxorubici

#path_in = path_in_BRCA_Paclitaxel
#path_out = path_out_BRCA_Paclitaxel

##path_in = path_in_BRCA_Docetaxel
##path_out = path_out_BRCA_Docetaxel
##
#path_in = path_in_BRCA_Cyclophosphamide
#path_out = path_out_BRCA_Cyclophosphamide

#  setting global path values for STAD data 
#path_in =  '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Stomach_cancer_fluorouracil'
#path_out = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Stomach_cancer_sp2_fluorouracil'


##########################

def generate_maps() :

    testexistencedirectory(path_in)
    testexistenceanscreatedirectory(path_out)

    directories = sorted([d for d in os.listdir(path_in) if os.path.isdir(os.path.join(path_in, d)) and d.startswith('TCGA')])

    printlist(directories, 'directories', 10)

    dir_to_files_map_up = {}
    dir_to_files_map_down = {}
    dir_to_files_map_up_logfc = {}
    dir_to_files_map_down_logfc = {}

    for dirname in directories :
        new_path = os.path.join(path_in, dirname)
        files = [f for f in os.listdir(new_path) if os.path.isfile(os.path.join(new_path, f)) and f.startswith('TCGA')]
     #   print(get_patient_name(dirname), files)
        dir_to_files_map_up[dirname] = [x for x in files if x.endswith('_UP_gene_list.txt')]
        dir_to_files_map_down[dirname] = [x for x in files if x.endswith('_DOWN_gene_list.txt')]
        dir_to_files_map_up_logfc[dirname] = [x for x in files if x.endswith('_UP_gene_logfc.txt')]
        dir_to_files_map_down_logfc[dirname] = [x for x in files if x.endswith('_DOWN_gene_logfc.txt')]

    printmap(dir_to_files_map_up, 'dir_to_files_map_up', 100)
    printmap(dir_to_files_map_down, 'dir_to_files_map_down', 100)
    printmap(dir_to_files_map_up_logfc, 'dir_to_files_map_up_logfc', 100)
    printmap(dir_to_files_map_down_logfc, 'dir_to_files_map_down_logfc', 100)

    return (dir_to_files_map_up, dir_to_files_map_down, dir_to_files_map_up_logfc, dir_to_files_map_down_logfc)


(dir_to_files_map_up, dir_to_files_map_down, dir_to_files_map_up_logfc, dir_to_files_map_down_logfc) =  generate_maps()

# 
    
CRC_irinotecan_patients = [
'TCGA-A6-2671-01A',
'TCGA-A6-2683-01A',
'TCGA-A6-5664-01A',
'TCGA-AA-3494-01A',
'TCGA-AA-3517-01A',
'TCGA-AA-3548-01A',
'TCGA-AA-3679-01A',
'TCGA-AA-3692-01A',
'TCGA-AA-3869-01A',
'TCGA-AA-3973-01A',
'TCGA-AA-A01D-01A',
'TCGA-AA-A01R-01A',
'TCGA-AA-A02F-01A',
'TCGA-AD-6965-01A',
'TCGA-AF-2689-01A',
'TCGA-AG-3584-01A',
'TCGA-AG-3731-01A',
'TCGA-AG-3887-01A',
'TCGA-AG-A02G-01A',
'TCGA-AY-4070-01A',
'TCGA-AZ-4682-01A',
'TCGA-AZ-6600-01A',
'TCGA-CM-4748-01A',
'TCGA-CM-5868-01A',
'TCGA-CM-6675-01A',
'TCGA-CM-6678-01A',
'TCGA-DC-6156-01A',
'TCGA-F4-6806-01A',
'TCGA-F5-6571-01A',
'TCGA-G4-6294-01A',
'TCGA-G4-6297-01A',
'TCGA-G4-6298-01A',
'TCGA-G4-6303-01A',
'TCGA-NH-A6GA-01A',
    ]


CRC_capecitabine_patients = [
'TCGA-4N-A93T-01',
'TCGA-A6-5662-01',
'TCGA-AA-3494-01',
'TCGA-AA-3844-01',
'TCGA-AA-3869-01',
'TCGA-AA-3955-01',
'TCGA-AA-3994-01',
'TCGA-AA-A01D-01',
'TCGA-AD-6889-01',
'TCGA-AD-6901-01',
'TCGA-AF-2687-01',
'TCGA-AF-2692-01',
'TCGA-AF-A56N-01',
'TCGA-AG-3732-01',
'TCGA-AG-3887-01',
'TCGA-AG-A01N-01',
'TCGA-AH-6544-01',
'TCGA-AY-4070-01',
'TCGA-AZ-4615-01',
'TCGA-AZ-4616-01',
'TCGA-CM-4743-01',
'TCGA-D5-6538-01',
'TCGA-DY-A1DC-01',
'TCGA-EI-7004-01',
'TCGA-F4-6460-01',
'TCGA-F4-6805-01',
'TCGA-F4-6806-01',
'TCGA-F4-6807-01',
'TCGA-F5-6571-01',
'TCGA-F5-6863-01',
'TCGA-G4-6297-01',
'TCGA-G4-6309-01',
'TCGA-NH-A50U-01',
'TCGA-NH-A8F7-01',
'TCGA-NH-A8F8-01',
'TCGA-QG-A5YV-01',
'TCGA-QG-A5YW-01',
'TCGA-QG-A5YX-01'
    ]


STAD_flu_patients_to_process = [
 'TCGA-3M-AB46-01A',
 'TCGA-BR-6457-01A',
 'TCGA-BR-6458-01A',
 'TCGA-BR-6706-01A',
 'TCGA-BR-6707-01A',
 'TCGA-BR-6710-01A',
 'TCGA-BR-6802-01A',
 'TCGA-BR-8373-01A',
 'TCGA-BR-8381-01A',
 'TCGA-BR-8590-01A',
 'TCGA-BR-8592-01A',
 'TCGA-BR-8690-01A',
 'TCGA-BR-A44T-01A',
 'TCGA-D7-6520-01A',
 'TCGA-F1-A448-01A',
 'TCGA-FP-7829-01A',
 'TCGA-FP-7998-01A',
 'TCGA-FP-A4BE-01A',
 'TCGA-HJ-7597-01A',
 'TCGA-MX-A666-01A',
 'TCGA-VQ-A8E0-01A',
 'TCGA-VQ-A8PD-01A',
 'TCGA-VQ-A8PP-01A',
 'TCGA-VQ-A91D-01A',
 'TCGA-VQ-A91E-01A',
 'TCGA-VQ-A91K-01A',
 'TCGA-VQ-A91N-01A',
 'TCGA-VQ-A91Q-01A',
 'TCGA-VQ-A91S-01A',
 'TCGA-VQ-A91V-01A',
 'TCGA-VQ-A91W-01A',
 'TCGA-VQ-A91Y-01A',
 'TCGA-VQ-A91Z-01A',
 'TCGA-VQ-A922-01A',
 'TCGA-VQ-A928-01A',
 'TCGA-VQ-A94O-01A',
 'TCGA-VQ-A94R-01A',
 'TCGA-VQ-AA64-01A',
 'TCGA-VQ-AA68-01A',
 'TCGA-VQ-AA6B-01A',
 'TCGA-VQ-AA6D-01A',
 'TCGA-VQ-AA6F-01A',
 'TCGA-VQ-AA6G-01A',
 'TCGA-VQ-AA6K-01A',
    ]


def handling_irinotecan_patients () :

    present_irinotecan_patients = set(irinotecan_patients) & set(directories)
    absent_irinotecan_patients = set(irinotecan_patients) - set(directories)
    #present_irinotecan_patients = (set(irinotecan_patients) & set(directories)) - set(['TCGA-AZ-6600-01A', 'TCGA-AA-3517-01A'])


    print('present data' , len(present_irinotecan_patients), 'absent data', len(absent_irinotecan_patients))
    printlist(list(present_irinotecan_patients), 'present_irinotecan_patients')
    printlist(list(absent_irinotecan_patients), 'absent_irinotecan_patients')


    present_capecitabine_patients = (set([get_patient_name(x) for x in capecitabine_patients]) & set([get_patient_name(x) for x in directories])) 

    print('present data' , len(present_capecitabine_patients))
    printlist(list(present_capecitabine_patients), 'present_capecitabine_patients')

    present_capecitabine_patients2 = [x + '-01A' for x in present_capecitabine_patients]

    printlist(list(present_capecitabine_patients2), 'present_capecitabine_patients2')


    # intersection of teh two lsits

    int_irintitecan_and_capecitabin = set(present_irinotecan_patients) & set(present_capecitabine_patients2)
    printlist(list(int_irintitecan_and_capecitabin), 'int_irintitecan_and_capecitabin')

    return present_irinotecan_patients
     

##patient_up_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_' + patient_id +  '_UP_gene_list.txt'
##patient_down_genes_file = '../../New_R_code_2025_04_15/DrugMerge_SP_Final/1-input_files/DEGs/Colorectal_cancer_sp/TCGA-COAD_' + patient_id +  '_DOWN_gene_list.txt'
##testexistencefile('.', patient_up_genes_file)
##testexistencefile('.', patient_down_genes_file)


def generate_new_up_down_lists(patient_id , sample_id, threshold_up, threshold_down, dump_to_files) :

    new_patient_up_genes_file = os.path.join(path_in, sample_id, dir_to_files_map_up[sample_id][0])
    new_patient_down_genes_file = os.path.join(path_in, sample_id, dir_to_files_map_down[sample_id][0])

    new_patient_up_genes_file_logfc = os.path.join(path_in, sample_id, dir_to_files_map_up_logfc[sample_id][0])
    new_patient_down_genes_file_logfc = os.path.join(path_in, sample_id, dir_to_files_map_down_logfc[sample_id][0])

    testexistencefile('.', new_patient_up_genes_file)
    testexistencefile('.', new_patient_down_genes_file)
    testexistencefile('.', new_patient_up_genes_file_logfc)
    testexistencefile('.', new_patient_down_genes_file_logfc)

    print('*******************************')
    acc_fc_down, _ = load_data('.', new_patient_down_genes_file_logfc, '\t', 'verbose')
    new_s_down = load_gene_list('.', new_patient_down_genes_file)
    
    printlist(acc_fc_down, 'acc_fc_down', 10)
    
    print('*******************************')
    acc_fc_up , _ = load_data('.', new_patient_up_genes_file_logfc, '\t', 'verbose')
    new_s_up = load_gene_list('.', new_patient_up_genes_file)

#    print(len(acc_fc_up), acc_fc_up[0], acc_fc_up[1])
    printlist(acc_fc_up, 'acc_fc_up', 10)

    print('*******************************')

# handing UP genes

    in_list_up = [x[0].split('.')[0] for x in acc_fc_up]
    printlist(in_list_up, 'in_list_up', 10)

    res_up = idmap.convert(in_list_up ,'ensg','symb')
    printlist(res_up, 'result in symb coding  for UP', 10)

    not_empty_symbols_up = [x for x in res_up if x != 'N/A']

    print('counting not N/A in up ', len(not_empty_symbols_up), 'set size ' ,  len(set(not_empty_symbols_up)))

    venndiagramsets('original symbols', new_s_up,  'converted symbols', set(not_empty_symbols_up))

    assert len(res_up) == len(in_list_up)

# filrering by FC 
    up_mapping = {}
    for i in range(len(res_up)) :
        symbol = res_up[i]
        fc_value = float(acc_fc_up[i][1])
        if symbol != 'N/A' :
            if symbol not in up_mapping.keys() and fc_value >= threshold_up :
                up_mapping[symbol] =  fc_value

    printmap(up_mapping, 'new up_mapping', 10)


# handing DOWN genes

    in_list_down = [x[0].split('.')[0] for x in acc_fc_down]
    printlist(in_list_down, 'in_list_down', 10)

    res_down = idmap.convert(in_list_down ,'ensg','symb')
    printlist(res_down, 'result in symb coding for DOWN', 10)
    not_empty_symbols_down = [x for x in res_down if x != 'N/A']

    print('counting not N/A in down ', len(not_empty_symbols_down), 'set size ' ,  len(set(not_empty_symbols_down)))

    venndiagramsets('original symbols', new_s_down,  'converted symbols', set(not_empty_symbols_down))

    assert len(res_down) == len(in_list_down)

# filrering by FC 
    down_mapping = {}
    for i in range(len(res_down)) :
        symbol = res_down[i]
        fc_value = float(acc_fc_down[i][1])
        if symbol != 'N/A' :
            if symbol not in down_mapping.keys() and fc_value <=  threshold_down :
                down_mapping[symbol] =  fc_value

    printmap(down_mapping, 'new down_mapping', 10)

    if dump_to_files == True :

        f_up = open(os.path.join(path_out, patient_id +  '_' + str(threshold_up)  + '_UP_gene_list.txt' ) , 'w')
        for x in up_mapping.keys() :
            f_up.write(x + '\n')
        f_up.close
        
        f_down = open(os.path.join(path_out, patient_id + '_' + str(threshold_down) + '_DOWN_gene_list.txt' ) , 'w')
        for x in down_mapping.keys() :
            f_down.write(x + '\n')
        f_down.close      


###################
#patient_id_1 = 'TCGA-AZ-6600'
#sample_id_1 = 'TCGA-AZ-6600-01A'

##patient_id_1 = 'TCGA-AA-3517'
##sample_id_1 = 'TCGA-AA-3517-01A'

          

#generate_new_up_down_lists(patient_id_1 , sample_id_1, 1.0, -1.0, False)
#generate_new_up_down_lists(patient_id_1 , sample_id_1, 1.5, -1.5, True)


##for x in present_irinotecan_patients :
##    sample_id_1 = x
##    patient_id_1 = get_patient_name(x)
##    generate_new_up_down_lists(patient_id_1 , sample_id_1, 1.0, -1.0, True)
    

##
##s_up = load_gene_list('.', patient_up_genes_file)
##new_s_up = load_gene_list('.', new_patient_up_genes_file)
##
##s_down = load_gene_list('.', patient_down_genes_file)

##
##venndiagramsets('s_up', s_up, 'new_s_up', new_s_up)
##venndiagramsets('s_down', s_down, 'new_s_down', new_s_down)


###########################
# a code to generate 'balanced' lists if genes up and down via the FC
# also call the conversion after the filtering.



def generate_new_up_down_lists_balanced(patient_id , sample_id, size_limit, dump_to_files) :

    new_patient_up_genes_file = os.path.join(path_in, sample_id, dir_to_files_map_up[sample_id][0])
    new_patient_down_genes_file = os.path.join(path_in, sample_id, dir_to_files_map_down[sample_id][0])

    new_patient_up_genes_file_logfc = os.path.join(path_in, sample_id, dir_to_files_map_up_logfc[sample_id][0])
    new_patient_down_genes_file_logfc = os.path.join(path_in, sample_id, dir_to_files_map_down_logfc[sample_id][0])

    testexistencefile('.', new_patient_up_genes_file)
    testexistencefile('.', new_patient_down_genes_file)
    testexistencefile('.', new_patient_up_genes_file_logfc)
    testexistencefile('.', new_patient_down_genes_file_logfc)

    print('*******************************')
    acc_fc_down, _ = load_data('.', new_patient_down_genes_file_logfc, '\t', 'verbose')
    acc_fc_up , _ = load_data('.', new_patient_up_genes_file_logfc, '\t', 'verbose')

    list_up = [(x[0].split('.')[0], float(x[1]))  for x in acc_fc_up]
    list_down = [(x[0].split('.')[0], float(x[1]))  for x in acc_fc_down]

    sorted_list_up = sorted(list_up, key=lambda tup: tup[1], reverse=True)
    printlist(sorted_list_up, 'sorted_list_up', 10)            
    sorted_list_down = sorted(list_down, key=lambda tup: tup[1], reverse=False)
    printlist(sorted_list_down, 'sorted_list_down', 10)

##    counter_up = 0
##    last_up_fc = None
##    out_up = []
##    for item in sorted_list_up :
##        print(item[0])
##        genename = idmap.convert(item[0] ,'ensg','symb')
##        print(genename)
##        if genename != 'N/A' :
##            out_up.append(genename)
##            counter_up +=1
##            last_up_fc = item[1]
##        if counter_up >= size_limit :
##            break
##    print('Up genes size  {}, last FC {}'.format(len(), last_up_fc))
    
    reduced_list_up = [x[0] for x in sorted_list_up[:size_limit] ]
    res_up = idmap.convert(reduced_list_up ,'ensg','symb')
    printlist(reduced_list_up, 'result in symb coding for UP', 10)
    not_empty_symbols_up = [x for x in res_up  if x != 'N/A']
    printlist(not_empty_symbols_up, 'Non N/A result in symb coding for UP', 10)
    last_up_index = min(size_limit -1, len(sorted_list_up) -1 )
    print('Last up FC value {}'.format(sorted_list_up[last_up_index]))
    
    reduced_list_down = [x[0] for x in sorted_list_down[:size_limit] ]
    res_down = idmap.convert(reduced_list_down ,'ensg','symb')
    printlist(reduced_list_down, 'result in symb coding for DOWN', 10)
    not_empty_symbols_down = [x for x in res_down  if x != 'N/A']
    printlist(not_empty_symbols_down, 'Non N/A result in symb coding for DOWN', 10)
    last_down_index = min(size_limit -1, len(sorted_list_down) -1 )
    print('Last down FC value {}'.format(sorted_list_down[last_down_index]))
      
    if dump_to_files == True :

        f_up = open(os.path.join(path_out, patient_id +  '_' + str(size_limit)  + '_UP_gene_list.txt' ) , 'w')
        for x in not_empty_symbols_up :
            f_up.write(x + '\n')
        f_up.close
        
        f_down = open(os.path.join(path_out, patient_id + '_' + str(size_limit) + '_DOWN_gene_list.txt' ) , 'w')
        for x in not_empty_symbols_down  :
            f_down.write(x + '\n')
        f_down.close      
            

#  exec
#generate_new_up_down_lists_balanced(patient_id_1 , sample_id_1, 4000, True)


global_size_limit = 4000



#present_irinotecan_patients = handling_irinotecan_patients ()
##
###for x in present_irinotecan_patients :
#for x in     set(['TCGA-AZ-6600-01A', 'TCGA-AA-3517-01A']) :
##for x in     ['TCGA-F5-6571-01A', 'TCGA-DC-6156-01A', 'TCGA-AG-A02G-01A', 'TCGA-AG-3887-01A', 'TCGA-AG-3731-01A', 'TCGA-AG-3584-01A']:
##for x in ['TCGA-AG-3584-01A'] :
##    sample_id_2 = x
##    patient_id_2 = get_patient_name(x)
##    generate_new_up_down_lists_balanced(patient_id_2 , sample_id_2, global_size_limit, True)
                

##sample_id_2 = 'TCGA-AZ-4682-01B'
##patient_id_2 = 'TCGA-AZ-4682'
##generate_new_up_down_lists_balanced(patient_id_2 , sample_id_2, global_size_limit, True)


##for x in ['TCGA-AG-3584-01A'] :
##    sample_id_2 = x
##    patient_id_2 = get_patient_name(x)
##    generate_new_up_down_lists_balanced(patient_id_2 , sample_id_2, global_size_limit, True)

def main_function_one(x) :
    
    sample_id_2 = x
    patient_id_2 = get_patient_name(x)
    generate_new_up_down_lists_balanced(patient_id_2 , sample_id_2, global_size_limit, True)
    

###############
    
def main_function_all(patients_to_process) : # all but the first one.   used for testing.. 
    
#    for x in directories[1:]: # all but the first one.
    for x in patients_to_process:
        sample_id_2 = x
        patient_id_2 = get_patient_name(x)
        generate_new_up_down_lists_balanced(patient_id_2 , sample_id_2, global_size_limit, True)
  

#################
# first BRCA fluoruraciil
##one_test = 'TCGA-BRCA_TCGA-A8-A0AB-01A'
##one_test ='TCGA-BRCA_TCGA-A8-A07C-01A'
##main_function_one(one_test)

#main_function_all(directories_fl )

#################
# first BRCA  Doxorubicin Hydrochloride
#one_test = 'TCGA-BRCA_TCGA-A2-A0CK-01A'
#main_function_one(one_test)

#main_function_all(dox_patients_to_process)

#################
# first BRCA  Paclitaxel
#main_function_all(pac_patients_to_process)


#################
# first BRCA  Docetaxel 
#main_function_all(doc_patients_to_process)


#################
# first BRCA  Cyclophosphamide
#main_function_all(cyc_patients_to_process)

##################################
# STAD patientes fluorouracil

#STAD_flu_patients_to_process2 = ['TCGA-STAD_' + x for x in STAD_flu_patients_to_process ]
#main_function_all(STAD_flu_patients_to_process2)


