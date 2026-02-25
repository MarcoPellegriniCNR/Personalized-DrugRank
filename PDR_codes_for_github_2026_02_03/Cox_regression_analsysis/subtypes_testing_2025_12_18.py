# code to read the subtypes for several types of cancer

from new_loader_routines import load_data

from generic_util import venndiagramsets, printmap, printlist, increment_lists, domain_restriction


sub_STAD_fluorouracil_only = [
'TCGA-BR-6710',
'TCGA-BR-6457',
'TCGA-3M-AB46',
'TCGA-BR-A44T',
'TCGA-FP-7829',
'TCGA-FP-A4BE',
'TCGA-MX-A666',
'TCGA-VQ-A8E0',
'TCGA-VQ-A8PD',
'TCGA-VQ-A8PP',
'TCGA-VQ-A91D',
'TCGA-VQ-A91E',
'TCGA-VQ-A91K',
'TCGA-VQ-A91N',
'TCGA-VQ-A91Q',
'TCGA-VQ-A91S',
'TCGA-VQ-A91W',
'TCGA-VQ-A91Y',
'TCGA-VQ-A922',
'TCGA-VQ-A928',
'TCGA-VQ-A94O',
'TCGA-VQ-A94R',
'TCGA-VQ-AA64',
'TCGA-VQ-AA68',
'TCGA-VQ-AA6B',
'TCGA-VQ-AA6D',
'TCGA-VQ-AA6F',
'TCGA-VQ-AA6G',
'TCGA-VQ-AA6K',
    ]

CR_Fluoruracil_only_STAD_pats = [
'TCGA-BR-6710',
'TCGA-BR-6457',
'TCGA-3M-AB46',
'TCGA-BR-A44T',
'TCGA-FP-7829',
'TCGA-FP-A4BE',
'TCGA-MX-A666',
'TCGA-VQ-A91E',
'TCGA-VQ-A91K',
'TCGA-VQ-A91Q',
'TCGA-VQ-A91S',
'TCGA-VQ-A91W',
'TCGA-VQ-AA64',
'TCGA-VQ-AA68',
'TCGA-VQ-AA6B',
'TCGA-VQ-AA6D',
]

PD_Fluoruracil_only_STAD_pats = [
'TCGA-VQ-A8E0',
'TCGA-VQ-A8PD',
'TCGA-VQ-A8PP',
'TCGA-VQ-A91D',
'TCGA-VQ-A91N',
'TCGA-VQ-A91Y',
'TCGA-VQ-A922',
'TCGA-VQ-A928',
'TCGA-VQ-A94O',
'TCGA-VQ-A94R',
'TCGA-VQ-AA6F',
'TCGA-VQ-AA6G',
'TCGA-VQ-AA6K',
]





def STAD_analyze_subtypes(pat_id_list) :

    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'

    CLINICAL_DATA_FILE_NAME = 'STAD_Subtype.full.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, '\t', 'verbose')   
    acc = res[0]
    header = res[1]

    patient_idx =  header.index('Patient ID')
    subtype_idx = header.index('Subtype')
    subtype_map = {}

    pat_set = set()
    for x in acc :
        pat_id = x[patient_idx]
        if pat_id == '' :
            continue
        subtype = x[subtype_idx]
        pat_set.add(pat_id)
        subtype_map[pat_id] = subtype
        
    venndiagramsets('File' , pat_set, 'fluoruracil_only list', set(pat_id_list))       

    res_map = domain_restriction(pat_id_list, subtype_map, policy='any')
    printmap(res_map, 'Restricted STAD subtypes map')
    return res_map

    
##########





###########
def analyze_maps_STAD(res_map, CR_list, PD_list) :

    tp=0
    tn=0
    fp = 0
    fn = 0

    for key in res_map.keys():

        if res_map[key] in ['STAD_CIN', 'STAD_EBV']  and key in CR_list:
            tn += 1
        elif res_map[key] in ['STAD_CIN', 'STAD_EBV']  and key in PD_list:
            fn += 1
        elif res_map[key] in ['STAD_GS', 'STAD_MSI']  and key in CR_list:
            fp  +=1           
        elif res_map[key] in ['STAD_GS', 'STAD_MSI']  and key in PD_list:
            tp += 1
        else :
            print('Key {} not used'.format(key))

    neg = tn +fp
    pos = tp + fn

    acc = float(tp+tn)/float(neg+pos) if neg+pos > 0 else 0.0
    ppv = float(tp)/float(tp + fp) if  tp + fp >0  else  0.0 
    npv = float(tn)/float(tn + fn) if tn + fn > 0 else 0.0
    print('tn, fn, fp, tp, sum', tn, fn, fp, tp, tn+fn+fp+tp)
    print('Acc {}, ppv {} NPV {}'.format(acc, ppv, npv))

res_map = STAD_analyze_subtypes(sub_STAD_fluorouracil_only)
#analyze_maps_STAD(res_map, CR_Fluoruracil_only_STAD_pats, PD_Fluoruracil_only_STAD_pats)


#########################
#    BRCA lists

sub_BRCA_paracitaxel_4 = [
'TCGA-OL-A5S0',
'TCGA-D8-A3Z6',
'TCGA-A7-A6VY',
'TCGA-OL-A5RV',
'TCGA-C8-A1HM',
'TCGA-BH-A0AZ',
'TCGA-EW-A2FR',
'TCGA-S3-A6ZG',
'TCGA-A7-A6VW',
'TCGA-EW-A1OY',
'TCGA-A7-A4SA',]


shared_GEO_brca_paclitaxel = ['TCGA-C8-A27A', 'TCGA-LL-A740', 'TCGA-A2-A3XS', 'TCGA-C8-A27B', 'TCGA-GM-A2DF', 'TCGA-GM-A3XL', 'TCGA-AC-A2BK', 'TCGA-GM-A2DM',
                                                      'TCGA-A7-A4SE', 'TCGA-AC-A2FF', 'TCGA-GM-A3XG', 'TCGA-OL-A5RW', 'TCGA-A2-A3XW', 'TCGA-GM-A3XN', 'TCGA-OL-A66I', 'TCGA-BH-A0DI',
                                                      'TCGA-GM-A2DH', 'TCGA-A2-A0CK', 'TCGA-A2-A3XY', 'TCGA-BH-A0B5', 'TCGA-GM-A2DB', 'TCGA-A2-A0EP', 'TCGA-E9-A6HE', 'TCGA-A2-A4S2',
                                                      'TCGA-MS-A51U', 'TCGA-LL-A7SZ', 'TCGA-A2-A0CR', 'TCGA-A2-A3XX', 'TCGA-GM-A2DN']



sub_class_2_BRCA = [
'TCGA-GM-A2DC',
'TCGA-EW-A424',
'TCGA-S3-A6ZF',
'TCGA-LL-A5YO',
'TCGA-OL-A6VO',
'TCGA-OL-A66P',
'TCGA-OL-A66O',
'TCGA-GM-A5PV',
'TCGA-LL-A441',
'TCGA-EW-A3U0',
'TCGA-AC-A3TM',
'TCGA-A7-A4SB',
'TCGA-A7-A4SD',
'TCGA-A7-A5ZW',
'TCGA-A7-A6VV',
'TCGA-A7-A6VX',
'TCGA-A8-A097',
'TCGA-BH-A18F',
'TCGA-BH-A18I',
'TCGA-BH-A8G0',
'TCGA-E2-A573',
'TCGA-XX-A899',
]


sub_BRCA_class_3_docetaxel  = sub_class_2_BRCA + [  'TCGA-BH-A18H', 'TCGA-A2-A0CK', 'TCGA-BH-A0DG', 'TCGA-BH-A0DV', 'TCGA-A2-A3XY', 'TCGA-AC-A2FK', 'TCGA-A2-A3Y0', 'TCGA-A2-A4RY', 'TCGA-A2-A3XT', 'TCGA-A2-A3XS', 'TCGA-AC-A3BB', 'TCGA-BH-A0DT', 'TCGA-E9-A1R2']


sub_BRCA_class_0_fluorouracil  = [
'TCGA-GM-A3NW',
'TCGA-C8-A26X',
'TCGA-E9-A3X8',
'TCGA-A8-A07G',
'TCGA-A8-A0AB',
'TCGA-AC-A2QI',
'TCGA-AC-A3YI',
'TCGA-E9-A3Q9',
'TCGA-C8-A8HR',
'TCGA-E9-A1NE',
'TCGA-E9-A1NI',
'TCGA-E9-A2JS',
]


sub_BRCA_class__Doxorubicin_Hydrochloride =[
'TCGA-E9-A6HE',
'TCGA-BH-A0B5',
'TCGA-BH-A0DI',
'TCGA-AC-A2FF',
'TCGA-A2-A0CR',
'TCGA-LL-A740',
'TCGA-A2-A4S2',
'TCGA-A2-A3XY',
'TCGA-A7-A4SE',
'TCGA-A2-A3XX',
'TCGA-A2-A3XW',
'TCGA-OL-A66I',
'TCGA-OL-A5RW',
'TCGA-A2-A3XS',
'TCGA-MS-A51U',
'TCGA-LL-A7SZ',
'TCGA-A2-A0EP',
'TCGA-AC-A2BK',
]

#from clinical_data_analysis2 import map_BRCA_doxorubicin_geo_results

#sub_BRCA_class_X_doxorubicin = list(map_BRCA_doxorubicin_geo_results.keys())


###########################


def BRCA_analyze_subtypes(pat_id_list, in_string) :

    CLINICAL_DATA_PATH_NAME = '../../Clinical_Data_Directory'

    CLINICAL_DATA_FILE_NAME = 'BRCA_Subtype.full.csv'

    res = load_data(CLINICAL_DATA_PATH_NAME, CLINICAL_DATA_FILE_NAME, '\t', 'verbose')   
    acc = res[0]
    header = res[1]

    patient_idx =  header.index('Patient ID')
    subtype_idx = header.index('Subtype')
    subtype_map = {}

    pat_set = set()
    for x in acc :
        pat_id = x[patient_idx]
        if pat_id == '' :
            continue
        subtype = x[subtype_idx]
        pat_set.add(pat_id)
        subtype_map[pat_id] = subtype
        
    venndiagramsets('File' , pat_set, in_string , set(pat_id_list))       

    res_map = domain_restriction(pat_id_list, subtype_map, policy='any')
    printmap(res_map, 'Restricted BRCA  subtypes map')
    return res_map

    
#######################

res_map_sub_BRCA_paracitaxel_4  = BRCA_analyze_subtypes(sub_BRCA_paracitaxel_4 + shared_GEO_brca_paclitaxel , 'sub_BRCA_paracitaxel_4')
res_map_sub_BRCA_docetaxel_3 = BRCA_analyze_subtypes(sub_BRCA_class_3_docetaxel, 'sub_BRCA_docetaxel_3')
res_map_sub_BRCA_fluorouracil_0  = BRCA_analyze_subtypes(sub_BRCA_class_0_fluorouracil , 'sub_BRCA_class_0_fluorouracil')
res_map_sub_BRCA_Doxorubicin_Hydrochloride_E_1  = BRCA_analyze_subtypes(sub_BRCA_class__Doxorubicin_Hydrochloride , 'res_map_sub_BRCA_Doxorubicin_Hydrochloride_E_1')
res_map_sub_BRCA_docetaxel_D_2 = BRCA_analyze_subtypes(sub_class_2_BRCA, 'sub_BRCA_class_2_D_docetaxel')
