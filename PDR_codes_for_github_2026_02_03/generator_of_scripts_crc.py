#  script to generate command lines for CRC patients data analysis

##
##
##
##docetaxel_pat_list = [ 'TCGA-BRCA_TCGA-A7-A4SD-01A',
##             'TCGA-BRCA_TCGA-AC-A3TM-01A',
##             'TCGA-BRCA_TCGA-A7-A4SF-01A',
##             'TCGA-BRCA_TCGA-OL-A66O-01A',
##             'TCGA-BRCA_TCGA-GM-A2DC-01A',
##             'TCGA-BRCA_TCGA-LL-A5YO-01A',
##             'TCGA-BRCA_TCGA-A7-A6VX-01A',
##             'TCGA-BRCA_TCGA-A2-A3XV-01A', 
##              'TCGA-BRCA_TCGA-EW-A424-01A',
##             'TCGA-BRCA_TCGA-BH-A0DD-01A',
##             'TCGA-BRCA_TCGA-A8-A08O-01A',
##             'TCGA-BRCA_TCGA-A7-A4SB-01A',
##             'TCGA-BRCA_TCGA-EW-A3U0-01A',
##             'TCGA-BRCA_TCGA-GM-A5PV-01A',
##             'TCGA-BRCA_TCGA-LL-A441-01A',
##             'TCGA-BRCA_TCGA-AC-A2FB-01A',
##             'TCGA-BRCA_TCGA-BH-A18F-01A',
##             'TCGA-BRCA_TCGA-BH-A0HY-01A',
##             'TCGA-BRCA_TCGA-A8-A097-01A',
##             'TCGA-BRCA_TCGA-BH-A18I-01A',
##             'TCGA-BRCA_TCGA-E2-A573-01A',
##             'TCGA-BRCA_TCGA-A7-A5ZW-01A',
##             'TCGA-BRCA_TCGA-A7-A6VV-01A',
##             'TCGA-BRCA_TCGA-BH-A8G0-01A',
##             'TCGA-BRCA_TCGA-LD-A7W5-01A',
##             'TCGA-BRCA_TCGA-XX-A899-01A',
##             'TCGA-BRCA_TCGA-S3-A6ZF-01A',
##             'TCGA-BRCA_TCGA-A2-A3XT-01A',
##             'TCGA-BRCA_TCGA-A2-A3XZ-01A',
##             'TCGA-BRCA_TCGA-OL-A66P-01A',
##             'TCGA-BRCA_TCGA-BH-A202-01A',
##             'TCGA-BRCA_TCGA-OL-A6VO-01A'
##]
##
###drug_name =  'docetaxel'
##
### parietns shareed with other drug lsits
##pat_list = {'TCGA-BRCA_TCGA-A2-A3XS-01A', 'TCGA-BRCA_TCGA-E9-A1R2-01A', 'TCGA-BRCA_TCGA-BH-A0DV-01A', 'TCGA-BRCA_TCGA-A2-A3Y0-01A', 'TCGA-BRCA_TCGA-A2-A0CK-01A', 'TCGA-BRCA_TCGA-A2-A3XY-01A', 'TCGA-BRCA_TCGA-GM-A2DI-01A', 'TCGA-BRCA_TCGA-A8-A07C-01A', 'TCGA-BRCA_TCGA-GM-A4E0-01A', 'TCGA-BRCA_TCGA-AC-A2FK-01A', 'TCGA-BRCA_TCGA-BH-A18H-01A', 'TCGA-BRCA_TCGA-AC-A3BB-01A', 'TCGA-BRCA_TCGA-A2-A4RY-01A', 'TCGA-BRCA_TCGA-BH-A0DG-01A', 'TCGA-BRCA_TCGA-GM-A2DA-01A', 'TCGA-BRCA_TCGA-BH-A0DT-01A'}
##
##########################################
##
##pactaxel_pat_list = ['TCGA-BRCA_TCGA-OL-A5S0-01A',
## 'TCGA-BRCA_TCGA-D8-A3Z6-01A',
## 'TCGA-BRCA_TCGA-A7-A6VY-01A',
## 'TCGA-BRCA_TCGA-OL-A5RV-01A',
## 'TCGA-BRCA_TCGA-C8-A1HM-01A',
## 'TCGA-BRCA_TCGA-BH-A0AZ-01A',
## 'TCGA-BRCA_TCGA-EW-A2FR-01A',
## 'TCGA-BRCA_TCGA-S3-A6ZG-01A',
## 'TCGA-BRCA_TCGA-A7-A6VW-01A',
## 'TCGA-BRCA_TCGA-EW-A1OY-01A',
## 'TCGA-BRCA_TCGA-A7-A4SA-01A']
##
###drug_name = 'paclitaxel'
##
###pat_list = {'TCGA-BRCA_TCGA-BH-A0B5-01A', 'TCGA-BRCA_TCGA-GM-A2DF-01A', 'TCGA-BRCA_TCGA-A2-A4S2-01A', 'TCGA-BRCA_TCGA-C8-A27B-01A', 'TCGA-BRCA_TCGA-AC-A2FF-01A', 'TCGA-BRCA_TCGA-GM-A3XG-01A', 'TCGA-BRCA_TCGA-A7-A4SE-01A', 'TCGA-BRCA_TCGA-BH-A0DI-01A', 'TCGA-BRCA_TCGA-A2-A3XY-01A', 'TCGA-BRCA_TCGA-A2-A0CK-01A', 'TCGA-BRCA_TCGA-GM-A2DM-01A', 'TCGA-BRCA_TCGA-GM-A3XN-01A', 'TCGA-BRCA_TCGA-GM-A2DH-01A', 'TCGA-BRCA_TCGA-GM-A2DN-01A', 'TCGA-BRCA_TCGA-LL-A740-01A', 'TCGA-BRCA_TCGA-OL-A66I-01A', 'TCGA-BRCA_TCGA-AC-A2BK-01A', 'TCGA-BRCA_TCGA-MS-A51U-01A', 'TCGA-BRCA_TCGA-GM-A2DB-01A', 'TCGA-BRCA_TCGA-C8-A27A-01A', 'TCGA-BRCA_TCGA-A2-A3XS-01A', 'TCGA-BRCA_TCGA-A2-A3XW-01A', 'TCGA-BRCA_TCGA-A2-A3XX-01A', 'TCGA-BRCA_TCGA-GM-A3XL-01A', 'TCGA-BRCA_TCGA-LL-A7SZ-01A', 'TCGA-BRCA_TCGA-OL-A5RW-01A', 'TCGA-BRCA_TCGA-A2-A0EP-01A', 'TCGA-BRCA_TCGA-E9-A6HE-01A', 'TCGA-BRCA_TCGA-A2-A0CR-01A'}
##
############################################
##doxorubicin_pat_list =[
##'TCGA-BRCA_TCGA-BH-A18H-01A', 
##'TCGA-BRCA_TCGA-BH-A0DI-01A', 
##'TCGA-BRCA_TCGA-AC-A2FK-01A', 
##'TCGA-BRCA_TCGA-LL-A740-01A', 
##'TCGA-BRCA_TCGA-A7-A4SE-01A', 
##'TCGA-BRCA_TCGA-A2-A3XS-01A', 
##'TCGA-BRCA_TCGA-E9-A6HE-01A', 
##'TCGA-BRCA_TCGA-BH-A0B5-01A',
##'TCGA-BRCA_TCGA-A2-A4RY-01A', 
##'TCGA-BRCA_TCGA-OL-A66I-01A', 
##'TCGA-BRCA_TCGA-LL-A73Z-01A', 
##'TCGA-BRCA_TCGA-A2-A3XX-01A', 
##'TCGA-BRCA_TCGA-MS-A51U-01A',
##'TCGA-BRCA_TCGA-A2-A3Y0-01A', 
##'TCGA-BRCA_TCGA-BH-A0BS-01A', 
##'TCGA-BRCA_TCGA-E9-A1R2-01A', 
##'TCGA-BRCA_TCGA-AC-A2FF-01A', 
##'TCGA-BRCA_TCGA-AC-A3BB-01A', 
##'TCGA-BRCA_TCGA-A2-A4S2-01A', 
##'TCGA-BRCA_TCGA-BH-A0DG-01A', 
##'TCGA-BRCA_TCGA-AC-A2FE-01A', 
##'TCGA-BRCA_TCGA-A2-A0EP-01A', 
##'TCGA-BRCA_TCGA-BH-A0C3-01A',
## 'TCGA-BRCA_TCGA-OL-A5RW-01A', 
## 'TCGA-BRCA_TCGA-A2-A3XY-01A', 
## 'TCGA-BRCA_TCGA-A2-A3XW-01A', 
## 'TCGA-BRCA_TCGA-BH-A0BT-01A', 
## 'TCGA-BRCA_TCGA-LL-A7SZ-01A', 
## 'TCGA-BRCA_TCGA-A2-A0CR-01A', 
## 'TCGA-BRCA_TCGA-BH-A0DV-01A', 
## 'TCGA-BRCA_TCGA-S3-A6ZH-01A', 
## 'TCGA-BRCA_TCGA-A2-A0CK-01A', 
## 'TCGA-BRCA_TCGA-BH-A0DT-01A', 
## 'TCGA-BRCA_TCGA-AC-A2BK-01A'
##    ]
##
###drug_name = 'doxorubicin'
##
###################################
# stad fluorouracil

##pat_list = [
### 'TCGA-3M-AB46-01A',
## 'TCGA-BR-6457-01A',
## 'TCGA-BR-6458-01A',
## 'TCGA-BR-6706-01A',
## 'TCGA-BR-6707-01A',
## 'TCGA-BR-6710-01A',
## 'TCGA-BR-6802-01A',
## 'TCGA-BR-8373-01A',
## 'TCGA-BR-8381-01A',
## 'TCGA-BR-8590-01A',
## 'TCGA-BR-8592-01A',
## 'TCGA-BR-8690-01A',
## 'TCGA-BR-A44T-01A',
## 'TCGA-D7-6520-01A',
## 'TCGA-F1-A448-01A',
## 'TCGA-FP-7829-01A',
## 'TCGA-FP-7998-01A',
## 'TCGA-FP-A4BE-01A',
## 'TCGA-HJ-7597-01A',
## 'TCGA-MX-A666-01A',
#### 'TCGA-VQ-A8E0-01A',
## 'TCGA-VQ-A8PD-01A',
## 'TCGA-VQ-A8PP-01A',
## 'TCGA-VQ-A91D-01A',
## 'TCGA-VQ-A91E-01A',
## 'TCGA-VQ-A91K-01A',
## 'TCGA-VQ-A91N-01A',
## 'TCGA-VQ-A91Q-01A',
## 'TCGA-VQ-A91S-01A',
## 'TCGA-VQ-A91V-01A',
## 'TCGA-VQ-A91W-01A',
## 'TCGA-VQ-A91Y-01A',
## 'TCGA-VQ-A91Z-01A',
## 'TCGA-VQ-A922-01A',
## 'TCGA-VQ-A928-01A',
## 'TCGA-VQ-A94O-01A',
## 'TCGA-VQ-A94R-01A',
## 'TCGA-VQ-AA64-01A',
## 'TCGA-VQ-AA68-01A',
## 'TCGA-VQ-AA6B-01A',
## 'TCGA-VQ-AA6D-01A',
## 'TCGA-VQ-AA6F-01A',
## 'TCGA-VQ-AA6G-01A',
## 'TCGA-VQ-AA6K-01A',
##    ]


###################################
# stad fluorouracil_only

##pat_list = [
##'TCGA-BR-6710',
##'TCGA-BR-6457',
##'TCGA-3M-AB46',
##'TCGA-BR-A44T',
##'TCGA-FP-7829',
##'TCGA-FP-A4BE',
##'TCGA-MX-A666',
##'TCGA-VQ-A91E',
##'TCGA-VQ-A91K',
##'TCGA-VQ-A91Q',
##'TCGA-VQ-A91S',
##'TCGA-VQ-A91W',
##'TCGA-VQ-AA64',
##'TCGA-VQ-AA68',
##'TCGA-VQ-AA6B',
##'TCGA-VQ-AA6D',
##'TCGA-VQ-A8E0',
##'TCGA-VQ-A8PD',
##'TCGA-VQ-A8PP',
##'TCGA-VQ-A91D',
##'TCGA-VQ-A91N',
##'TCGA-VQ-A91Y',
##'TCGA-VQ-A922',
##'TCGA-VQ-A928',
##'TCGA-VQ-A94O',
##'TCGA-VQ-A94R',
##'TCGA-VQ-AA6F',
##'TCGA-VQ-AA6G',
##'TCGA-VQ-AA6K',
##]


###########################################
# CRC irinotecan 


##pat_list = [
###CR/PR
##'TCGA-AA-A02F',
##'TCGA-AA-A01R',
##'TCGA-AA-3869',
##'TCGA-F5-6571',
##'TCGA-AA-3517',
##'TCGA-AA-3548',
##'TCGA-AG-3731',
##'TCGA-F4-6806',
##'TCGA-AG-A02G',
###PD
##'TCGA-AZ-6600',
##'TCGA-A6-2671',
##'TCGA-AA-A01D',
##'TCGA-NH-A6GA',
##'TCGA-A6-5664',
##'TCGA-AG-3584',
###SD
##'TCGA-AA-3973',
##'TCGA-CM-4748',
##    ]


#################################################
#  CRC  oxaliplatin

pat_list = [
#PD/SD
'TCGA-NH-A8F8',
'TCGA-AF-6136',
'TCGA-A6-A566',
'TCGA-AH-6544',
'TCGA-NH-A50U',
'TCGA-NH-A6GB',
'TCGA-F5-6702',
#CR
'TCGA-AH-6903',
'TCGA-A6-A5ZU',
'TCGA-A6-A56B',
'TCGA-SS-A7HO',
'TCGA-AH-6643',
'TCGA-AY-A8YK',
'TCGA-AF-A56L',
'TCGA-NH-A50V',
]



    




#pat_list_2 = [x.split('_')[1]  for x in pat_list  ]
pat_list_2 = [x.split('_')[0]  for x in pat_list  ]
pat_list_3 = [x.split('-')[:3]  for x in pat_list_2]
pat_list_4 = ['-'.join(x)  for x in pat_list_3]
print(pat_list_4)


def gen_enrcichR_scripts(in_list, database) :
    for x in in_list : 

        print('./cl_My_script_colorectal_pipeline_2_functional.R  {} up 10 > log_{}_up_functional_{}.txt &'.format(x,x, database))
        print('sleep 3')
        print('./cl_My_script_colorectal_pipeline_2_functional.R  {} down 10 > log_{}_down_functional_{}.txt &'.format(x,x, database))
        print('\n')
        print('sleep 3')
        print('\n')
  
   
#run
#gen_enrcichR_scripts(pat_list_4, 'KEGG_2021_Human')


#drug_name = 'flurouracil'

########################

def gen_pdrugmerge_scripts_old(in_list) :

    date = '2025_09_22'
    pvals = [0.05, 0.01, 0.001]

    for x in in_list :
        print('\n\n-------------------------- Processing {} ----------------------'.format(x))
        for p in pvals :
            print('\n-------------------------\n\n')
            print('python3 cl_new_loader_routines_modulation.py {}  GEO_stad_{}  {}   >  log_new_pipeline_v3_{}_{}_all_p{}_FC4000_GEO_stad_{}.txt &'.format(x, drug_name, p,x,date,p, drug_name))
            print('\n')
            print('grep @@ log_new_pipeline_v3_{}_{}_all_p{}_FC4000_GEO_stad_{}.txt'.format(x, date, p, drug_name))


######################


def gen_pdrugmerge_scripts(in_list) :

    date = '2025_09_22'
    pvals = [0.05, 0.01, 0.001]

    for x in in_list :
        print('\n\n-------------------------- Processing {} ----------------------'.format(x))
        for p in pvals :
 #           print('\n-------------------------\n\n')
            print('python3 cl_new_loader_routines_modulation.py {}  GEO_stad_{}  {}   >  log_new_pipeline_v3_{}_{}_all_p{}_FC4000_GEO_stad_{}.txt &'.format(x, drug_name, p,x,date,p, drug_name))
 
    for x in in_list :
        print('\n\n-------------------------- Processing {} ----------------------'.format(x))
        for p in pvals :
             print('\n-------------------------\n\n')
             print('grep @@ log_new_pipeline_v3_{}_{}_all_p{}_FC4000_GEO_stad_{}.txt'.format(x, date, p, drug_name))



######################
#gen_pdrugmerge_scripts_(pat_list_4)


######################
# genration of drug functional commands


def gen_pdrugmerge_scripts_functional(in_list) :

    date = '2025_12_04'
#    database = 'Reactome_Pathways_2024'
#    database = 'GO_Biological_Process_2025'
    database = 'KEGG_2021_Human'

    ref_drug = 'fluorouracil'
#    ref_drug = 'irinotecan'

    for x in in_list :
        print('\n\n-------------------------- Processing {} ----------------------'.format(x))
        print('python3 cl_drug_functional_analysis.py   CRC {}   {}  {}  0.0001  GEO  0.01  >  log_{}_drugs_functional_rankings_{}_{}-1.txt  &'.format(ref_drug, x, database, x, database,date))
        print('grep "@@@@@@ For drug fluorouracil"  log_{}_drugs_functional_rankings_{}_{}-1.txt'.format(x, database, date))

 
##    for x in in_list :
##        print('\n\n-------------------------- Processing {} ----------------------'.format(x))
       


########
# run
#gen_pdrugmerge_scripts_funtional(pat_list_4)

