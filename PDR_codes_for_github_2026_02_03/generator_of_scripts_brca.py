# generator of comands for analysizing BRCA patients


docetaxel_pat_list = [ 'TCGA-BRCA_TCGA-A7-A4SD-01A',
             'TCGA-BRCA_TCGA-AC-A3TM-01A',
             'TCGA-BRCA_TCGA-A7-A4SF-01A',
             'TCGA-BRCA_TCGA-OL-A66O-01A',
             'TCGA-BRCA_TCGA-GM-A2DC-01A',
             'TCGA-BRCA_TCGA-LL-A5YO-01A',
             'TCGA-BRCA_TCGA-A7-A6VX-01A',
             'TCGA-BRCA_TCGA-A2-A3XV-01A', 
              'TCGA-BRCA_TCGA-EW-A424-01A',
             'TCGA-BRCA_TCGA-BH-A0DD-01A',
             'TCGA-BRCA_TCGA-A8-A08O-01A',
             'TCGA-BRCA_TCGA-A7-A4SB-01A',
             'TCGA-BRCA_TCGA-EW-A3U0-01A',
             'TCGA-BRCA_TCGA-GM-A5PV-01A',
             'TCGA-BRCA_TCGA-LL-A441-01A',
             'TCGA-BRCA_TCGA-AC-A2FB-01A',
             'TCGA-BRCA_TCGA-BH-A18F-01A',
             'TCGA-BRCA_TCGA-BH-A0HY-01A',
             'TCGA-BRCA_TCGA-A8-A097-01A',
             'TCGA-BRCA_TCGA-BH-A18I-01A',
             'TCGA-BRCA_TCGA-E2-A573-01A',
             'TCGA-BRCA_TCGA-A7-A5ZW-01A',
             'TCGA-BRCA_TCGA-A7-A6VV-01A',
             'TCGA-BRCA_TCGA-BH-A8G0-01A',
             'TCGA-BRCA_TCGA-LD-A7W5-01A',
             'TCGA-BRCA_TCGA-XX-A899-01A',
             'TCGA-BRCA_TCGA-S3-A6ZF-01A',
             'TCGA-BRCA_TCGA-A2-A3XT-01A',
             'TCGA-BRCA_TCGA-A2-A3XZ-01A',
             'TCGA-BRCA_TCGA-OL-A66P-01A',
             'TCGA-BRCA_TCGA-BH-A202-01A',
             'TCGA-BRCA_TCGA-OL-A6VO-01A'
]

#drug_name =  'docetaxel'

# parietns shareed with other drug lists
pat_list = {'TCGA-BRCA_TCGA-A2-A3XS-01A',
            'TCGA-BRCA_TCGA-E9-A1R2-01A',
            'TCGA-BRCA_TCGA-BH-A0DV-01A',
            'TCGA-BRCA_TCGA-A2-A3Y0-01A',
            'TCGA-BRCA_TCGA-A2-A0CK-01A',
            'TCGA-BRCA_TCGA-A2-A3XY-01A',
            'TCGA-BRCA_TCGA-GM-A2DI-01A',
            'TCGA-BRCA_TCGA-A8-A07C-01A',
            'TCGA-BRCA_TCGA-GM-A4E0-01A',
            'TCGA-BRCA_TCGA-AC-A2FK-01A',
            'TCGA-BRCA_TCGA-BH-A18H-01A',
            'TCGA-BRCA_TCGA-AC-A3BB-01A',
            'TCGA-BRCA_TCGA-A2-A4RY-01A',
            'TCGA-BRCA_TCGA-BH-A0DG-01A',
            'TCGA-BRCA_TCGA-GM-A2DA-01A',
            'TCGA-BRCA_TCGA-BH-A0DT-01A'}

########################################

pactaxel_pat_list = ['TCGA-BRCA_TCGA-OL-A5S0-01A',
 'TCGA-BRCA_TCGA-D8-A3Z6-01A',
 'TCGA-BRCA_TCGA-A7-A6VY-01A',
 'TCGA-BRCA_TCGA-OL-A5RV-01A',
 'TCGA-BRCA_TCGA-C8-A1HM-01A',
 'TCGA-BRCA_TCGA-BH-A0AZ-01A',
 'TCGA-BRCA_TCGA-EW-A2FR-01A',
 'TCGA-BRCA_TCGA-S3-A6ZG-01A',
 'TCGA-BRCA_TCGA-A7-A6VW-01A',
 'TCGA-BRCA_TCGA-EW-A1OY-01A',
 'TCGA-BRCA_TCGA-A7-A4SA-01A']

#drug_name = 'paclitaxel'

#pat_list = {'TCGA-BRCA_TCGA-BH-A0B5-01A', 'TCGA-BRCA_TCGA-GM-A2DF-01A', 'TCGA-BRCA_TCGA-A2-A4S2-01A', 'TCGA-BRCA_TCGA-C8-A27B-01A', 'TCGA-BRCA_TCGA-AC-A2FF-01A', 'TCGA-BRCA_TCGA-GM-A3XG-01A', 'TCGA-BRCA_TCGA-A7-A4SE-01A', 'TCGA-BRCA_TCGA-BH-A0DI-01A', 'TCGA-BRCA_TCGA-A2-A3XY-01A', 'TCGA-BRCA_TCGA-A2-A0CK-01A', 'TCGA-BRCA_TCGA-GM-A2DM-01A', 'TCGA-BRCA_TCGA-GM-A3XN-01A', 'TCGA-BRCA_TCGA-GM-A2DH-01A', 'TCGA-BRCA_TCGA-GM-A2DN-01A', 'TCGA-BRCA_TCGA-LL-A740-01A', 'TCGA-BRCA_TCGA-OL-A66I-01A', 'TCGA-BRCA_TCGA-AC-A2BK-01A', 'TCGA-BRCA_TCGA-MS-A51U-01A', 'TCGA-BRCA_TCGA-GM-A2DB-01A', 'TCGA-BRCA_TCGA-C8-A27A-01A', 'TCGA-BRCA_TCGA-A2-A3XS-01A', 'TCGA-BRCA_TCGA-A2-A3XW-01A', 'TCGA-BRCA_TCGA-A2-A3XX-01A', 'TCGA-BRCA_TCGA-GM-A3XL-01A', 'TCGA-BRCA_TCGA-LL-A7SZ-01A', 'TCGA-BRCA_TCGA-OL-A5RW-01A', 'TCGA-BRCA_TCGA-A2-A0EP-01A', 'TCGA-BRCA_TCGA-E9-A6HE-01A', 'TCGA-BRCA_TCGA-A2-A0CR-01A'}

##########################################
doxorubicin_pat_list =[
'TCGA-BRCA_TCGA-BH-A18H-01A', 
'TCGA-BRCA_TCGA-BH-A0DI-01A', 
'TCGA-BRCA_TCGA-AC-A2FK-01A', 
'TCGA-BRCA_TCGA-LL-A740-01A', 
'TCGA-BRCA_TCGA-A7-A4SE-01A', 
'TCGA-BRCA_TCGA-A2-A3XS-01A', 
'TCGA-BRCA_TCGA-E9-A6HE-01A', 
'TCGA-BRCA_TCGA-BH-A0B5-01A',
'TCGA-BRCA_TCGA-A2-A4RY-01A', 
'TCGA-BRCA_TCGA-OL-A66I-01A', 
'TCGA-BRCA_TCGA-LL-A73Z-01A', 
'TCGA-BRCA_TCGA-A2-A3XX-01A', 
'TCGA-BRCA_TCGA-MS-A51U-01A',
'TCGA-BRCA_TCGA-A2-A3Y0-01A', 
'TCGA-BRCA_TCGA-BH-A0BS-01A', 
'TCGA-BRCA_TCGA-E9-A1R2-01A', 
'TCGA-BRCA_TCGA-AC-A2FF-01A', 
'TCGA-BRCA_TCGA-AC-A3BB-01A', 
'TCGA-BRCA_TCGA-A2-A4S2-01A', 
'TCGA-BRCA_TCGA-BH-A0DG-01A', 
'TCGA-BRCA_TCGA-AC-A2FE-01A', 
'TCGA-BRCA_TCGA-A2-A0EP-01A', 
'TCGA-BRCA_TCGA-BH-A0C3-01A',
 'TCGA-BRCA_TCGA-OL-A5RW-01A', 
 'TCGA-BRCA_TCGA-A2-A3XY-01A', 
 'TCGA-BRCA_TCGA-A2-A3XW-01A', 
 'TCGA-BRCA_TCGA-BH-A0BT-01A', 
 'TCGA-BRCA_TCGA-LL-A7SZ-01A', 
 'TCGA-BRCA_TCGA-A2-A0CR-01A', 
 'TCGA-BRCA_TCGA-BH-A0DV-01A', 
 'TCGA-BRCA_TCGA-S3-A6ZH-01A', 
 'TCGA-BRCA_TCGA-A2-A0CK-01A', 
 'TCGA-BRCA_TCGA-BH-A0DT-01A', 
 'TCGA-BRCA_TCGA-AC-A2BK-01A'
    ]

#drug_name = 'doxorubicin'

#################################

pat_list_2 = [x.split('_')[1]  for x in pat_list  ]

pat_list_3 = [x.split('-')[:3]  for x in pat_list_2]

pat_list_4 = ['-'.join(x)  for x in pat_list_3]


print(pat_list_4)


def gen_enrichR_scripts(in_list) :
    for x in in_list : 

        print('./cl_My_script_breast_pipeline_2.R  {} up 3 > log_{}_up.txt &'.format(x,x))
        print('./cl_My_script_breast_pipeline_2.R  {} down 3 > log_{}_down.txt &'.format(x,x))
  
   
#run
#gen_enrichR_scripts(pat_list_4)

def gen_pdrugmerge_scripts(in_list) :

    date = '2025_08_03'
    pvals = [0.05, 0.01, 0.001]

    for x in in_list :
        print('\n\n-------------------------- Processing {} ----------------------'.format(x))
        for p in pvals :
            print('\n-------------------------\n\n')
            print('python3 cl_new_loader_routines_modulation.py {}  GEO_brca_{}  {}   >  log_new_pipeline_v3_{}_{}_all_p{}_FC4000_GEO_brca_{}.txt &'.format(x, drug_name, p,x,date,p, drug_name))
            print('\n')
            print('grep @@ log_new_pipeline_v3_{}_{}_all_p{}_FC4000_GEO_brca_{}.txt'.format(x, date, p, drug_name))


# run
#gen_pdrugmerge_scripts(pat_list_4)

