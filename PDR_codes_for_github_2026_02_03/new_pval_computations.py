#### code to perform pvaleu computations jsut on the basis of teh lenth of the data and hits lists

import scipy.stats

parameters_map = {}

def new_compute_zscores_pval(size_data, size_hits, inprecision, inrhr) :

    assert size_data >0
    assert size_hits > 0
    assert size_data >= size_hits

    key = (size_data, size_hits)
    if key not  in parameters_map.keys() :        
        res = new_generate_distibution_parameters(size_data, size_hits)
        parameters_map[key] = res

    precmean = parameters_map[key][0][0]
    precstd = parameters_map[key][0][1]
    rhrmean= parameters_map[key][1][0]
    rhrstd = parameters_map[key][1][1]

    z_prec = (inprecision - precmean)/precstd
    p_prec = scipy.stats.norm.sf(abs(z_prec))*2
    z_rhr = (inrhr - rhrmean)/ rhrstd
    p_rhr = scipy.stats.norm.sf(abs(z_rhr))*2

    return ((z_prec,p_prec),(z_rhr,p_rhr))

    
######################
# generate the mans and std for the chosen size parametes

from randomizedmodel import generate_random

def new_generate_distibution_parameters(size_data, size_hits):

    print ('Generating random distributions for pair : ' ,  size_data, size_hits)

#    outmap = {}

#    for x in datamap.keys() :
#    record = datamap[x]
#    numitems = record[0]
#    numhits = record[6]
    numitems = size_data
    numhits =  size_hits
    repeatflag = True
    iteration = 0
    while repeatflag :    
        (precat20_parameters1, rhr_parameters1) = generate_random(numitems, numhits, 10000)
        (precat20_parameters2, rhr_parameters2) = generate_random(numitems, numhits, 100000)
        err1 = abs(precat20_parameters1[0] - precat20_parameters2[0])
        err2 = abs(rhr_parameters1[0] - rhr_parameters2[0])        
        if  err1 > 0.01 or err2 > 0.01 :
            iteration += 1
            (precat20_parameters1, rhr_parameters1) = (precat20_parameters2, rhr_parameters2)
            (precat20_parameters2, rhr_parameters2) = generate_random(numitems, numhits, 100000 + (iteration*10000))

        elif  err1 <= 0.01 and err2 <= 0.01 :
            repeatflag = False
        elif iteration > 10 :
            print ("Warning : Not reached convergence after ", iteration, "iterations for ", size_data, size_hits)
            

##    outmap[x] = (precat20_parameters2, rhr_parameters2 , (err1,err2))
##    print (x, (precat20_parameters2, rhr_parameters2), (err1,err2))
    out = (precat20_parameters2, rhr_parameters2 , (err1,err2))
    print ('For () {},{} results is {}'.format(size_data, size_hits, out))
    print ("********************************************" )      

    return out    
