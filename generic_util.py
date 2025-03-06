# generic utilities

import math

###########################################################
# generic utilities

def printmap(inmap, explanationstring, limit = None) :

    print ("*************************************************************")
    print (explanationstring, "size", len(inmap), 'with limit at', limit)
    counter = 0
    key_list = [x for x in inmap.keys()]
    key_list_sorted = sorted(key_list)
    for x in key_list_sorted :
        print (x, inmap[x])
        counter += 1
        if limit != None :
            if counter  > limit :
                break

def printlist(inlist, explanationstring, limit = None) :

    print ("*************************************************************")
    print (explanationstring, "size", len(inlist), 'with limit at', limit)
    counter = 0
    for x in range(len(inlist)) :
        print (x, inlist[x])
        counter += 1
        if limit != None :
            if counter  > limit :
                break

def venndiagramdrugs(strd1, d1, strd2, d2) :

    set1 = set(d1.keys())
    set2 = set(d2.keys())

    print ("In ", strd1 , " not in ", strd2 , set1 - set2)
    print ("In ", strd2 , " not in ", strd1 , set2 - set1)
    print ("In both ",  set1 & set2)
    

def minuslog10(y) :
    if y == 0.0 :
        return 1000
    else :
        return -math.log10(y)


########################################
# this function incremetn a couter associated to keys in a mapping

def increment_count(mapping, key, first_value = 0) :

    if key in mapping.keys() :
        mapping[key] += 1
    else :
        mapping[key] = first_value


#################################
# function ti restrict the domain of a dictionary
# subset policy assumes all keys in the list are in the dictionary, if not rise an error
#

def domain_restriction(restriction_list, in_dict, policy='subset') :

    assert len(restriction_list) == len(set(restriction_list)) # cehecks for absence of duplicate keys

    out_dict = {}
    for x in restriction_list :
        if x in in_dict :
            out_dict[x] = in_dict[x]
        else:
            if policy=='subset':
                print('Error key not in map domain: {}'.format(x))
                assert False
    
    return out_dict


    
