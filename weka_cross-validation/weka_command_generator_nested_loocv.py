
import os
import sys
from os import listdir
from os.path import isfile, join
import argparse

mypath_1 = './autoweka_inputs'
#mypath_1 = "./seconda_generazione_merged_files_arff"
onlyfiles_1 = sorted([f for f in listdir(mypath_1) if isfile(join(mypath_1, f))])
print(onlyfiles_1)

NUMINDEXEDINPUTFILES = len(onlyfiles_1)
TRAINFILENAMES = onlyfiles_1
print('NUMINDEXEDINPUTFILES ', NUMINDEXEDINPUTFILES)

AUTOWEKADIRECTORY_INPUTS = os.path.join('.', 'autoweka_inputs')
AUTOWEKADIRECTORY_OUTPUTS = os.path.join('.', 'autoweka_nested_loocv_outputs') # cahnged output filder



evaluator_map = dict([        
                           ("cfssubseteval", "weka.attributeSelection.CfsSubsetEval -P 1 -E 1"),
                           ("correlation", "weka.attributeSelection.CorrelationAttributeEval"),
                           ("gainratio", "weka.attributeSelection.GainRatioAttributeEval"),
                           ("wrapperj48", "weka.attributeSelection.ClassifierAttributeEval -execution-slots 1 -B weka.classifiers.trees.J48 -F 5 -T 0.01 -R 1 -E DEFAULT -- -C 0.25 -M 2"),
                           ("wrapperj48se", 'weka.attributeSelection.ClassifierSubsetEval -B weka.classifiers.trees.J48 -T -H "Click to set hold out or test instances" -E DEFAULT -- -C 0.25 -M 2'),
                        ])


search_map = dict([
                           ("bestfirst", "weka.attributeSelection.BestFirst -D 1 -N 5"),
                           ("greedystepwise", "weka.attributeSelection.GreedyStepwise -T -1.7976931348623157E308 -N -1 -num-slots 1"),
                           ("ranker100", "weka.attributeSelection.Ranker -T -1.7976931348623157E308 -N 100"),

                        ])

assocaition_pairs = dict([ ('cfs_best', ("cfssubseteval","bestfirst")),
                               ('cfs_greedy', ("cfssubseteval","greedystepwise")),
                               ('corr_ranker', ("correlation","ranker100")),
                               ('gain_ranker', ("gainratio","ranker100")),
                               ('j48_ranker', ("wrapperj48","ranker100")),
                               ('j48_greedy', ("wrapperj48se","greedystepwise")),
                          ]  )


#####################

def count_pats(in_path, in_filename) :

    f = open(os.path.join(in_path, in_filename), 'r')

    line = f.readline()
    while not line.startswith('@DATA'):
        line =  f.readline()
    assert line.startswith('@DATA')                          
    line = f.readline()
    counter = 0
    while line != '' :
        counter += 1
        line = f.readline()

    print("number patients ", counter, in_filename)
    f.close()

    return counter


###########################################################

def mainfunction2(fileindex, action, effect, feature, machine, time_limit, fold_number) :

    if action == 'genmodel'   and feature == None :

        trainfilename  =  TRAINFILENAMES[fileindex]
        wekamodellog   = trainfilename.replace('.arff', '.modellog')
        wekapredictionslog =  trainfilename.replace('.arff', '.predictions.CSV')

        fulltrainname = os.path.join(AUTOWEKADIRECTORY_INPUTS ,trainfilename)
        fullwekamodellog=os.path.join(AUTOWEKADIRECTORY_OUTPUTS ,wekamodellog)
        fullwekapredictionslog = os.path.join(AUTOWEKADIRECTORY_OUTPUTS ,wekapredictionslog)
        print ( 'fulltrainname',  fulltrainname)
        print ('fullwekamodellog', fullwekamodellog)
        print('fullwekapredictionslog', fullwekapredictionslog)

        num_patients = count_pats(AUTOWEKADIRECTORY_INPUTS ,trainfilename)

        if fold_number != None :
            assert fold_number < num_patients ,  f"fold number greater than number of patietns  {num_patients }"
        else :
            fold_number = num_patients

        
#        s1 = "java -classpath ../weka/weka-3-8-3/weka.jar weka.Run "
        s1 = "java -classpath ~/weka/weka-3-8-7/weka.jar weka.Run "
        
 #       s2 = "-match-anywhere Auto "
        s2 =  " weka.classifiers.meta.AutoWEKAClassifier "
        
        s3 = " -t "+fulltrainname+" "

        s31 = '  -x  ' + str(fold_number) + '  -s 1234  -k  '
        
#        s4 = " -timeLimit 120 "
        s4 = " -timeLimit " + str(time_limit)  +  " "
        
        s5 = " -metric errorRate  "
        s6 = " -memLimit 3072 "
 #       s7 = " -d "+fullmodelname+" "
        s7 = " -seed 123 "
#        s8 = " -no-cv "
        s8 = ' -classifications "weka.classifiers.evaluation.output.prediction.CSV -distribution -decimals 12 -file ' + fullwekapredictionslog + '"'
        
        s9 = " > "+fullwekamodellog+" 2>&1  & "

        commandstring = s1 +s2 +s3 +s31 + s4 + s5 + s6 + s7 + s8 +s9


    if action == 'genmodel'   and feature != None :

        trainfilename  =  TRAINFILENAMES[fileindex]
        wekamodellog   = trainfilename.replace('.arff', '_' + feature + '.modellog')
        wekapredictionslog =  trainfilename.replace('.arff',  '_' + feature +  '.predictions.CSV')

        fulltrainname = os.path.join(AUTOWEKADIRECTORY_INPUTS ,trainfilename)
        fullwekamodellog=os.path.join(AUTOWEKADIRECTORY_OUTPUTS ,wekamodellog)
        fullwekapredictionslog = os.path.join(AUTOWEKADIRECTORY_OUTPUTS ,wekapredictionslog)
        print ( 'fulltrainname',  fulltrainname)
        print ('fullwekamodellog', fullwekamodellog)
        print('fullwekapredictionslog', fullwekapredictionslog)

        num_patients = count_pats(AUTOWEKADIRECTORY_INPUTS ,trainfilename)

        if fold_number != None :
            assert fold_number < num_patients ,  f"fold number greater than number of patietns  {num_patients }"
        else :
            fold_number = num_patients

        evaluationcomm= evaluator_map[assocaition_pairs[feature][0]]
        searchcomm=     search_map[assocaition_pairs[feature][1]]

        s1 = 'java -classpath ~/weka/weka-3-8-7/weka.jar weka.Run weka.classifiers.meta.AttributeSelectedClassifier '
        s2 = "  -t "+fulltrainname+" "
        s3 =   '  -x  ' + str(fold_number) + '  -s 1234  -k  '
        s31 = ' -classifications "weka.classifiers.evaluation.output.prediction.CSV -distribution -decimals 12 -file ' + fullwekapredictionslog + '"'       
        s4 =  " -E " + "'" + evaluationcomm +"' "
        s5 =  " -S " + "'" + searchcomm +"' "
        s6 = ' -W weka.classifiers.meta.AutoWEKAClassifier -- '
        s7 =   " -timeLimit " + str(time_limit)  +  " "
        s8 = " -metric errorRate  "
        s9 = " -memLimit 3072 "
        s10 = " -seed 123 "
        s11 = " > "+fullwekamodellog+" 2>&1  & "

        commandstring = s1 +s2 +s3 + s31 +s4 + s5 + s6 + s7 + s8 + s9 +s10 + s11

    elif action == 'eval_auc' :

        trainfilename  =  TRAINFILENAMES[fileindex]
        if feature == None :
            auc_eval_log_name = trainfilename.replace('.arff', '.auc')
            wekapredictionslog =  trainfilename.replace('.arff', '.predictions.CSV')
        else :
            auc_eval_log_name   = trainfilename.replace('.arff', '_' + feature + '.auc')
            wekapredictionslog =  trainfilename.replace('.arff',  '_' + feature +  '.predictions.CSV')
            
       

        fulltrainname = os.path.join(AUTOWEKADIRECTORY_INPUTS ,trainfilename)
        fullwekapredictionslog = os.path.join(AUTOWEKADIRECTORY_OUTPUTS ,wekapredictionslog)
        full_auc_eval_log_name=os.path.join(AUTOWEKADIRECTORY_OUTPUTS , auc_eval_log_name)

        print ( 'fulltrainname',  fulltrainname)
        print('fullwekapredictionslog', fullwekapredictionslog)
        print (' full_auc_eval_log _name',  full_auc_eval_log_name)


        num_patients = count_pats(AUTOWEKADIRECTORY_INPUTS ,trainfilename)


        s1 = "python3 nested_loocv_auc.py  "
        s2 = " --predictions " +  fullwekapredictionslog + " "
        s3 = " --arff " +  fulltrainname +  " "
        s31 = " --num_pat " + str(num_patients) + " "
        s4 = " --positive PD "
        s5 = " --output " +  full_auc_eval_log_name + " "

        commandstring = s1 + s2 + s3 + s31 + s4 + s5


    elif action == 'eval_conf_int' :

        trainfilename  =  TRAINFILENAMES[fileindex]
        if feature == None :
            auc_eval_log_name = trainfilename.replace('.arff', '.auc')
            auc_ci_eval_log_name = trainfilename.replace('.arff', '.ci_auc')
        else :
            auc_eval_log_name   = trainfilename.replace('.arff', '_' + feature + '.auc')
            auc_ci_eval_log_name = trainfilename.replace('.arff',  '_' + feature + '.ci_auc')
            
        
        full_auc_eval_log_name=os.path.join(AUTOWEKADIRECTORY_OUTPUTS , auc_eval_log_name)
        full_auc_ci_eval_log_name = os.path.join(AUTOWEKADIRECTORY_OUTPUTS ,auc_ci_eval_log_name)

        print (' full_auc_eval_log _name',  full_auc_eval_log_name)
        print(" full_auc_ci_eval_log_name", full_auc_ci_eval_log_name)
      
        s1 = "python3 delong_and_bootstrap_CI_computation.py "
        s2 = " --input " + full_auc_eval_log_name + " "
        s3 = " > " + full_auc_ci_eval_log_name +" 2>&1  & "

        commandstring = s1 + s2  + s3

# print the string
    if effect == 'print' :
        print ("Command to be executed:\n")
        print ( commandstring)
    
# execute the string
    if effect == 'execute' :
        print ("Sending to excution:\n")
        os.system("echo "+ "Begin: " + commandstring)
        os.system(commandstring)



##########################################################
#

def mainfunction():
    parser = argparse.ArgumentParser(description='Parse input param for the command to launch autoweka')

    # index of the file to be used
    parser.add_argument('fileindex',
                        action='store',
                        type=int,
                        choices=[x for x in range(NUMINDEXEDINPUTFILES)],
                        help='Index of the data set from zero to ' + str(NUMINDEXEDINPUTFILES -1))

    parser.add_argument('machine',
                        action='store',
#                        choices=["laptop", "pc64", "elmer"],
                        choices=["elmer"],
 #                       default = 'elmer',
                        help='Defines the directory  holding the input files in individual machines')


    parser.add_argument('action',
                        action='store',
  #                      required = True,
  #                       choices=['genfsdatatest2', 'genfsdata', 'genmodel','testmodel'],
  #                      choices=['genmodel','testmodel'],
                        choices=['genmodel', 'eval_auc', 'eval_conf_int'],
                        default=None,
                        help='Actions possible: apply feature selection, generate a model, apply a model')

    parser.add_argument(
        "--aw_time_limit",
#        required=True,
        type = int,
        default = 30,
        help="Time in minutes of one call to autoweka"
        )

    parser.add_argument(
        "--fold_number",
#        required=True,
        type = int,
        default = None,
        help="k for k-fold external validatiin, ouned by teh number of patients. Deault None = nubmer of patients"
        )


    parser.add_argument('--effect',
                        action='store',
                        choices=['print', 'execute'],
                        default='print',
                        help='effect possible: print command sting or execute it')


    parser.add_argument('--feature',
                        action='store',
                        choices=['cfs_best', 'cfs_greedy','corr_ranker','gain_ranker', 'j48_ranker', 'j48_greedy'],
                        default=None,
                        help='Select the feature selection mode')


    args = parser.parse_args()

    



    print( "fileindex", args.fileindex)
    print( "machine", args.machine)
    print( "action", args.action)
    print ('time limit for AW', args.aw_time_limit)
    print('fold number ', args.fold_number)
    print( "effect", args.effect)
    print( "feature", args.feature)

    mainfunction2(args.fileindex, args.action, args.effect, args.feature, args.machine, args.aw_time_limit, args.fold_number)     
 

###############################################
# main

if __name__ == "__main__":
    mainfunction()
    
