# -*- coding: cp1252 -*-
# This codes take all files in a specified directory that have a epcifed prefix
# and extract selected rows with a specified prefix.
# we donot do the comelex paseing of the file_summary.py version that was too
# linekd to a specific formatting of he lines.

import os
import time
import datetime

##from file_handling_util import testexistencedirectory, testexistencefile, testexistenceanscreatedirectory
##from generic_util import printmap, printlist


##########################################
# get file properti such as time of creation, in the script.


def listfilestatus(path):

    filenames = []
#    path = os.path.join("C:\\","temp")
    for filename in os.listdir(path):
        info = os.stat(os.path.join(path,filename))
    #    print(filename, info.st_mtime)
        mtime = info.st_mtime
    #    print(filename, info)
        timestamp_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d-%H:%M')
        print (filename, timestamp_str)
        filenames.append(filename)
    return filenames

###
# tesing
if __name__ == "__main__":
    listfilestatus('.')
    listfilestatus('../Data_sets/SINGLE_PATIENT_Drug_modulation_data/Drug_modulation_data_C&P/colorectal-cancer-cp')
    
