#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 12:36:48 2017

@author: hhuang2
"""

from utils import preprocessing as pp
import glob
import re
#import sqlite3 as sql

__author__ = "Hu Huang"
__copyright__ = "Copyright 2018, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

##
BMTcaseInfo_fp = "../../rawData/SG41-52_caseID.csv"
logfile_fp = "../Output/All_Cases_check/logfiles_01152018/"

for groupID in range(41,53):
    all_files = glob.glob("../../rawData/HLA_whole_gene/SG"+str(groupID)+"/reformatted/*.xls")
    for fp_file in all_files:
        logfile_name = logfile_fp + re.sub(".xls", "_log.txt", fp_file.split("/")[-1])
        seq_obj = pp.load_seq_file_newFormat(fp_file, BMTcaseInfo_fp, logfile_name)
        pp.saveAsSQLdb(seq_obj, "../Output/All_Cases_check/", "All", fp_file)
        
