#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 12:36:48 2017

@author: hhuang2
"""

from utils import preprocessing as pp
import glob
import re
import sqlite3 as sql

__author__ = "Hu Huang"
__copyright__ = "Copyright 2018, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

#st_group = 'sg41'
BMTcaseInfo_fp = "../../rawData/SG41-52_caseID.csv"
logfile_fp = "../Output/SG41_52/2018/IMGTv3310/logfiles_01252018/"
for groupID in range(41,53):
    all_files = glob.glob("../../rawData/HLA_whole_gene/SG"+str(groupID)+"/reformatted/*.xls")
    for fp_file in all_files:
        logfile_name = logfile_fp + re.sub(".xls", "_log.txt", fp_file.split("/")[-1])
        seq_obj = pp.load_seq_file_newFormat(fp_file, BMTcaseInfo_fp, logfile_name, version = "3310")
        pp.saveAsSQLdb(seq_obj, "../Output/SG41_52/2018/IMGTv3310/originalDB_01252018/", "SG41_52")



    
## extract only cases with BMT IDs
    
all_DB_files = glob.glob("../Output//SG41_52/2018/IMGTv3310//originalDB_01252018/*.db")

for db_file in all_DB_files:
    conn = sql.connect(db_file) # automatically creates a file if doesn't exist
    cursor = conn.cursor()
    #cursor.execute('SELECT * FROM OriginalSeqs')
    
    available_records = []
    for row in cursor.execute('SELECT * FROM OriginalSeqs'):
        if row[0].isnumeric():
            available_records.append(row)  ## Remove empty BMT IDs
            
    conn.close()
    locus = db_file.split("_")[5]
    
    print("Number of Available individuals at locus " + locus + " is " + str(len(available_records)))
    print("Number of Available cases at locus "+locus+" is "+ str(len(available_records)/4))
    
    new_db_fp = "../Output/SG41_52/2018/IMGTv3310/AvailDB/SG41_52_HLA_"+ locus +"_avail.db"
    conn = sql.connect(new_db_fp)
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS OriginalSeqs
                   (BMT_caseID text, NMDP_ID text, DRtype text, 
                    Audit text, Active text, Comment text,
                    HLATyping text, PS text, Block1 text, Block2 text, File_Path text)''')
    cursor.executemany('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?,?,?)', available_records)
    conn.commit()
    conn.close()

### Check IDs with extra sequences
all_DB_files = glob.glob("../Output/SG39//*.db")
db_file = all_DB_files[3]
conn = sql.connect(db_file) # automatically creates a file if doesn't exist
cursor = conn.cursor()
    #cursor.execute('SELECT * FROM OriginalSeqs')
    
irregular_records = []
for row in cursor.execute('SELECT * FROM OriginalSeqs'):
    if row[-1].isalpha():
        irregular_records.append(row)
            
conn.close()
