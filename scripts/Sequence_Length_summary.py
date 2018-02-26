#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 09:31:24 2018

@author: hhuang2
"""

import glob
import sqlite3 as sql
import csv
#from utils import IMGTdbIO

five_loci = ['A', 'B', 'C', 'DRB1', 'DQB1']

ID_file = "../Output/All_Case_IDs/IMGTv3310/Ten_ten_tx_seq.csv"
ID_list = ()
with open(ID_file, 'r') as csv_file:
    reader = csv.reader(csv_file)
    for row in reader:
        if row[0] != 'BMT':
            ID_list += (row[0],) # 18329

All_Db = "../Output/All_Case_IDs/IMGTv3310/fiveLoci_pairedSeq_Matching_stats_0128.db"

locus = 'DQB1'
db_39 = "../Output/SG39/2018/SG39_DRpairs/SG39_HLA_"+locus+"_paired.db"
db_41_52 = "../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpairs/SG41_52_HLA_"+locus+"_paired.db"

output_db = "../Output/All_Case_IDs/IMGTv3310/ten_ten_seqsLength_stats.db"

con1 = sql.connect(All_Db)
cur1 = con1.cursor()
Group = 'BatchGroup'

conn_save = sql.connect(output_db)
c = conn_save.cursor()
c.execute('''CREATE TABLE SequenceLengths_DQB1
          (BMT_caseID text, Donor_PS1 int,  Donor_PS2 int, Recipient_PS1 int, Recipient_PS2 int 
          )''')
#c.execute('''CREATE TABLE SequenceLengths_DRB1
##          (BMT_caseID text, Donor_PS1_Block1 int,  Donor_PS2_Block2 int, 
#          Recipient_PS1_Block1 int, Recipient_PS2_Block2 int 
#          )''')
for ID in ID_list: 
    ID = (ID,)
    cur1.execute('SELECT ' + Group + ' FROM MatchStats WHERE BMT_caseID = ?', ID)
    BatchGroup = cur1.fetchone()
    if BatchGroup != None:
        if BatchGroup[0] == 'SG39':
            db_fp = db_39
        else:
            db_fp = db_41_52
        
        con2 = sql.connect(db_fp)
        cur2 = con2.cursor()
        SeqsTypes = 'PS, Donor, Recipient'
        cur2.execute('SELECT ' + SeqsTypes + ' FROM OriginalSeqs WHERE BMT_caseID = ?', ID)
        Seqs = cur2.fetchall()
        
        con2.close()
    
        record = ID + (len(Seqs[0][1]), len(Seqs[1][1]), len(Seqs[0][2]), len(Seqs[1][2]), )
        c.execute('INSERT INTO SequenceLengths_DQB1 VALUES (?,?,?,?,?)', record)
    
conn_save.commit()

conn_save.close()

con1.close()

#AlgnFiles = glob.glob(singleMismatch_fp+'CaseID_'+key+'_Locus_'+locus+'*.pkl')