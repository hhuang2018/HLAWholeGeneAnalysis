#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
     
"""
import glob
import sqlite3 as sql
from utils import phase_block_check as ps
from utils import IMGTdbIO

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

###############
# Class I
###############
#locus = "A"
all_DB_files = glob.glob("../Output/SG41_52/AvailDB/*.db")
db_file = all_DB_files[3] ## 1:A, 3:C, 4:B

locus = db_file.split('_')[4]

conn = sql.connect(db_file) # automatically creates a file if doesn't exist
conn.row_factory = sql.Row  # Each row is a dictionary: {colNames: Value}
cursor = conn.cursor()
#cursor.execute('SELECT * FROM OriginalSeqs')

All_caseID_cursor = cursor.execute('SELECT BMT_caseID FROM OriginalSeqs')
All_caseIDs_tuples = All_caseID_cursor.fetchall()

All_caseIDs = [item[0] for item in All_caseIDs_tuples]
All_caseIDs = list(set(All_caseIDs))

print(db_file+" has " + str(len(All_caseIDs)) + " cases.")

BMTcaseInfo_fp = "../../rawData/SG41-52_caseID.csv"
BMT_IDtable = ps.readBMTinfo(BMTcaseInfo_fp)

logf = open("../Output/SG41_52/logfiles_preprocessing/HLA_"+ locus +"_preprocess.txt", "w")
logf.write(db_file+" has " + str(len(All_caseIDs)) + " cases.\n\n")
available_records = {}
counter = 0
for BMTcase in All_caseIDs:
    query_c = cursor.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID=?', (BMTcase, ))
    all_query = query_c.fetchall()
    
    Donor = []
    Recipient = []
    for line in all_query:
                
        if line['DRtype'] == 'D':
            Donor.append({'HLATyping': line['HLATyping'], 'Block1':line['Block1'], 'PS': line['PS'], 'Block2': line['Block2']})
        elif line['DRtype'] == 'R':
            Recipient.append({'HLATyping': line['HLATyping'], 'Block1': line['Block1'], 'PS': line['PS'], 'Block2': line['Block2']})
    
    if len(Donor)>2 :
        #print(BMTcase + " has more donor sequences. Please double check!")
        logf.write(BMTcase + " has more donor sequences. Please double check!\n")
        counter += 1
        #break
    if len(Recipient)>2 :
        #print(BMTcase + " has more recipient sequences. Please double check!")
        logf.write(BMTcase + " has more recipient sequences. Please double check!\n")
        counter += 1
        #break
        
    if len(Donor)==1 and (Donor[0]['PS'] == 'PS3' or Donor[0]['Block2'] != ''):
        #print(BMTcase + " has more recipient sequences. Please double check!")
        logf.write(BMTcase + " ONLY has targeted sequences. Please double check!\n")
        counter += 1
        Donor = []
        #break
    if len(Recipient)==1 and (Recipient[0]['PS'] == 'PS3' or Recipient[0]['Block2'] != ''):
        #print(BMTcase + " has more recipient sequences. Please double check!")
        logf.write(BMTcase + " ONLY has targeted sequences. Please double check!\n")
        counter += 1
        Recipient = []
        #break
    if len(Donor) != 0 and len(Recipient) != 0:
        if BMTcase not in available_records.keys():
                available_records[BMTcase] = {'Active': line['Active'], 'Audit':line['Audit'], 'Comment': line['Comment'], 'QC': ''}

        if Recipient[0]['HLATyping'] == Donor[0]['HLATyping']: # first pair 1 vs 1
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[0]['HLATyping'], 'Recipient': Recipient[0]['Block1'], 'Donor': Donor[0]['Block1']}
            if Recipient[0]['Block2'].isalpha() or Donor[0]['Block2'].isalpha(): 
                if Recipient[0]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq'
                    
                if  Donor[0]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
            else: 
                available_records[BMTcase]['QC'] = 'PASS'
        
            if Recipient[1]['HLATyping'] == Donor[1]['HLATyping']: # second pair
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[1]['Block1']}
                if Recipient[1]['Block2'].isalpha() or Donor[1]['Block2'].isalpha(): 
                    if Recipient[1]['Block2'].isalpha():
                        available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq'
                        
                    if  Donor[1]['Block2'].isalpha():
                        available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
                else: 
                    available_records[BMTcase]['QC'] = 'PASS'
            else:
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'] + "+" + Donor[1]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[1]['Block1']}
                available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
                
        elif Recipient[0]['HLATyping'] == Donor[1]['HLATyping']: # first pair: 1 vs 2
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[0]['HLATyping'], 'Recipient': Recipient[0]['Block1'], 'Donor': Donor[1]['Block1']}
            if Recipient[0]['Block2'].isalpha() or Donor[1]['Block2'].isalpha(): 
                if Recipient[0]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq;'
                
                if  Donor[1]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
            else: 
                available_records[BMTcase]['QC'] = 'PASS'
            
            if Recipient[1]['HLATyping'] == Donor[0]['HLATyping']: # second pair
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[0]['Block1']}
                if Recipient[1]['Block2'].isalpha() or Donor[0]['Block2'].isalpha(): 
                    if Recipient[1]['Block2'].isalpha():
                        available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq;'
                        
                    if  Donor[0]['Block2'].isalpha():
                        available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
                else: 
                    available_records[BMTcase]['QC'] = 'PASS'
            else:
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'] + "+" + Donor[0]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[0]['Block1']}
                available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
        
        ###
        elif Recipient[1]['HLATyping'] == Donor[0]['HLATyping']: # first pair: 2 vs 1
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[0]['Block1']}
            if Recipient[1]['Block2'].isalpha() or Donor[0]['Block2'].isalpha(): 
                if Recipient[1]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq;'
                
                if  Donor[0]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
            else: 
                available_records[BMTcase]['QC'] = 'PASS'
            
            if Recipient[0]['HLATyping'] == Donor[1]['HLATyping']: # second pair
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[0]['HLATyping'], 'Recipient': Recipient[0]['Block1'], 'Donor': Donor[1]['Block1']}
                if Recipient[0]['Block2'].isalpha() or Donor[1]['Block2'].isalpha(): 
                    if Recipient[0]['Block2'].isalpha():
                        available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq;'
                        
                    if  Donor[1]['Block2'].isalpha():
                        available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
                else: 
                    available_records[BMTcase]['QC'] = 'PASS'
            else:
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[0]['HLATyping'] + "+" + Donor[1]['HLATyping'], 'Recipient': Recipient[0]['Block1'], 'Donor': Donor[1]['Block1']}
                available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
        
        ###
        elif Recipient[1]['HLATyping'] == Donor[1]['HLATyping']: # first pair: 2 vs 2
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[1]['Block1']}
            if Recipient[1]['Block2'].isalpha() or Donor[1]['Block2'].isalpha(): 
                if Recipient[1]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = '; R:Unexpected Block2 Seq;'
                
                if  Donor[1]['Block2'].isalpha():
                    available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + '; D:Unexpected Block2 Seq'
            else: 
                available_records[BMTcase]['QC'] = 'PASS'
            
            available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[0]['HLATyping']+ "+" +Donor[0]['HLATyping'], 'Recipient': Recipient[0]['Block1'], 'Donor': Donor[0]['Block1']}
            available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
    
        else:
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[0]['HLATyping']+ "+" +Donor[0]['HLATyping'], 'Recipient': Recipient[0]['Block1'], 'Donor': Donor[0]['Block1']}
            available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
            
            available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping']+ "+" +Donor[1]['HLATyping'], 'Recipient': Recipient[1]['Block1'], 'Donor': Donor[1]['Block1']}
            available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
    else:
        #print('Case ID ' + BMTcase + ' is missing matching donor or recipient info.' )
        
        BMTtb_index = BMT_IDtable['BMTcase'].index(BMTcase)
        if len(Donor) == 0:
            logf.write('Case ID\t' + BMTcase + ': Donor Info (DID: ' + BMT_IDtable['NMDP_DID'][BMTtb_index] + ' )is missing.\n')
            counter += 1
        if len(Recipient) == 0:
            logf.write('Case ID\t' + BMTcase + ': Recipient Info(RID '+ BMT_IDtable['NMDP_RID'][BMTtb_index] +' ) is missing.\n')
            counter += 1
logf.write('\nTotal missing values:\t' + str(counter) + '\n')
logf.close()

conn.close()

output = "../Output/SG41_52/SG41_52_DRpairs/"
for BMTcase, case_SeqInfo in available_records.items():
    filename = output + "SG41_52_HLA_" + locus + "_paired.db"
    # original sequence table
    conn = sql.connect(filename) # automatically creates a file if doesn't exist
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS OriginalSeqs
                   (BMT_caseID text, Audit text, Active text, Comment text,
                    QC text, HLATyping text, PS text, Donor text, Recipient text)''')
    BMT_caseID = str(BMTcase)
    Audit = str(case_SeqInfo["Audit"])
    Active = str(case_SeqInfo["Active"])
    Comment = str(case_SeqInfo["Comment"])
    QC = str(case_SeqInfo['QC'])
    
    # PS1
    PS_set = ['PS1', 'PS2']
    for PS in PS_set:
        HLATyping = case_SeqInfo[PS]['HLATyping']
        Donor = case_SeqInfo[PS]['Donor']
        Recipient = case_SeqInfo[PS]['Recipient']
        record = (BMT_caseID, Audit, Active, Comment, QC, HLATyping, PS, Donor, Recipient, )
        cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?)', record)
        conn.commit()

conn.close()

fname = output + 'SG41_52_HLA_' + locus + '_paired'
IMGTdbIO.save_dict2pickle(available_records, fname)

#aa = IMGTdbIO.load_pickle2dict(fname, output)


################# 
# Class II
#################
all_DB_files = glob.glob("../Output/SG41_52/AvailDB/*.db")
db_file = all_DB_files[5] ## 0: DPB1 2:DRB1 5:DQB1

locus = db_file.split('_')[4]

conn = sql.connect(db_file) # automatically creates a file if doesn't exist
conn.row_factory = sql.Row  # Each row is a dictionary: {colNames: Value}
cursor = conn.cursor()
#cursor.execute('SELECT * FROM OriginalSeqs')

All_caseID_cursor = cursor.execute('SELECT BMT_caseID FROM OriginalSeqs')
All_caseIDs_tuples = All_caseID_cursor.fetchall()

All_caseIDs = [item[0] for item in All_caseIDs_tuples]
All_caseIDs = list(set(All_caseIDs))

print(db_file+" has " + str(len(All_caseIDs)) + " cases.")

BMTcaseInfo_fp = "../../rawData/SG41-52_caseID.csv"
BMT_IDtable = ps.readBMTinfo(BMTcaseInfo_fp)

available_records = {}
logf = open("../Output/SG39/HLA_"+ locus +"_preprocess.txt", "w")
logf.write(db_file+" has " + str(len(All_caseIDs)) + " cases.\n\n")

counter = 0
for BMTcase in All_caseIDs:
    query_c = cursor.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID=?', (BMTcase, ))
    all_query = query_c.fetchall()
    
    Donor = []
    Recipient = []
    for line in all_query:
                
        if line['DRtype'] == 'D':
            Donor.append({'HLATyping': line['HLATyping'], 'Block1':line['Block1'], 'PS': line['PS'], 'Block2': line['Block2']})
        elif line['DRtype'] == 'R':
            Recipient.append({'HLATyping': line['HLATyping'], 'Block1': line['Block1'], 'PS': line['PS'], 'Block2': line['Block2']})
    
    if len(Donor)>2 :
        #print(BMTcase + " has more donor sequences. Please double check!")
        logf.write(BMTcase + " has more donor sequences. Please double check!\n")
        counter += 1
        #break
    if len(Recipient)>2 :
        #print(BMTcase + " has more recipient sequences. Please double check!")
        logf.write(BMTcase + " has more recipient sequences. Please double check!\n")
        counter += 1
        #break
    
    if len(Donor)==1 and (Donor[0]['PS'] == 'PS3' or Donor[0]['Block1'] == 'No sequence for nomakes'):
        #print(BMTcase + " has more recipient sequences. Please double check!")
        logf.write(BMTcase + " ONLY has targeted sequences. Please double check!\n")
        counter += 1
        Donor = []
        #break
    if len(Recipient)==1 and (Recipient[0]['PS'] == 'PS3' or Recipient[0]['Block1'] == 'No sequence for nomakes'):
        #print(BMTcase + " has more recipient sequences. Please double check!")
        logf.write(BMTcase + " ONLY has targeted sequences. Please double check!\n")
        counter += 1
        Recipient = []
        
    if len(Donor) == 2 and len(Recipient) == 2:
        if BMTcase not in available_records.keys():
                available_records[BMTcase] = {'Active': line['Active'], 'Audit':line['Audit'], 'Comment': line['Comment'], 'QC': ''}

        if Recipient[0]['HLATyping'] == Donor[0]['HLATyping']: # first pair 1 vs 1
            
            if Recipient[0]['Block2'].isalpha():
                Rec_seq = Recipient[0]['Block1'] + '----------' + Recipient[0]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[0]['Block1']
            if Donor[0]['Block2'].isalpha():
                Don_seq = Donor[0]['Block1'] + '----------' + Donor[0]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[0]['Block1']
                
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[0]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
            available_records[BMTcase]['QC'] = 'PASS'
        
            ## second pair       
            if Recipient[1]['Block2'].isalpha():
                Rec_seq = Recipient[1]['Block1'] + '----------' + Recipient[1]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[1]['Block1']
            if Donor[1]['Block2'].isalpha():
                Don_seq = Donor[1]['Block1'] + '----------' + Donor[1]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[1]['Block1']
                
            if Recipient[1]['HLATyping'] == Donor[1]['HLATyping']: # second pair
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
                available_records[BMTcase]['QC'] = 'PASS'
            else:
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'] + "+" + Donor[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
                available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
                
        elif Recipient[0]['HLATyping'] == Donor[1]['HLATyping']: # first pair: 1 vs 2
            
            if Recipient[0]['Block2'].isalpha():
                Rec_seq = Recipient[0]['Block1'] + '----------' + Recipient[0]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[0]['Block1']
                
            if Donor[1]['Block2'].isalpha():
                Don_seq = Donor[1]['Block1'] + '----------' + Donor[1]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[1]['Block1']
                
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[0]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
            
            available_records[BMTcase]['QC'] = 'PASS'
            
            # second pair
            if Recipient[1]['Block2'].isalpha():
                    Rec_seq = Recipient[1]['Block1'] + '----------' + Recipient[1]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[1]['Block1']
                
            if Donor[0]['Block2'].isalpha():
                Don_seq = Donor[0]['Block1'] + '----------' + Donor[0]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[0]['Block1']
                
            if Recipient[1]['HLATyping'] == Donor[0]['HLATyping']: # second pair
         
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
                available_records[BMTcase]['QC'] = 'PASS'
                
            else:
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping'] + "+" + Donor[0]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
                available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
        
        ###
        elif Recipient[1]['HLATyping'] == Donor[0]['HLATyping']: # first pair: 2 vs 1
            
            if Recipient[1]['Block2'].isalpha():
                Rec_seq = Recipient[1]['Block1'] + '----------' + Recipient[1]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[1]['Block1']
                
            if Donor[0]['Block2'].isalpha():
                Don_seq = Donor[0]['Block1'] + '----------' + Donor[0]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[0]['Block1']
            
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
           
            available_records[BMTcase]['QC'] = 'PASS'
            
            # second pair
            if Recipient[0]['Block2'].isalpha():
                Rec_seq = Recipient[0]['Block1'] + '----------' + Recipient[0]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[0]['Block1']
                
            if Donor[1]['Block2'].isalpha():
                Don_seq = Donor[1]['Block1'] + '----------' + Donor[1]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[1]['Block1']
                
            if Recipient[0]['HLATyping'] == Donor[1]['HLATyping']: # second pair
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[0]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
                available_records[BMTcase]['QC'] = 'PASS'
                
            else:
                available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[0]['HLATyping'] + "+" + Donor[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
                available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
        
        ###
        elif Recipient[1]['HLATyping'] == Donor[1]['HLATyping']: # first pair: 2 vs 2
            
            if Recipient[1]['Block2'].isalpha():
                Rec_seq = Recipient[1]['Block1'] + '----------' + Recipient[1]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[1]['Block1']
                
            if Donor[1]['Block2'].isalpha():
                Don_seq = Donor[1]['Block1'] + '----------' + Donor[1]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[1]['Block1']
                
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
            available_records[BMTcase]['QC'] = 'PASS'
            
            # second pair
            if Recipient[0]['Block2'].isalpha():
                Rec_seq = Recipient[0]['Block1'] + '----------' + Recipient[0]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[0]['Block1']
                
            if Donor[0]['Block2'].isalpha():
                Don_seq = Donor[0]['Block1'] + '----------' + Donor[0]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[0]['Block1']
                
            available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[0]['HLATyping']+ "+" +Donor[0]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
            available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
    
        else: # both pairs don't match
            if Recipient[0]['Block2'].isalpha():
                Rec_seq = Recipient[0]['Block1'] + '----------' + Recipient[0]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[0]['Block1']
                
            if Donor[0]['Block2'].isalpha():
                Don_seq = Donor[0]['Block1'] + '----------' + Donor[0]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[0]['Block1']
            available_records[BMTcase]['PS1'] = {'HLATyping': Recipient[0]['HLATyping']+ "+" +Donor[0]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
            available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
            
            if Recipient[1]['Block2'].isalpha():
                Rec_seq = Recipient[1]['Block1'] + '----------' + Recipient[1]['Block2'] # 10 '-' between two blocks
            else:
                Rec_seq = Recipient[1]['Block1']
                
            if Donor[1]['Block2'].isalpha():
                Don_seq = Donor[1]['Block1'] + '----------' + Donor[1]['Block2'] # 10 '-' between two blocks
            else:
                Don_seq = Donor[1]['Block1']
            available_records[BMTcase]['PS2'] = {'HLATyping': Recipient[1]['HLATyping']+ "+" +Donor[1]['HLATyping'], 'Recipient': Rec_seq, 'Donor': Don_seq}
            available_records[BMTcase]['QC'] = available_records[BMTcase]['QC'] + "; D-R pair GL-string doesn't match."
    
    elif len(Donor) == 0 or len(Recipient) == 0:
        #print('Case ID ' + BMTcase + ' is missing matching donor or recipient info.' )
        
        BMTtb_index = BMT_IDtable['BMTcase'].index(BMTcase)
        if len(Donor) == 0:
            logf.write('Case ID\t' + BMTcase + ': Donor Info (DID: ' + BMT_IDtable['NMDP_DID'][BMTtb_index] + ' )is missing.\n')
            counter += 1
        if len(Recipient) == 0:
            logf.write('Case ID\t' + BMTcase + ': Recipient Info(RID '+ BMT_IDtable['NMDP_RID'][BMTtb_index] +' ) is missing.\n')
            counter += 1
    elif len(Donor) == 1 or len(Recipient) == 1:
        
        BMTtb_index = BMT_IDtable['BMTcase'].index(BMTcase)
        if len(Donor) == 1:
            logf.write('Case ID\t' + BMTcase + ': Donor (DID: ' + BMT_IDtable['NMDP_DID'][BMTtb_index] + ' ) has only one PS.\n')
            counter += 1
        if len(Recipient) == 1:
            logf.write('Case ID\t' + BMTcase + ': Recipient (RID '+ BMT_IDtable['NMDP_RID'][BMTtb_index] +' ) has only one PS.\n')
            counter += 1
    else:
        logf.write('Case ID\t' + BMTcase + ' has irregular number of sequences. Please double check.\n')
        counter += 1
        

logf.write('\nTotal missing values:\t' + str(counter) + '\n')
logf.close()

conn.close()

output = "../Output/SG41_52/SG41_52_DRpairs/"
for BMTcase, case_SeqInfo in available_records.items():
    filename = output + "SG41_52_HLA_" + locus + "_paired.db"
    # original sequence table
    conn = sql.connect(filename) # automatically creates a file if doesn't exist
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE IF NOT EXISTS OriginalSeqs
                   (BMT_caseID text, Audit text, Active text, Comment text,
                    QC text, HLATyping text, PS text, Donor text, Recipient text)''')
    BMT_caseID = str(BMTcase)
    Audit = str(case_SeqInfo["Audit"])
    Active = str(case_SeqInfo["Active"])
    Comment = str(case_SeqInfo["Comment"])
    QC = str(case_SeqInfo['QC'])
    
    # PS1
    PS_set = ['PS1', 'PS2']
    for PS in PS_set:
        HLATyping = case_SeqInfo[PS]['HLATyping']
        Donor = case_SeqInfo[PS]['Donor']
        Recipient = case_SeqInfo[PS]['Recipient']
        record = (BMT_caseID, Audit, Active, Comment, QC, HLATyping, PS, Donor, Recipient, )
        cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?)', record)
        conn.commit()

conn.close()

fname = output + 'SG41_52_HLA_' + locus + '_paired'
IMGTdbIO.save_dict2pickle(available_records, fname)