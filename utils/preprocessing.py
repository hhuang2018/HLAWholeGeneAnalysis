#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
     
"""
#from collections import defaultdict
#import re
# import pickle # import csv
import sys
from utils import phase_block_check
import sqlite3
import csv
from xlrd import open_workbook

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

###############################################################################
## Preprocessing functions
###############################################################################
def isDonor(NMDP_ID):
    """
    Check if the NMDP_ID is a donor or a recipient
    Donor: 4-4-1
    Recipient: 3-3-1
    """
    if NMDP_ID.index("-") == 4:
        return("D")
    else:
        return("R")
    
def readBMTinfo(fp, header = True):
    """
    Read BMT case information into a list
    header - flag true if the table has a header; then will remove the header
    """
    # fp = "../../rawData/SG39_caseID.csv"
    
    caseIDs = {"BMTcase":[], "NMDP_DID":[], "NMDP_RID":[], "D_Audit":[], "D_Active":[], "D_Comment":[], 
               "R_Audit":[], "R_Active":[], "R_Comment":[]}
    with open(fp, 'r') as f:
        if header:
            next(f) # skip header
        reader=csv.reader(f)
        for line in reader:
           if line[0] not in caseIDs["BMTcase"]: # new record
               if line[2] == "D": # donor
                   caseIDs["BMTcase"].append(line[0])
                   caseIDs["NMDP_DID"].append(line[1])
                   caseIDs["D_Audit"].append(line[4])
                   caseIDs["D_Active"].append(line[5])
                   caseIDs["D_Comment"].append(line[6])
                   caseIDs["NMDP_RID"].append(" ")
                   caseIDs["R_Audit"].append(" ")
                   caseIDs["R_Active"].append(" ")
                   caseIDs["R_Comment"].append(" ")
                   
               elif line[2] == "R": # Recipient
                   caseIDs["BMTcase"].append(line[0])
                   caseIDs["NMDP_RID"].append(line[1])
                   caseIDs["R_Audit"].append(line[4])
                   caseIDs["R_Active"].append(line[5])
                   caseIDs["R_Comment"].append(line[6])
                   caseIDs["NMDP_DID"].append(" ")
                   caseIDs["D_Audit"].append(" ")
                   caseIDs["D_Active"].append(" ")
                   caseIDs["D_Comment"].append(" ")
                   
           elif line[2] == "R": # if the new record is recipient
               index = caseIDs["BMTcase"].index(line[0])
               caseIDs["NMDP_RID"][index] = line[1]
               caseIDs["R_Audit"][index] = line[4]
               caseIDs["R_Active"][index] = line[5]
               caseIDs["R_Comment"][index] = line[6]
               
           elif line[2] == "D": # new donor record
               index = caseIDs["BMTcase"].index(line[0])
               caseIDs["NMDP_DID"][index] = line[1]
               caseIDs["D_Audit"][index] = line[4]
               caseIDs["D_Active"][index] = line[5]
               caseIDs["D_Comment"][index] = line[6]
    # find un matched IDs           
    index_D = [i for i, x in enumerate(caseIDs["NMDP_DID"]) if x == " "]  # missing 2 donors
    [caseIDs["BMTcase"][i] for i in index_D]
    index_R = [i for i, x in enumerate(caseIDs["NMDP_RID"]) if x == " "]  # missing 0 donors
    [caseIDs["BMTcase"][i] for i in index_R]
    
    index_NAudit = [i for i, x in enumerate(caseIDs["D_Audit"]) if x != "Y"]      # 151 donors
    index_NAudit.append([i for i, x in enumerate(caseIDs["R_Audit"]) if x != "Y"]) # 152 individuals
    
    index_NActive = [i for i, x in enumerate(caseIDs["D_Active"]) if x != "Y"] # 88 donors
    index_NActive.append([i for i, x in enumerate(caseIDs["R_Active"]) if x != "Y"]) # 89 indiduals
    
    return(caseIDs)
    
def isTableExist(cursor, tableName):
   
    cursor.execute("""
                   SELECT name 
                   FROM sqlite_master 
                   WHERE type='table' AND name=?;
                   """, (tableName, ))

    exists = bool(cursor.fetchone())
    
    return(exists)


def load_seq_file(fp, BMTcaseInfo_fp, log_file, file_format = "xls"):
    """
    Read txt format sequence data from HML, and covert into Dictionary structure
    Plus alignment
    """
    
 #   if path.exists(fp):
    if file_format == "txt":  ## Tab-delimited text files
        #fp = "data/test_data.txt"    # for test
        txt_lines = open(fp).readlines()
        header = txt_lines[0].rstrip().split("\t")
        seq_dict = []
        for LineIndex in range(1, len(txt_lines)):
            temp_line_content = txt_lines[LineIndex].rstrip().split("\t")
            temp_dict = {}
            for ind, Items in enumerate(header):
                temp_dict[Items] = temp_line_content[ind]
            seq_dict.append(temp_dict)
                
    elif  file_format == "csv": ## CSV format files -- need to check this segment
        with open(fp, 'r') as f:
            header = next(f) # first line header
            header = header.rstrip().split(",")
            seq_dict = []
            reader = csv.reader(f)
            for line in reader:
                temp_line_content = txt_lines[LineIndex].rstrip().split("\t")
                temp_dict = {}
                for ind, Items in enumerate(header):
                    temp_dict[Items] = temp_line_content[ind]
                seq_dict.append(temp_dict)
                    
    elif file_format == "xls": ## Excel (xls) format files
       # fp = "../../rawData/xls/File_6.xls"    # for test
        # file 4:  17,016 lines; file 5:   6,770 lines; file 6: 41,516 lines
        # flle 7:  34,026 lines; file 8:  48,914 lines; file 9: 27,020 lines
        # file 10: 29,864 lines; file 11: 36,898 lines
        
        wb = open_workbook(fp)
        seq_dict = []
        for s in wb.sheets():
            #print 'Sheet:',s.name
            for row in range(s.nrows):
                if row == 0: # header
                    header = [s.cell(row,col).value for col in range(s.ncols)]
                else:
                    temp_dict = {}
                    for col in range(s.ncols):
                        temp_dict[header[col]] = s.cell(row,col).value
                    seq_dict.append(temp_dict)      
    print("The file has " + str(len(seq_dict)) + "Lines of records.\n")
    
    #BMTcaseInfo_fp = "../../rawData/SG39_caseID.csv"
    CaseIDs = readBMTinfo(BMTcaseInfo_fp)  ## 3608 cases total
    
    # num_seqs = len(seq_dict)
    new_seq_dict = {}
    #### build original sequence dictionary
    for items in seq_dict:  
        #if len(new_seq_dict) == 0:  # first item in the dictionary
        #    new_seq_dict[items.get('NMDP_DID/_RID')] = {items.get('Locus'): {'GLstring':[items.get('GL-string')], 
        #                'Sequence':[items.get('Sequence')]}}
            #new_seq_dict.append(new_item)
            #counter += 1
        ## BMT case ID
        if items.get('NMDP_DID/_RID') in CaseIDs["NMDP_DID"]: # donor
            BMTcase = CaseIDs["BMTcase"][CaseIDs["NMDP_DID"].index(items.get('NMDP_DID/_RID'))]
            DRtype = "D"
            Audit = CaseIDs["D_Audit"][CaseIDs["NMDP_DID"].index(items.get('NMDP_DID/_RID'))]
            Active = CaseIDs["D_Active"][CaseIDs["NMDP_DID"].index(items.get('NMDP_DID/_RID'))]
            Comment = CaseIDs["D_Comment"][CaseIDs["NMDP_DID"].index(items.get('NMDP_DID/_RID'))]
        elif items.get('NMDP_DID/_RID') in CaseIDs["NMDP_RID"]: # recipient
            BMTcase = CaseIDs["BMTcase"][CaseIDs["NMDP_RID"].index(items.get('NMDP_DID/_RID'))]
            DRtype = "R"
            Audit = CaseIDs["R_Audit"][CaseIDs["NMDP_RID"].index(items.get('NMDP_DID/_RID'))]
            Active = CaseIDs["R_Active"][CaseIDs["NMDP_RID"].index(items.get('NMDP_DID/_RID'))]
            Comment = CaseIDs["R_Comment"][CaseIDs["NMDP_RID"].index(items.get('NMDP_DID/_RID'))]
        else:
            BMTcase = " "
            DRtype = " "
            Audit = " "
            Active = " "
            Comment = " "
        ##
        if items.get('NMDP_DID/_RID') in list(new_seq_dict.keys()): # if existing ID, then append
            if items.get('Locus') in list(new_seq_dict[items.get('NMDP_DID/_RID')].keys()): ## If the Locus is already exists, then append
                new_seq_dict[items.get('NMDP_DID/_RID')][items.get('Locus')]['GLstring'].append(items.get('GL-string'))
                new_seq_dict[items.get('NMDP_DID/_RID')][items.get('Locus')]['Sequence'].append(items.get('Sequence'))
                new_seq_dict[items.get('NMDP_DID/_RID')][items.get('Locus')]['block'].append(items.get('Block'))
                new_seq_dict[items.get('NMDP_DID/_RID')][items.get('Locus')]['phase'].append(items.get('Phase'))
            else: # if not a new loucs, then add a new record
                # temp_locus = {items.get('Locus'): {'GLstring':[items.get('GL-string')], 'Sequence':[items.get('Sequence')]}}
                new_seq_dict[items.get('NMDP_DID/_RID')][items.get('Locus')] = {'GLstring':[items.get('GL-string')], 'phase':[items.get('Phase')], 'block':[items.get('Block')], 'Sequence':[items.get('Sequence')]}
        else: # if it's a new ID, then add a new record
            '''new_item = {'NMDP_ID': items.get('NMDP_DID/_RID'),
                        'Locus': {items.get('Locus'): {'GLstring':[items.get('GL-string')], 
                        'Sequence':[items.get('Sequence')]}}}
            new_seq_dict.append(new_item)
            counter += 1'''
            new_seq_dict[items.get('NMDP_DID/_RID')] = {"BMTcase": BMTcase, "DRtype": DRtype, "Audit": Audit, "Active": Active, "Comment": Comment,
                                                         items.get('Locus'): {'GLstring':[items.get('GL-string')], 
                                  'phase':[items.get('Phase')], 'block':[items.get('Block')],
                                  'Sequence':[items.get('Sequence')]}}
    
    ### Correct phase set, and merge block
    corrected_seq_table = {}
    
    sys.stdout = open(log_file,'wt')
    
    for individual_ID, individual_seq in new_seq_dict.items():
        loci = list(individual_seq.keys())
        loci.remove("BMTcase")
        loci.remove("DRtype")
        loci.remove("Audit")
        loci.remove("Active")
        loci.remove("Comment")
        for locus in loci:

            #print(locus)
            if individual_seq[locus]['GLstring'][0] != "":
                corrected_typing = phase_block_check.check_seq_typing(individual_seq[locus]['GLstring'], individual_seq[locus]['Sequence'], individual_ID)
 
                if individual_ID in list(corrected_seq_table.keys()): # if existing ID, then append
                    corrected_seq_table[individual_ID][locus] = corrected_typing# "GLstring"], 'Sequence', 'phase', 'block'
                else: # if it's a new ID, then add a new record
                    corrected_seq_table[individual_ID] = {locus: corrected_typing}
            else: # missing GLstrings
                print("<><><><><><><><>\n"+"NMDP_ID: " + individual_ID + " at locus " +locus + " is missing proper GL-strings.\n<><><><><><><><>\n")
        corrected_seq_table[individual_ID]["BMTcase"] = individual_seq["BMTcase"]
        corrected_seq_table[individual_ID]["DRtype"] = individual_seq["DRtype"]
        corrected_seq_table[individual_ID]["Audit"] = individual_seq["Audit"]
        corrected_seq_table[individual_ID]["Active"] = individual_seq["Active"]
        corrected_seq_table[individual_ID]["Comment"] = individual_seq["Comment"]
    print("Corrected Table has " + str(len(corrected_seq_table)) + " ID records.\n")
    # file 4:   614 IDs; file 5:   386 IDs;  file 6: 2,011* IDs; file 7:  921* IDs
    # file 8: 2,961 IDs; file 9: 1,638 IDs; file 10: 1,463 IDs;  file 11: 140* IDs
        
        
    # TODO: Build sequence database - by BMTcaseID, NMDP_ID; different locus saved in different files
    # saveAsSQLdb
    
    # TODO: Build and search Class II gene database. (-DRB1, DPB1, DQB1)
    # TODO: 
    return(corrected_seq_table) 

def saveAsSQLdb(seq_obj, output):
    """
    Save record in a SQL database; each locus have one db file; 
    seq_obj: corrected sequence table
    In each db file, table1 - originalSequences; table2 - Exon/Intron; table3 - translated Protein sequence
    """
    # output = "../Output/"
    # BMTcaseInfo_fp = "../../rawData/SG39_caseID.csv"
    # keys: BMTcase, NMDP_DID, NMDP_RID, D_Audit, D_Active, D_comment, R_Audit, R_Active, R_comment
    for individual_ID, individual_seq in seq_obj.items():
        loci = list(individual_seq.keys())
        loci.remove("BMTcase")
        loci.remove("DRtype")
        loci.remove("Audit")
        loci.remove("Active")
        loci.remove("Comment")
        for locus in loci:
            filename = output + "SG39_HLA_" + locus + "_originalTB.db"
            
            # original sequence table
            conn = sqlite3.connect(filename) # automatically creates a file if doesn't exist
            cursor = conn.cursor()
            cursor.execute('''CREATE TABLE IF NOT EXISTS OriginalSeqs
                           (BMT_caseID text, NMDP_ID text, DRtype text, 
                           Audit text, Active text, Comment text,
                           HLATyping text, PS text, Block1 text, Block2 text)''')
            BMT_caseID = str(individual_seq["BMTcase"])
            NMDP_ID = str(individual_ID)
            DRtype = str(individual_seq["DRtype"]) # isDonor(NMDP_ID)
            Audit = str(individual_seq["Audit"])
            Active = str(individual_seq["Active"])
            Comment = str(individual_seq["Comment"])
            
            cursor.execute('SELECT count(*) FROM OriginalSeqs WHERE NMDP_ID=?', (NMDP_ID, ))
            record_temp = cursor.fetchone() 
            if(record_temp[0] == 0): # if there is no record of this ID, then insert the record
                for PhaseID in individual_seq[locus].keys():
                    if len(individual_seq[locus][PhaseID]['blockIDs']) == 1:
                        HLATyping = str(individual_seq[locus][PhaseID]['GLstring'])
                        Block1 = str(individual_seq[locus][PhaseID]['Sequence'][0])
                        PS = str(PhaseID)
                        Block2 = ""
                        record = (BMT_caseID, NMDP_ID, DRtype, Audit, Active, Comment, HLATyping, PS, Block1, Block2,)
                        cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?,?)', record)
                        conn.commit()
                    elif len(individual_seq[locus][PhaseID]['blockIDs']) == 2:
                        for ind in range(2):
                            PS = str(PhaseID)
                            HLATyping = str(individual_seq[locus][PhaseID]['GLstring'])
                            if individual_seq[locus][PhaseID]['blockIDs'][ind] == 1:
                                Block1 = str(individual_seq[locus][PhaseID]['Sequence'][ind])
                            else:
                                Block2 = str(individual_seq[locus][PhaseID]['Sequence'][ind])
                        record = (BMT_caseID, NMDP_ID, DRtype, Audit, Active, Comment, HLATyping, PS, Block1, Block2,)
                        cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?,?)', record)
                        conn.commit()
                    else:
                        Block1 = ''
                        Block2 = ''
                        for ind in range(2):
                            PS = str(PhaseID)
                            HLATyping = str(individual_seq[locus][PhaseID]['GLstring'])
                            if individual_seq[locus][PhaseID]['blockIDs'][ind] == 1:
                                Block1 = Block1 + '*****' + str(individual_seq[locus][PhaseID]['Sequence'][ind])
                            else:
                                Block2 = Block2 + '*****' +str(individual_seq[locus][PhaseID]['Sequence'][ind])
                        record = (BMT_caseID, NMDP_ID, DRtype, Audit, Active, Comment, HLATyping, PS, Block1, Block2,)
                        cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?,?)', record)
                        conn.commit()
            conn.close()
 
def saveExonIntronInfo():
    
    # Extract Exon/Intron 
    conn2 = sqlite3.connect(filename) # automatically creates a file if doesn't exist
    cursor2 = conn2.cursor()
    if locus in ['A', 'C']: # 8 exons, 7 introns, 5'-UTR, 3'-UTR
        cursor2.execute('''CREATE TABLE IF NOT EXISTS ExonIntron
                        (BMT_caseID text, NMDP_ID text, DRtype text, 
                         Audit text, Active text, Comment text,
                         HLATyping text, PS text, 
                         five_prime_UTR text, Exon1 text, Intron1 text, Exon2 text, Intron2 text,
                         Exon3 text, Intron3 text, Exon4 text, Intron4 text, Exon5 text, Intron5 text,
                         Exon6 text, Intron6 text, Exon7 text, Intron7 text, Exon8 text, three_prime_UTR text)''')
    elif locus == 'B': # 7 exons, 6 introns, 5'-UTR, 3'-UTR  
        cursor2.execute('''CREATE TABLE IF NOT EXISTS ExonIntron
                        (BMT_caseID text, NMDP_ID text, DRtype text, 
                         Audit text, Active text, Comment text,
                         HLATyping text, PS text, 
                         five_prime_UTR text, Exon1 text, Intron1 text, Exon2 text, Intron2 text,
                         Exon3 text, Intron3 text, Exon4 text, Intron4 text, Exon5 text, Intron5 text,
                         Exon6 text, Intron6 text, Exon7 text, three_prime_UTR text)''')
    elif locus == 'DQB1': # Intron1-Exon2-Intron2-Exon3-Intron3-Exon4-Intron4
        cursor2.execute('''CREATE TABLE IF NOT EXISTS ExonIntron
                        (BMT_caseID text, NMDP_ID text, DRtype text,
                         Audit text, Active text, Comment text,
                         HLATyping text, PS text, 
                         Intron1 text, Exon2 text, Intron2 text,
                         Exon3 text, Intron3 text, Exon4 text, Intron4 text)''')
    elif locus in ['DPB1', 'DRB1']: # Intron1-Exon2-Intron2 and intron2-Exon3-Intron3
        cursor2.execute('''CREATE TABLE IF NOT EXISTS ExonIntron
                        (BMT_caseID text, NMDP_ID text, DRtype text, 
                         Audit text, Active text, Comment text,
                         HLATyping text, PS text, 
                         Intron1 text, Exon2 text, Intron2 text,
                         Exon3 text, Intron3 text)''')

    cursor2.execute('SELECT count(*) FROM ExonIntron WHERE NMDP_ID=?', (NMDP_ID, ))
    record_temp2 = cursor2.fetchone() 
    if(record_temp2[0] == 0):
        
        ExonIntronSeqs_ps = {}
        
        if len(individual_seq[locus]['phase']) == 1: ## homozygous
            ExonIntronSeqs_ph1 = findExonIntron(HLATypings = individual_seq[locus]['GLstring'][0], 
                                                Sequence = individual_seq[locus]['Sequence'], 
                                                Blocks = individual_seq[locus]['block'])
            ExonIntronSeqs_ps["PS1"] = ExonIntronSeqs_ph1 
            ExonIntronSeqs_ps["PS2"] = ExonIntronSeqs_ph1 
        elif len(individual_seq[locus]['phase'])== 2: ## two phases
            
            for ps in range(2):
                index_ph = [indexI for indexI, value in enumerate(individual_seq[locus]["phase"]) if value == str(ps + 1)]
                
                #   ExonIntronSeqs_ps["PS" + str(ps + 1)] = findExonIntron(HLATypings = individual_seq[locus]['GLstring'][ps], 
                #                                                         Sequence = [individual_seq[locus]['Sequence'][x] for x in index_ph], 
                #                                                         Blocks = individual_seq[locus]['block'][x] for x in index_ph)
                    
            for PhaseID in range(2):  
            
                if locus in ['A', 'C']: # 8 exons, 7 introns, 5'-UTR, 3'-UTR
                    record = (BMT_caseID, NMDP_ID, DRtype, 
                              Audit, Active, Comment, HLATyping, PS, 
                              ExonIntronSeqs["five_prime_UTR"], ExonIntronSeqs["Exon1"],
                              ExonIntronSeqs["Intron1"], ExonIntronSeqs["Exon2"],
                              ExonIntronSeqs["Intron2"], ExonIntronSeqs["Exon3"],
                              ExonIntronSeqs["Intron3"], ExonIntronSeqs["Exon4"],
                              ExonIntronSeqs["Intron4"], ExonIntronSeqs["Exon5"],
                              ExonIntronSeqs["Intron5"], ExonIntronSeqs["Exon6"],
                              ExonIntronSeqs["Intron6"], ExonIntronSeqs["Exon7"],
                              ExonIntronSeqs["Intron7"], ExonIntronSeqs["Exon8"],
                              ExonIntronSeqs["three_prime_UTR"],) ### 25 items
                    cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?)', record)
        
                elif locus == 'B': # 7 exons, 6 introns, 5'-UTR, 3'-UTR  
                    record = (BMT_caseID, NMDP_ID, DRtype, 
                              Audit, Active, Comment, HLATyping, PS, 
                              ExonIntronSeqs["five_prime_UTR"], ExonIntronSeqs["Exon1"],
                              ExonIntronSeqs["Intron1"], ExonIntronSeqs["Exon2"],
                              ExonIntronSeqs["Intron2"], ExonIntronSeqs["Exon3"],
                              ExonIntronSeqs["Intron3"], ExonIntronSeqs["Exon4"],
                              ExonIntronSeqs["Intron4"], ExonIntronSeqs["Exon5"],
                              ExonIntronSeqs["Intron5"], ExonIntronSeqs["Exon6"],
                              ExonIntronSeqs["Intron6"], ExonIntronSeqs["Exon7"],
                              ExonIntronSeqs["three_prime_UTR"],) ## 23 items
                    cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?, ?,?,?)', record)
                elif locus == 'DQB1': # Intron1-Exon2-Intron2-Exon3-Intron3-Exon4-Intron4
                    record = (BMT_caseID, NMDP_ID, DRtype, 
                              Audit, Active, Comment, HLATyping, PS, 
                              ExonIntronSeqs["Intron1"], ExonIntronSeqs["Exon2"],
                              ExonIntronSeqs["Intron2"], ExonIntronSeqs["Exon3"],
                              ExonIntronSeqs["Intron3"], ExonIntronSeqs["Exon4"],
                              ExonIntronSeqs["Intron4"],)  ## 15 items
                    cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?, ?,?,?,?,?, ?,?,?,?,?)', record)
                elif locus in ['DPB1', 'DRB1']: # Intron1-Exon2-Intron2 and intron2-Exon3-Intron3
                    record = (BMT_caseID, NMDP_ID, DRtype, 
                              Audit, Active, Comment, HLATyping, PS, 
                              ExonIntronSeqs["Intron1"], ExonIntronSeqs["Exon2"],
                              ExonIntronSeqs["Intron2"], ExonIntronSeqs["Exon3"],
                              ExonIntronSeqs["Intron3"], ) ## 13 items
                    cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?, ?,?,?,?,? ,?,?,?)', record)
        
                conn2.commit()
        conn2.close()
