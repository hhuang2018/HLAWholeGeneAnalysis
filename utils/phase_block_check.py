#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequence validation and PS/BLOCK correction

HLA-A,-B,-C: 2 PS per sample, 1 Block per PS



Created on Fri Sep  1 15:30:32 2017

@author: hhuang2
"""
from utils import IMGTdbIO
import re
import csv
import difflib


__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

def isDonor(NMDP_ID):
    """
    Check if the NMDP_ID is a donor or a recipient
    Donor: 4-4-1
    Recipient: 3-3-1
    
    With the function readBMTinfo(), this function became obsolete. 
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
    
## sql database operation
def isTableExist(cursor, tableName):
   
    cursor.execute("""
                   SELECT name 
                   FROM sqlite_master 
                   WHERE type='table' AND name=?;
                   """, (tableName, ))

    exists = bool(cursor.fetchone())
    
    return(exists)

### Roll-up HLA gl-string
def rollupHLAtyping(HLAtyping, rmPrefix = True):
    """"
    Roll up 4-field gl-string into two-field representative gl-string; remove "HLA-" prefix
    """
    if rmPrefix:
        HLAtyping = re.sub("HLA-", "", HLAtyping)
        
    tplist = HLAtyping.split("/")
    rollupTypings = [tp.split(":")[0]+":"+tp.split(":")[1] for tp in tplist]
    rollupTypings = list(set(rollupTypings))
    
    if len(rollupTypings) == 1:
        rollupTypings = rollupTypings[0]

    return(rollupTypings)
    
def check_oneBlock_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID):
    '''
    For one block one phase sequence
    '''
    Locus =  tplist[0].split("*")[0] 
    
    ARS0seq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
    ARS1seq = IMGTdbIO.readIMGTsql(tplist[1], field='Exon2, Exon3')
        
    if seq_count > 2:
        print("Please check the ID: " + ID + " Locus " + Locus + "! More sequences than expected.")
          
    QueryTyping = {}
    for seq_item in unique_Query:
        
        if ARS0seq != ARS1seq: #  if the two types have different ARS regions
        
            if ARS0seq[0] in seq_item  and ARS0seq[1] in seq_item: # the first type
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS1"]['Sequence'].append(seq_item)
                    QueryTyping["PS1"]['blockIDs'].append(2)
            elif ARS1seq[0] in seq_item and ARS1seq[1] in seq_item: # second type
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS2"]['Sequence'].append(seq_item)
                    QueryTyping["PS2"]['blockIDs'].append(2)
            else: 
                if "PS3" not in QueryTyping.keys():
                        QueryTyping["PS3"] = {"GLstring": unique_HLATyping_list, "Sequence": [seq_item], "blockIDs": [1]}
                else:
                        QueryTyping["PS3"]['Sequence'].append(seq_item)
                        QueryTyping["PS3"]['blockIDs'].append(2)
                print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
        
        else: #  if the two types have the same ARS regions
            ARS0seq1456 = IMGTdbIO.readIMGTsql(tplist[0], field='Exon1, Exon4, Exon5, Exon6')
            ARS1seq1456 = IMGTdbIO.readIMGTsql(tplist[1], field='Exon1, Exon4, Exon5, Exon6')
            if ARS0seq1456 != ARS1seq1456:
                if ARS0seq1456[0] in seq_item  and ARS0seq1456[1] in seq_item and ARS0seq1456[2] in seq_item and ARS0seq1456[3] in seq_item: # the first type
                    if "PS1" not in QueryTyping.keys():
                        QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                    else:
                        QueryTyping["PS1"]['Sequence'].append(seq_item)
                        QueryTyping["PS1"]['blockIDs'].append(2)
                
                elif ARS1seq1456[0] in seq_item  and ARS1seq1456[1] in seq_item and ARS1seq1456[2] in seq_item and ARS1seq1456[3] in seq_item: # second type
                    if "PS2" not in QueryTyping.keys():
                        QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                    else:
                        QueryTyping["PS2"]['Sequence'].append(seq_item)
                        QueryTyping["PS2"]['blockIDs'].append(2)
                else: 
                    if "PS3" not in QueryTyping.keys():
                        QueryTyping["PS3"] = {"GLstring": unique_HLATyping_list, "Sequence": [seq_item], "blockIDs": [1]}
                    else:
                        QueryTyping["PS3"]['Sequence'].append(seq_item)
                        QueryTyping["PS3"]['blockIDs'].append(2)
                    print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
            else:
                ARS0seq7 = IMGTdbIO.readIMGTsql(tplist[0], field='Exon7')
                ARS1seq7 = IMGTdbIO.readIMGTsql(tplist[1], field='Exon7')
                if ARS0seq7 != ARS1seq7:
                    if ARS0seq7[0] in seq_item: # the first type
                        if "PS1" not in QueryTyping.keys():
                            QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                        else:
                            QueryTyping["PS1"]['Sequence'].append(seq_item)
                            QueryTyping["PS1"]['blockIDs'].append(2)
                
                    elif ARS1seq7[0] in seq_item: # second type
                        if "PS2" not in QueryTyping.keys():
                            QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                        else:
                            QueryTyping["PS2"]['Sequence'].append(seq_item)
                            QueryTyping["PS2"]['blockIDs'].append(2)
                    else: 
                        QueryTyping["PS3"] ={"GLstring": unique_HLATyping_list, "Sequence": [seq_item], "blockIDs": [1]}
                        print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
                else:## all 8 exons are the same
                    if "PS1" not in QueryTyping.keys():
                        QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                    elif "PS2" not in QueryTyping.keys(): 
                        QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                    else: 
                        QueryTyping["PS1"]['Sequence'].append(seq_item)
                        QueryTyping["PS1"]['blockIDs'].append(2)
                    print(ID + ": The sequence at Locus " + Locus + " two typings have exactly the same Exon sequences. Cannot distinguish by Exons.")

    if "PS1" in QueryTyping.keys() and "PS2" not in QueryTyping.keys(): ## Homozygous
            QueryTyping["PS2"] = QueryTyping["PS1"]
            
    return(QueryTyping)

def check_twoBlock_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID):
    '''
    Two blocks one phase sequences
    '''
    
    Locus =  tplist[0].split("*")[0] 
    ARS0seq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
    ARS1seq = IMGTdbIO.readIMGTsql(tplist[1], field='Exon2, Exon3')
        
    if seq_count > 4:
        print("Please check the ID: " + ID + " Locus " + Locus + "! More sequences than expected.")
    
    QueryTyping = {}
    for seq_item in unique_Query:
        if ARS0seq[0] in seq_item: # the first type; block 1; exon2
            if "PS1" not in QueryTyping.keys():
                QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
            else: # altered block order
                QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item, QueryTyping["PS1"]["Sequence"][0]], "blockIDs": [1,2]}
        elif ARS0seq[1] in seq_item: # the first type; block 2; exon3
            if "PS1" not in QueryTyping.keys():
                QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [2]}
            else:
                QueryTyping["PS1"]['Sequence'].append(seq_item)
                QueryTyping["PS1"]['blockIDs'].append(2)
                
        elif ARS1seq[0] in seq_item: # second type; block 1; exon2
            if "PS2" not in QueryTyping.keys():
                QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
            else:
                QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item, QueryTyping["PS2"]["Sequence"][0]], "blockIDs": [1,2]}
        elif ARS1seq[1] in seq_item: # second type; block2; exon3
            if "PS2" not in QueryTyping.keys():
                QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [2]}
            else:
                QueryTyping["PS2"]['Sequence'].append(seq_item)
                QueryTyping["PS2"]['blockIDs'].append(2)
        else: 
            QueryTyping["PS3"] ={"GLstring": unique_HLATyping_list, "Sequence": [seq_item], "blockIDs": [1]}
            print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
        
    if "PS1" in QueryTyping.keys() and "PS2" not in QueryTyping.keys(): ## Homozygous
        QueryTyping["PS2"] = QueryTyping["PS1"]
    
    return(QueryTyping)
    
def check_DQB102_Block_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID):
    '''
    Two blocks one phase sequences
    '''
    
    Locus =  tplist[0].split("*")[0] 
    ARS0seq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
    ARS1seq = IMGTdbIO.readIMGTsql(tplist[1], field='Exon2, Exon3')
        
    serotype = [tp.split(":")[0] for tp in tplist]
    
    if seq_count > 3:
        print("Please check the ID: " + ID + " Locus " + Locus + ", have heterozygotic DQB1*02 types or have more sequences than expected.")
    
    QueryTyping = {}
    for seq_item in unique_Query:
        # PS1
        if ARS0seq[0] in seq_item: # PS1 Exon 2
            if serotype[0] == "DQB1*02":  # DQB1*02 - 2 blocks
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                else: # altered block order
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item, QueryTyping["PS1"]["Sequence"][0]], "blockIDs": [1,2]}
            else:  # non-DQB1 - 1 block
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS1"]['Sequence'].append(seq_item)
                    QueryTyping["PS1"]['blockIDs'].append(2)
                    
        elif ARS0seq[1] in seq_item: # PS1 Exon 3
            if serotype[0] == "DQB1*02":  # DQB1*02 - 2 blocks
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [2]}
                else:
                    QueryTyping["PS1"]['Sequence'].append(seq_item)
                    QueryTyping["PS1"]['blockIDs'].append(2)
            else:  # non-DQB1 - 1 block
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS1"]['Sequence'].append(seq_item)
                    QueryTyping["PS1"]['blockIDs'].append(2)
                    
        ## PS2
        elif ARS1seq[0] in seq_item: # PS2 Exon 2
            if serotype[0] == "DQB1*02":  # DQB1*02 - 2 blocks
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                else: # altered block order
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item, QueryTyping["PS2"]["Sequence"][0]], "blockIDs": [1,2]}
            else:  # non-DQB1 - 1 block
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS2"]['Sequence'].append(seq_item)
                    QueryTyping["PS2"]['blockIDs'].append(2)
                    
        elif ARS1seq[1] in seq_item: # PS2 Exon 3
            if serotype[0] == "DQB1*02":  # DQB1*02 - 2 blocks
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [2]}
                else:
                    QueryTyping["PS2"]['Sequence'].append(seq_item)
                    QueryTyping["PS2"]['blockIDs'].append(2)
            else:  # non-DQB1 - 1 block
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS2"]['Sequence'].append(seq_item)
                    QueryTyping["PS2"]['blockIDs'].append(2)

        else: 
            QueryTyping["PS3"] ={"GLstring": unique_HLATyping_list, "Sequence": [seq_item], "blockIDs": [1]}
            print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
        
    if "PS1" in QueryTyping.keys() and "PS2" not in QueryTyping.keys(): ## Homozygous
        QueryTyping["PS2"] = QueryTyping["PS1"]
    
    return(QueryTyping)

def check_seq_typing(HLAtypings, Query, ID):
    '''
    Check sequence typing one by one. 
    A,B,C - one long sequence;
    DRB1, DPB1, DQB1*02: 2 blocks
    DQB1 the rest: one long sequence 
    
    Only need to check ARS region
    
    Input:
        HLAtyping - Both HLAtyping list; ambiguous ones only list the first one
        Query - one block sequence
        
    Output: 
        queryTyping: {"PS1": {"GLstring": HLAtyping1, "Sequence": Sequence1, "blockIDs": [1, 2]}, 
                      "PS2": {"GLstring": HLAtyping2, "Sequence": Sequence2, "blockIDs": [1, 2]}}
        blockIDs - corresponding to the sequence order.
    '''
    
    Locus =  HLAtypings[0].split("*")[0]
    Locus =  re.sub("HLA-", "", Locus)
    
    unique_Query = list(set(Query))
    
    unique_HLATyping = list(set(HLAtypings))
    unique_HLATyping_list = sum([items.split("+") for items in unique_HLATyping], []) # reduce(lambda x, y: x + y, [items.split("+") for items in HLAtypings], []) 

    unique_HLATyping_list = [re.sub("HLA-", "", tp) for tp in unique_HLATyping_list]
    ## remove ambiguous typings, only take the first one (more common type)
    # tplist = [re.sub("HLA-", "", tp) for tp in typing_list]
    tplist = [tp.split("/")[0] for tp in unique_HLATyping_list]
    serotype = [tp.split(":")[0] for tp in tplist]
    
    seq_count = len(unique_Query)
    
    QueryTyping = {}
    if Locus in ["A", "B", "C"]:
        '''ARS0seq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
        ARS1seq = IMGTdbIO.readIMGTsql(tplist[1], field='Exon2, Exon3')
        
        if seq_count > 2:
            print("Please check the ID: " + ID + " Locus " + Locus + "! More sequences than expected.")
          
        for seq_item in unique_Query:
            if ARS0seq[0] in seq_item  and ARS0seq[1] in seq_item: # the first type
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": seq_item, "blockIDs": [1]}
                else:
                    QueryTyping["PS1"]['Sequence'].append(seq_item)
                    QueryTyping["PS1"]['blockIDs'].append(2)
            elif ARS1seq[0] in seq_item and ARS1seq[1] in seq_item: # second type
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": seq_item, "blockIDs": [1]}
                else:
                    QueryTyping["PS2"]['Sequence'].append(seq_item)
                    QueryTyping["PS2"]['blockIDs'].append(2)
            else: 
                QueryTyping ={"Typing": "NotMatching", "index": -1, "Block": "NA"}
                print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
        
        if "PS1" in QueryTyping.keys() and "PS2" not in QueryTyping.keys(): ## Homozygous
            QueryTyping["PS2"] = QueryTyping["PS1"]
        '''
        if seq_count < len(unique_HLATyping_list):
            unique_HLATyping_list = list(set(unique_HLATyping_list))

        if seq_count == len(unique_HLATyping_list) and seq_count == 1: ## Homozygous
            unique_Query.append(unique_Query[0])
            unique_HLATyping_list.append(unique_HLATyping_list[0])
            seq_count = len(unique_Query)
        
        QueryTyping = check_oneBlock_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID)
        
    elif Locus in ["DRB1", "DPB1"]:
        ''' ARS0seq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
        ARS1seq = IMGTdbIO.readIMGTsql(tplist[1], field='Exon2, Exon3')
        
        if seq_count > 4:
            print("Please check the ID: " + ID + " Locus " + Locus + "! More sequences than expected.")
        
        QueryTyping = {}
        for seq_item in unique_Query:
            if ARS0seq[0] in seq_item: # the first type; block 1; exon2
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [1]}
                else: # altered block order
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item, QueryTyping["PS1"]["Sequence"]], "blockIDs": [1,2]}
            elif ARS0seq[1] in seq_item: # the first type; block 2; exon3
                if "PS1" not in QueryTyping.keys():
                    QueryTyping["PS1"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [2]}
                else:
                    QueryTyping["PS1"]['Sequence'].append(seq_item)
                    QueryTyping["PS1"]['blockIDs'].append(2)
            elif ARS1seq[0] in seq_item: # second type; block 1; exon2
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[1], "Sequence": [seq_item], "blockIDs": [1]}
                else:
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item, QueryTyping["PS2"]["Sequence"]], "blockIDs": [1,2]}
            elif ARS1seq[1] in seq_item: # second type; block2; exon3
                if "PS2" not in QueryTyping.keys():
                    QueryTyping["PS2"] = {"GLstring": unique_HLATyping_list[0], "Sequence": [seq_item], "blockIDs": [2]}
                else:
                    QueryTyping["PS2"]['Sequence'].append(seq_item)
                    QueryTyping["PS2"]['blockIDs'].append(2)
            else: 
                QueryTyping ={"Typing": "NotMatching", "index": -1, "Block": "NA"}
                print(ID + ": The sequence at Locus " + Locus + " doesn't match to either of the Typings")
        
        if "PS1" in QueryTyping.keys() and "PS2" not in QueryTyping.keys(): ## Homozygous
            QueryTyping["PS2"] = QueryTyping["PS1"]
        '''
        if seq_count == 2: ## Homozygous
            unique_Query.append(unique_Query[0])
            unique_Query.append(unique_Query[1])
            unique_HLATyping_list.append(unique_HLATyping_list[0])
            seq_count = len(unique_Query)
            
        QueryTyping = check_twoBlock_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID)
       
    elif Locus == "DQB1":
        
       if "DQB1*02" not in serotype:
           if seq_count < len(unique_HLATyping_list):
               unique_HLATyping_list = list(set(unique_HLATyping_list))
           if seq_count == len(unique_HLATyping_list) and seq_count == 1: ## Homozygous
               unique_Query.append(unique_Query[0])
               unique_HLATyping_list.append(unique_HLATyping_list[0])
               seq_count = len(unique_Query)
               
           QueryTyping = check_oneBlock_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID)
       else: 
           QueryTyping = check_DQB102_Block_seq(seq_count, tplist, unique_Query, unique_HLATyping_list, ID)
       
    return(QueryTyping)
       