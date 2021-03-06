#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
     
"""
from os import path, makedirs
from collections import defaultdict
import re
# import pickle # import csv
import sys
from utils import IMGTtools, IMGTdbIO
from Bio import SeqIO, pairwise2, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC #, generic_dna, generic_protein
#from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import subprocess 
# from StringIO import StringIO # Python 2
# try:
#    from StringIO import StringIO
# except ImportError:
#    from io import StringIO   # Python 3
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
## IMGT database related functions
###############################################################################

def find_IMGT_alignment(filename, HLAtyping, header_line = 8):
    """
    Find the aligned sequence from the database file
    """
    # Read alignment file
    alignment_list = IMGTdbIO.read_IMGT_alignment(filename, header_line)
    
    HLAtyping_pattern = re.compile(re.escape(HLAtyping))  ## pattern - find HLA typing
    
    #######  read alignement line
    aligned_seq = "" 
    for lines in alignment_list:
        if HLAtyping in HLAtyping_pattern.findall(lines): ## find the exact matching
            aligned_seq += re.sub(" ", "", re.sub(re.escape(HLAtyping), "", lines.rstrip()))
            break
        
    ## check if the HLA typing has the exact matching alignment
    if aligned_seq == "": ### if the HLA typing is not found in the list
        ##### find the
        aligned_typing_seq = re.sub("\.", "-", aligned_seq)  
        
    else:
        
        aligned_typing_seq = re.sub("\.", "-", aligned_seq)
    
    return(aligned_typing_seq)

def find_IMGT_sequence(file_fp, HLAtyping, locus):
    """
    Find the unaligned sequence from the database file
    """
    # Read fasta file
    filename = file_fp + locus + "_gen.fasta"
    seqs = list(SeqIO.parse(filename, "fasta"))
    
    HLAtyping = re.sub("HLA-", "", HLAtyping)
    HLAtyping_pattern = re.compile(re.escape(HLAtyping))
    index = 0
    for seq in seqs:
        if HLAtyping_pattern.search(seq.description):
            break;
        else: 
            index += 1
    return(seqs[index])        


def save_load_IMGTdb(HLA_locus = "A", fp='data/', mode = "r"):
    """
    """
    fname = 'HLA_' + HLA_locus + '_IMGTdb'
    if mode == "r" and path.exists(fp+fname): # if exist, load the data
        IMGTdbIO.load_IMGTdb(fname, fp)
    elif mode == "w":    # if not, then save the data
        IMGTdbIO.IMBTdb_2_dict(HLA_gene = "A", input_fp = "../IMGTHLA/")
    
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

def correct_block_typing(ID, 
                         bIndex,
                         bHLAtypings, 
                         bQuery, 
                         bPhases, 
                         bBlocks, 
                         locus):
    '''
    1. merge the same block with differnt intron lengths
    2. 
    '''
    
    #ID =ID
    #bIndex = index_ph1 
    #bHLAtypings = typing_list 
    #bQuery = query
    ##bPhases = phases
    #bBlocks = blocks
    #locus = locus
    
    serotype = [tp.split(":")[0] for tp in bHLAtypings]
    
    ARS0seq = IMGTdbIO.readIMGTsql(bHLAtypings[0], field='Exon2, Exon3')
    ARS1seq = IMGTdbIO.readIMGTsql(bHLAtypings[1], field='Exon2, Exon3')
    ps = int(bPhases[bIndex[0]])  # phase ID
    
    corrected_blocks = {}
    if locus in ['DRB1', 'DPB1'] or "DQB1*02" in serotype: # two blocks
    
        b1_index = [b_i for b_i in bIndex if bBlocks[b_i] == '1']
        if len(b1_index) == 1:
            if ARS0seq[0] in bQuery[b1_index[0]] or ARS1seq[0] in bQuery[b1_index[0]]: 
                bBlocks[b1_index[0]] = '1'
            elif ARS0seq[1] in bQuery[b1_index[0]] or ARS1seq[1] in bQuery[b1_index[0]]:
                bBlocks[b1_index[0]] = '2'
            
        else: # duplicates for each block
            
            if bQuery[b1_index[0]] in bQuery[b1_index[1]]: # the latter is the longer version
                if ARS0seq[0] in bQuery[b1_index[1]] or ARS1seq[0] in bQuery[b1_index[1]]: 
                    bBlocks[b1_index[1]] = '1'
                elif ARS0seq[1] in bQuery[b1_index[1]] or ARS1seq[1] in bQuery[b1_index[1]]:
                    bBlocks[b1_index[1]] = '2'
                bPhases.pop(b1_index[0])
                bBlocks.pop(b1_index[0])
                bQuery.pop(b1_index[0])
                b1_index.pop(0)
                
            elif bQuery[b1_index[1]] in bQuery[b1_index[0]]: # the latter is the longer version
                if ARS0seq[0] in bQuery[b1_index[0]] or ARS1seq[0] in bQuery[b1_index[0]]: 
                    bBlocks[b1_index[0]] = '1'
                elif ARS0seq[1] in bQuery[b1_index[0]] or ARS1seq[1] in bQuery[b1_index[0]]:
                    bBlocks[b1_index[0]] = '2'
                bPhases.pop(b1_index[1])
                bBlocks.pop(b1_index[1])
                bQuery.pop(b1_index[1])
                b1_index.pop(1)
        
        b2_index = [b_i for b_i in bIndex if bBlocks[b_i] == '2']
        if len(b2_index) == 1:
            if ARS0seq[0] in bQuery[b2_index[0]] or ARS1seq[0] in bQuery[b2_index[0]]: 
                bBlocks[b2_index[0]] = '1'
            elif ARS0seq[1] in bQuery[b2_index[0]] or ARS1seq[1] in bQuery[b2_index[0]]:
                bBlocks[b2_index[0]] = '2'
            
        else: # duplicates for each block
            
            if bQuery[b2_index[0]] in bQuery[b2_index[1]]: # the latter is the longer version
                if ARS0seq[0] in bQuery[b2_index[1]] or ARS1seq[0] in bQuery[b2_index[1]]: 
                    bBlocks[b2_index[1]] = '1'
                elif ARS0seq[1] in bQuery[b2_index[1]] or ARS1seq[1] in bQuery[b2_index[1]]:
                    bBlocks[b2_index[1]] = '2'
                bPhases.pop(b2_index[0])
                bBlocks.pop(b2_index[0])
                bQuery.pop(b2_index[0])
                b2_index.pop(0)
                
            elif bQuery[b2_index[1]] in bQuery[b2_index[0]]: # the latter is the longer version
                if ARS0seq[0] in bQuery[b2_index[0]] or ARS1seq[0] in bQuery[b2_index[0]]: 
                    bBlocks[b2_index[0]] = '1'
                elif ARS0seq[1] in bQuery[b2_index[0]] or ARS1seq[1] in bQuery[b2_index[0]]:
                    bBlocks[b2_index[0]] = '2'
                bPhases.pop(b2_index[1])
                bBlocks.pop(b2_index[1])
                bQuery.pop(b2_index[1])
                b2_index.pop(1)
        
        index_ps = [indexI for indexI, value in enumerate(bPhases) if value == ps]
        b1_index = [b_i for b_i in index_ps if bBlocks[b_i] == '1']
        b2_index = [b_i for b_i in index_ps if bBlocks[b_i] == '2']
        if b1_index < b2_index:
            temp_v = bQeury[b1_index]
            bQeury[b1_index] = bQeury[b2_index]
            bQeury[b2_index] = temp_v
            
            temp_v = bBlocks[b1_index]
            bBlocks[b1_index] = bBlocks[b2_index]
            bBlocks[b2_index] = temp_v
        
    else: 
        if bQuery[bIndex[0]] in bQuery[bIndex[1]]: # the latter is longer sequence
            bBlocks[bIndex[1]] = '1'
            bPhases.pop(bIndex[0])
            bBlocks.pop(bIndex[0])
            bQuery.pop(bIndex[0])
            
        elif bQuery[bIndex[1]] in bQuery[bIndex[0]]: # the latter is longer sequence
            bBlocks[bIndex[0]] = '1'
            bPhases.pop(bIndex[1])
            bBlocks.pop(bIndex[1])
            bQuery.pop(bIndex[1])
            
        else: ## if they are not contain each other then more complicated
            print("$$$$$$$$$$$$$$$$\n" + ID + " <> " + locus + " The sequences need manual check! \n" + "$$$$$$$$$$$$$$$$\n")
        
    corrected_blocks["HLATyping"] = bHLAtypings
    corrected_blocks["phases"] = bPhases
    corrected_blocks["blocks"] = bBlocks
    corrected_blocks["sequence"] = bQuery
        
    return(corrected_blocks)
    
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
        
    
def correct_phase_typing(ID, HLAtypings, query, phases, blocks, locus):
    """
    Check phase set typings; if they are homozygous, then only align to the reference sequence;
    else use IMGT/HLA database to align and get the correct typing
    """
    # test
    #ID = individual_ID
    #HLAtypings = individual_seq[locus]['GLstring']
    #query = individual_seq[locus]['Sequence']
    #phases = individual_seq[locus]['phase']
    #blocks = individual_seq[locus]['block']
    #blocks = [str(int(block)) for block in blocks]
    ## remove duplicates
    keep_items = [True] * len(query)
    for index1 in range(2):  # Phase ID
        index3 = [index2 for index2, value in enumerate(phases) if value == str(index1+1)] ## Block IDs for same phase ID
        phase_blocks = [blocks[ind] for ind in index3]
        unique_phase = list(set(phase_blocks))
        for block_ID in unique_phase:
            if phase_blocks.count(block_ID) > 1:
                dup_ind = [ind for ind, value in enumerate(phase_blocks) if value == block_ID]
                if len(dup_ind) == 2: # If there is only one duplicate
                    if query[index3[dup_ind[0]]] == query[index3[dup_ind[1]]]:
                        # then it's duplicated sequences
                        keep_items[index3[dup_ind[1]]] = False
                    else: # then need to further investigate
                        # ToDO merge blocks
                        print(ID + " <> " + locus + "<> Phase " + str(index1+1) + "<> block " + block_ID + ": The two Blocks are not identical!")
                    
                elif len(dup_ind) > 2: # more than two identical records
                    for i in range(len(dup_ind)):
                        for j in range(i + 1, len(dup_ind)):
                            if keep_items[index3[dup_ind[j]]] and query[index3[dup_ind[i]]] == query[index3[dup_ind[j]]]: 
                                keep_items[index3[dup_ind[j]]] = False
                            elif keep_items[index3[dup_ind[j]]]: # then need to further investigate
                                print(ID + " <> " +locus + "<> Phase " + str(index1+1) + "<> block " + block_ID + ": The the Blocks are not identical!")
                            #     print(ID + " <> " + locus + ": More than 2 duplicated items!")
                         # End of inner for-loop
                     # End of outer for-loop for check duplicates
                # End of If: checking identical items
           # End of If: check blocks that has the same Phase
        # End of middle for-loop: index3 - duplicated phase IDs
    # End of outer for-loop: index1 - phase IDs, 1 or 2
    new_HLAtyping = [HLAtypings[idx] for idx, v in enumerate(keep_items) if v == True]
    new_query = [query[idx] for idx, v in enumerate(keep_items) if v == True]
    new_phases = [phases[idx] for idx, v in enumerate(keep_items) if v == True]
    new_blocks = [blocks[idx] for idx, v in enumerate(keep_items) if v == True]
    
    HLAtypings = new_HLAtyping
    query = new_query
    phases = new_phases
    blocks = new_blocks
    
    HLAtypings = list(set(HLAtypings)) # unique typing list
    typing_list = sum([items.split("+") for items in HLAtypings], []) # reduce(lambda x, y: x + y, [items.split("+") for items in HLAtypings], []) 
    #this is supposed to be two typings for each locus.
 
    new_hlaTypings = typing_list
    if len(typing_list) == 1:  # Homozygous
        new_hlaTypings += new_hlaTypings 
    elif len(new_hlaTypings) > 2:  # more than 2 typings
        isAmbiguousTypes =[ "/" in item for item in new_hlaTypings]
        ambiguousIndex = [i for i, x in enumerate(isAmbiguousTypes) if x] # ambiguous type index
        unambiguousIndex = [i for i, x in enumerate(isAmbiguousTypes) if not x] # unambiguous index
        if len(ambiguousIndex) > 0:
            for ind in ambiguousIndex:
                for ind2 in unambiguousIndex:
                    if new_hlaTypings[ind2] in new_hlaTypings[ind].split("/"):
                        new_hlaTypings[ind] = new_hlaTypings[ind2]
        else: # non-ambiguous 
           print("$$$$$$$$$$$$$$$$\n " + ID + " <> " + locus + " : Multiple different typings for the same location" + "\n$$$$$$$$$$$$$$$$\n ")
   
    typing_list = new_hlaTypings
    typing_list = [re.sub("HLA-", "", tp) for tp in typing_list]
    ## remove ambiguous typings, only take the first one (more common type)
    # tplist = [re.sub("HLA-", "", tp) for tp in typing_list]
    tplist = [tp.split("/")[0] for tp in typing_list]
    #query_seqs = [Seq(q, generic_dna) for q in query]
    
    correct_HLAtypings = []
    # Quick and dirty way: only check Exons 2 and 3:
    if locus in ["A", "B", "C"]: # Class I, then only check Exons 2 and 3
    
        index_ph1 = [indexI for indexI, value in enumerate(phases) if value == str(1)]
        index_ph2 = [indexII for indexII, value in enumerate(phases) if value == str(2)]
        
        blocks = [itm + 1 for itm in range(len(index_ph1)) + range(len(index_ph2))]
            
        if rollupHLAtyping(tplist[0]) != rollupHLAtyping(tplist[1]): 
            #if the two fields are different, then need to check only ARS region
            
            ARSseq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
            
            if ARSseq[0] in query[0] and ARSseq[1] in query[0]: # it's correct phase typing
                # do nothing
                correct_HLAtypings = typing_list
            elif ARSseq[0] in query[1] and ARSseq[1] in query[1]: # swapped cases
                correct_HLAtypings = [typing_list[1], typing_list[0]]
            else: # if neither the case, then it's the wrong typing
                correct_HLAtypings = []
                print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
            
        else: # check if the two are the same typings, otherwise check other exons
    
            if tplist[0] == tplist[1]: # if the same typings, then do nothing. 
                ARSseq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
                
                if ARSseq[0] not in query[0] or ARSseq[1] not in query[0]: 
                    # if the ARS doesn't match, then it's the wrong typing
                    correct_HLAtypings = []
                    print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
                else:  # corect typing
                    correct_HLAtypings = typing_list    
                
                if len(query) == 1:
                    query += query
                    blocks += blocks
                    phases.append(u'2')
        
            else: # if they are different in the third or the fourth fields, then check other exons
                ARSseq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7')
                if ARSseq[0] in query[0] and ARSseq[1] in query[0] and ARSseq[2] in query[0] and ARSseq[3] in query[0] and ARSseq[4] in query[0] and ARSseq[5] in query[0] and ARSseq[6] in query[0]: # it's correct phase typing
                    # do nothing
                    correct_HLAtypings = typing_list
                elif ARSseq[0] in query[1] and ARSseq[1] in query[1] and ARSseq[2] in query[1] and ARSseq[3] in query[1] and ARSseq[4] in query[1] and ARSseq[5] in query[1] and ARSseq[6] in query[1]: # swapped cases
                    correct_HLAtypings = [typing_list[1], typing_list[0]]
                else: # if neither the case, then it's the wrong typing
                    correct_HLAtypings = []
                    print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
    elif locus in ["DRB1", "DPB1"]:  # two blocks
        ARSseq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
        index_ph1 = [indexI for indexI, value in enumerate(phases) if value == str(1)]
        index_ph2 = [indexII for indexII, value in enumerate(phases) if value == str(2)]
        
        blocks = [itm + 1 for itm in range(len(index_ph1)) + range(len(index_ph2))]
        
        if len(index_ph1) > 1 and len(index_ph2) > 1: # two blocks for each allele
            if ARSseq[0] in query[index_ph1[0]] or ARSseq[0] in query[index_ph1[1]]: # correct case
                correct_HLAtypings = typing_list
                if ARSseq[0] not in query[index_ph1[0]]: # block1- exon2, block2-exon3; if not, swap
                    blocks[index_ph1[0]] = '2'
                    blocks[index_ph1[1]] = '1'
            elif ARSseq[0] in query[index_ph2[0]] or ARSseq[0] in query[index_ph2[1]]:# swapped case
                correct_HLAtypings = [typing_list[1], typing_list[0]]
                if ARSseq[0] not in query[index_ph2[0]]: # block1- exon2, block2-exon3; if not, swap
                    blocks[index_ph2[0]] = '2'
                    blocks[index_ph2[1]] = '1'
            else: # if neither the case, then it's the wrong typing
                correct_HLAtypings = []
                print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
        elif (len(index_ph1) == 1 and len(index_ph2) > 1) or (len(index_ph1) > 1 and len(index_ph2) == 1): 
            print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n One phase set has only one sequence block")
            if ARSseq[0] in query[index_ph1[0]]: # correct case
                correct_HLAtypings = typing_list
            elif ARSseq[0] in query[index_ph2[0]]:# swapped case
                correct_HLAtypings = [typing_list[1], typing_list[0]]
            else: # if neither the case, then it's the wrong typing
                correct_HLAtypings = []
                print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
        elif len(index_ph1) == 1 and len(index_ph2) == 1:
            print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n Both phase sets have only one sequence block")
            if ARSseq[0] in query[index_ph1[0]]: # correct case
                correct_HLAtypings = typing_list
            elif ARSseq[0] in query[index_ph2[0]]:# swapped case
                correct_HLAtypings = [typing_list[1], typing_list[0]]
            else: # if neither the case, then it's the wrong typing
                correct_HLAtypings = []
                print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
                
    elif locus in ["DQB1"]:
        
        serotype = [tp.split(":")[0] for tp in typing_list]
        tplist_unique = list(set(typing_list))
        
        ARSseq = IMGTdbIO.readIMGTsql(typing_list[0], field='Exon2, Exon3')
        
        index_ph1 = [indexI for indexI, value in enumerate(phases) if value == str(1)]
        if len(index_ph1) > 1: # several blocks
        # merge them into one-- most likely have different length of  introns
        # ToDo: check blocks and merge or correct typing-block.
            if "DQB1*02" not in serotype[0]: 
                temp_list = typing_list
                typing_list[0] = temp_list[1]
                typing_list[1] = temp_list[0]
            
            corrected_blocks = correct_block_typing(ID, index_ph1, typing_list, query, phases, blocks, locus)
            
            typing_list = corrected_blocks["HLATyping"]
            phases = corrected_blocks["phases"] 
            blocks = corrected_blocks["blocks"] 
            query = corrected_blocks["sequence"] 
            
        index_ph2 = [indexII for indexII, value in enumerate(phases) if value == str(2)]
        if len(index_ph2) > 1:
            # merge them into one-- most likely have different length of  introns
            if "DQB1*02" not in serotype[1]: 
                temp_list = typing_list
                typing_list[0] = temp_list[1]
                typing_list[1] = temp_list[0]
            
            corrected_blocks = correct_block_typing(ID, index_ph2, typing_list, query, phases, blocks, locus)
            
            typing_list = corrected_blocks["HLATyping"]
            phases = corrected_blocks["phases"] 
            blocks = corrected_blocks["blocks"] 
            query = corrected_blocks["sequence"]
            
        index_ph1 = [indexI for indexI, value in enumerate(phases) if value == str(1)]
        index_ph2 = [indexII for indexII, value in enumerate(phases) if value == str(2)]
        
        blocks = [itm + 1 for itm in range(len(index_ph1)) + range(len(index_ph2))]

        if len(index_ph1) == 1 and len(index_ph2) == 1:
            if ARSseq[0] in query[index_ph1[0]] and ARSseq[1] in query[index_ph1[0]]: #and ARSseq[2] in query[index_ph1[0]]: # correct case
                correct_HLAtypings = typing_list
            elif ARSseq[0] in query[index_ph2[0]] and ARSseq[1] in query[index_ph2[0]]: # and ARSseq[2] in query[index_ph2[0]]:# swapped case
                correct_HLAtypings = [typing_list[1], typing_list[0]]
            else: # if neither the case, then it's the wrong typing
                correct_HLAtypings = []
                print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
        else:
            if(len(index_ph1) == 1 and len(index_ph2) > 1) or (len(index_ph1) > 1 and len(index_ph2) == 1): 
                if ARSseq[0] in query[index_ph1[0]]: # correct case
                    correct_HLAtypings = typing_list
                elif ARSseq[0] in query[index_ph2[0]]:# swapped case
                    correct_HLAtypings = [typing_list[1], typing_list[0]]
                else: # if neither the case, then it's the wrong typing
                    correct_HLAtypings = []
                    print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")
            else:
                if ARSseq[0] in query[index_ph1[0]] or ARSseq[0] in query[index_ph1[1]]: # correct case
                    correct_HLAtypings = typing_list
                elif ARSseq[0] in query[index_ph2[0]] or ARSseq[0] in query[index_ph2[1]]:# swapped case
                    correct_HLAtypings = [typing_list[1], typing_list[0]]
                else: # if neither the case, then it's the wrong typing
                    correct_HLAtypings = []
                    print("##********\n" + "ID: "+ ID + " Locus:" + locus + "\n The typings do not match the sequences. Please check the Typing!\n*******##\n")

    corrected_phaseSet = {"GLstring": correct_HLAtypings, "Sequence": query, "phase": phases, "block": blocks}
  
    return(corrected_phaseSet)



def load_seq_file(fp, BMTcaseInfo_fp, file_format = "txt"):
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
        fp = "../../rawData/xls/File_6.xls"    # for test
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
        len(seq_dict)
    BMTcaseInfo_fp = "../../rawData/SG39_caseID.csv"
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
    for individual_ID, individual_seq in new_seq_dict.iteritems():
        loci = list(individual_seq.keys())
        loci.remove("BMTcase")
        loci.remove("DRtype")
        loci.remove("Audit")
        loci.remove("Active")
        loci.remove("Comment")
        for locus in loci:

            #print(locus)
            if individual_seq[locus]['GLstring'][0] != "":
                corrected_typing = correct_phase_typing(individual_ID, individual_seq[locus]['GLstring'], individual_seq[locus]['Sequence'], 
                                                        individual_seq[locus]['phase'], individual_seq[locus]['block'], locus)
 
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
    len(corrected_seq_table)
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
    for individual_ID, individual_seq in seq_obj.iteritems():
        loci = list(individual_seq.keys())
        loci.remove("BMTcase")
        loci.remove("DRtype")
        loci.remove("Audit")
        loci.remove("Active")
        loci.remove("Comment")
        for locus in loci:
            filename = output + "SG39_HLA_" + locus + ".db"
            
            # original sequence table
            conn = sqlite3.connect(filename) # automatically creates a file if doesn't exist
            cursor = conn.cursor()
            cursor.execute('''CREATE TABLE IF NOT EXISTS OriginalSeqs
                           (BMT_caseID text, NMDP_ID text, DRtype text, 
                           Audit text, Active text, Comment text,
                           HLATyping text, PS text, Sequence text)''')
            BMT_caseID = unicode(individual_seq["BMTcase"])
            NMDP_ID = unicode(individual_ID)
            DRtype = unicode(individual_seq["DRtype"]) # isDonor(NMDP_ID)
            Audit = unicode(individual_seq["Audit"])
            Active = unicode(individual_seq["Active"])
            Comment = unicode(individual_seq["Comment"])
            
            cursor.execute('SELECT count(*) FROM OriginalSeqs WHERE NMDP_ID=?', (NMDP_ID, ))
            record_temp = cursor.fetchone() 
            if(record_temp[0] == 0): # if there is no record of this ID, then insert the record
            
                index_ph1 = [indexI for indexI, value in enumerate(individual_seq[locus]["phase"]) if value == str(1)]
                index_ph2 = [indexII for indexII, value in enumerate(individual_seq[locus]["phase"]) if value == str(2)]
                
                for PhaseID in range(2):  
                    try:
                        HLATyping = unicode(individual_seq[locus]['GLstring'][PhaseID])
                    except IndexError:
                        HLATyping = unicode("")
                        
                    PS = unicode(PhaseID + 1)
                    
                    if PS == "1":
                        if len(index_ph1) == 1:  
                            Sequence = unicode(individual_seq[locus]['Sequence'][PhaseID])
                        elif len(index_ph1) == 2:## two blocks
                            Sequence = unicode(individual_seq[locus]['Sequence'][index_ph1[0]] + "NNNNNNNNNNNNNNNNNNNN" + individual_seq[locus]['Sequence'][index_ph1[1]])
                    else:
                        if len(index_ph2) == 1:  
                            Sequence = unicode(individual_seq[locus]['Sequence'][PhaseID])
                        elif len(index_ph2) == 2:## two blocks
                            Sequence = unicode(individual_seq[locus]['Sequence'][index_ph2[0]] + "NNNNNNNNNNNNNNNNNNNN" + individual_seq[locus]['Sequence'][index_ph2[1]])
                    record = (BMT_caseID, NMDP_ID, DRtype, Audit, Active, Comment, HLATyping, PS, Sequence,)
                    
                    cursor.execute('INSERT INTO OriginalSeqs VALUES (?,?,?,?,?,?,?,?,?)', record)
                    conn.commit()
            conn.close()
            
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
                
                        ExonIntronSeqs_ps["PS" + str(ps + 1)] = findExonIntron(HLATypings = individual_seq[locus]['GLstring'][ps], 
                                                                              Sequence = [individual_seq[locus]['Sequence'][x] for x in index_ph], 
                                                                              Blocks = [individual_seq[locus]['block'][x] for x in index_ph])
                    
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
 
    
def findExonIntron(HLATypings, Sequence, Blocks):
    """
    Find Exon and Intron sequence part from query based on the reference
    One typing; 
    one or two blocks
    """
        
    IMGTglstrings = re.sub("HLA-", "", HLATypings)
    tplist = IMGTglstrings.split("/")
    #query_seqs = [Seq(q, generic_dna) for q in query]
    locus = tplist[0].split("*")[0]
    if locus in ['A', 'C']: ## # 8 exons, 7 introns, 5'-UTR, 3'-UTR
        Exonfields = ['Exon'+str(index+1) for index in range(8)]
        fieldsNum = 17 # 8+ 7 + 1 + 1
        # 0: five_prime_UTR
        # odd: Exon
        # even: Intron
        # 16: three_prime_UTR
        ExonSeq = IMGTdbIO.readIMGTsql(tplist[0], field = ", ".join(Exonfields))
        
        # query_exons = []
        # intron_seqs = []
        annotated_seq = {}
        if len(Blocks) == 1:
            query_seq = Sequence[0]
            intron_start = 0
            for Exon_index in range(8):
                if ExonSeq[Exon_index] in query_seq:
                    if Exon_index == 0: # Exon 1 
                        annotated_seq["five_prime_UTR"] = query_seq[0:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq['Exon' + str(Exon_index + 1)] = ExonSeq[Exon_index]
                        # intron_start = query_seq.index(ExonSeq[Exon_index]) + len(ExonSeq[Exon_index])
                        query_seq = re.sub(query_seq[0:query_seq.index(ExonSeq[Exon_index])+len(ExonSeq[Exon_index])], "", query_seq)
                    elif Exon_index == 7: # Exon 8
                        # annotated_seq["three_prime_UTR"] = query_seq[intron_start:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq["Intron"+ str(Exon_index)] = query_seq[0:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq['Exon' + str(Exon_index + 1)] = ExonSeq[Exon_index]
                        annotated_seq["three_prime_UTR"] = re.sub(query_seq[0:query_seq.index(ExonSeq[Exon_index])+len(ExonSeq[Exon_index])], "", query_seq)
                    else: # middle exons
                        annotated_seq["Intron"+ str(Exon_index)] = query_seq[0:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq['Exon' + str(Exon_index + 1)] = ExonSeq[Exon_index]
                        #intron_start = query_seq.index(ExonSeq[Exon_index]) + len(ExonSeq[Exon_index])
                        query_seq = re.sub(query_seq[0:query_seq.index(ExonSeq[Exon_index])+len(ExonSeq[Exon_index])], "", query_seq)
                
                else: ## ToDO: Finish ExonExtraction() function
                    print("ARS doens't match. ")
                     if Exon_index == 0: # Exon 1 
                        # annotated_seq["five_prime_UTR"] = query_seq[0:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq['Exon' + str(Exon_index + 1)] = ExonExtraction(ExonSeq[Exon_index], query_seq)
                        # intron_start = query_seq.index(ExonSeq[Exon_index]) + len(ExonSeq[Exon_index])
                        query_seq = re.sub(query_seq[0:query_seq.index(ExonSeq[Exon_index])+len(ExonSeq[Exon_index])], "", query_seq)
                    elif Exon_index == 7: # Exon 8
                        # annotated_seq["three_prime_UTR"] = query_seq[intron_start:query_seq.index(ExonSeq[Exon_index])]
                        # annotated_seq["Intron"+ str(Exon_index)] = query_seq[0:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq['Exon' + str(Exon_index + 1)] = ExonExtraction(ExonSeq[Exon_index], query_seq)
                        # annotated_seq["three_prime_UTR"] = re.sub(query_seq[0:query_seq.index(ExonSeq[Exon_index])+len(ExonSeq[Exon_index])], "", query_seq)
                    else: # middle exons
                        # annotated_seq["Intron"+ str(Exon_index)] = query_seq[0:query_seq.index(ExonSeq[Exon_index])]
                        annotated_seq['Exon' + str(Exon_index + 1)] = ExonExtraction(ExonSeq[Exon_index], query_seq)
                        # intron_start = query_seq.index(ExonSeq[Exon_index]) + len(ExonSeq[Exon_index])
                        # query_seq = re.sub(query_seq[0:query_seq.index(ExonSeq[Exon_index])+len(ExonSeq[Exon_index])], "", query_seq)

                        
                        
        else: # 2 blocks
            print("")
                
    elif locus == 'B': # 7 exons, 6 introns, 5'-UTR, 3'-UTR  
        # fields = ['Exon'+str(index+1) for index in range(7)]
        fieldsNum = 15 # 7 + 6 + 1 + 1
        # 0: five_prime_UTR
        # odd: Exon
        # even: Intron
        # 14: three_prime_UTR
    elif locus == 'DQB1': # Intron1-Exon2-Intron2-Exon3-Intron3-Exon4-Intron4
        # fields = ['Exon'+str(index+2) for index in range(3)]
        fieldsNum = 7 # 3 + 4 
        # odd: Intron
        # even: Exon
    elif locus in ['DRB1', 'DPB1']:  # Intron1-Exon2-Intron2 and intron2-Exon3-Intron3
        # fields = ['Exon'+str(index+2) for index in range(2)]
        fieldsNum = 3 # 2 + 3 
        # odd: Intron
        # even: Exon
    else:
        fields = "*"
        
    ARSseq = {}
    # Quick and dirty way: only check Exons 2 and 3:
    if field
    #if the two fields are different, then need to check only ARS region
        for Seq_index in range(len(Sequence)):
            
            ARSseq = IMGTdbIO.readIMGTsql(tplist[Seq_index], field = ", ".join(fields))
        
            query_ARS = []
            intron_seqs = []
            query_seq = Sequence[Seq_index]
            for ARS_index in range(len(fields)):
                if ARSseq[ARS_index] in query_seq:
                    intron_seqs.append(query_seq[0:query_seq.index(ARSseq[ARS_index])])
                    query_ARS.append(ARSseq[ARS_index])
                    query_seq = query_seq[query_seq.index(ARSseq[ARS_index]):len(query_seq)]
                else:
                    print("ARS doens't match. ")
            # the last intron
            if Sequence[Seq_index].index(ARSseq[ARS_index])+len(ARSseq[ARS_index]) <= len(Sequence[Seq_index]):
                 intron_seqs.append(Sequence[Seq_index][Sequence[Seq_index].index(ARSseq[ARS_index])+len(ARSseq[ARS_index]):len(Sequence[Seq_index])])
                
            for ARS_index in range(len(fields)):
                if fields[ARS_index] == "Exon1":
                    intronName = 'fivePrimeUTR'
                elif locus in ['A', 'C'] and fields[ARS_index] == "Exon8" :
                    intronName = 'fivePrimeUTR'
                Annotated_seqs['PS1'][intronName]
                        
        
        if ARSseq[0] in query[0] and ARSseq[1] in query[0]: # it's correct phase typing
            # do nothing
            correct_HLAtypings = typing_list
        elif ARSseq[0] in query[1] and ARSseq[1] in query[1]: # swapped cases
            correct_HLAtypings = [typing_list[1], typing_list[0]]
        else: # if neither the case, then it's the wrong typing
            correct_HLAtypings = []
            print("The typings do not match the sequences. Please check the Typing!\n")
            
    else: # check if the two are the same typings, otherwise check other exons
    
        if tplist[0] == tplist[1]: # if the same typings, then do nothing. 
            ARSseq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon2, Exon3')
            
            if ARSseq[0] not in query[0] or ARSseq[1] not in query[0]: 
                # if the ARS doesn't match, then it's the wrong typing
                correct_HLAtypings = []
                print("The typings do not match the sequences. Please check the Typing!\n")
            else:  # corect typing
                correct_HLAtypings = typing_list              
        
        else: # if they are different in the third or the fourth fields, then check other exons
            ARSseq = IMGTdbIO.readIMGTsql(tplist[0], field='Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7')
            if ARSseq[0] in query[0] and ARSseq[1] in query[0] and ARSseq[2] in query[0] and ARSseq[3] in query[0] and ARSseq[4] in query[0] and ARSseq[5] in query[0] and ARSseq[6] in query[0]: # it's correct phase typing
                # do nothing
                correct_HLAtypings = typing_list
            elif ARSseq[0] in query[1] and ARSseq[1] in query[1] and ARSseq[2] in query[1] and ARSseq[3] in query[1] and ARSseq[4] in query[1] and ARSseq[5] in query[1] and ARSseq[6] in query[1]: # swapped cases
                correct_HLAtypings = [typing_list[1], typing_list[0]]
            else: # if neither the case, then it's the wrong typing
                correct_HLAtypings = []
                print("The typings do not match the sequences. Please check the Typing!\n")
    
    

def align_sequence(ID, SeqObj, method = "Muscle"):
    """
    Algin the query to the reference sequences. Used to verify the typing and phase correction,
    and exon/intron parsing
    """
    ## for test
    ID = individual_ID
    SeqObj = individual_seq
    
    if method == "Muscle":
        muscle_cmd = MuscleCommandline(clwstrict = True)
       # muscle_cmd2 = MuscleCommandline(input = "../data/test_locus_A.fasta",  clw = True)
       # stdout, stderr = muscle_cmd2()
        child = subprocess.Popen(str(muscle_cmd),
                                 stdin = subprocess.PIPE,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE,
                                 universal_newlines = True,
                                 shell = (sys.platform!="win32"))
        
        #align = AlignIO.read(StringIO(stdout), "clustal")
        #print(align)
        loci = list(individual_seq.keys())
        loci.remove("BMTcase")
        loci.remove("DRtype")
        loci.remove("Audit")
        loci.remove("Active")
        loci.remove("Comment")
        
        for locus in loci:
           
            HLAtypings = list(set(SeqObj[locus]["GLstring"])) # unique typing list
            typing_list = sum([items.split("+") for items in HLAtypings], [])
            
            ## remove ambiguous typings, only take the first one (more common type)
            tplist = [re.sub("HLA-", "", tp) for tp in typing_list]
            #tplist = [tp.split("/")[0] for tp in typing_list]
            #query_seqs = [Seq(q, generic_dna) for q in query]
            
            
            
            reference = SeqRecord(seq = Seq(re.sub("-", "", IMGTdbIO.readIMGTsql(tplist[0], field='AlignedGenomSeq')[0]), IUPAC.ExtendedIUPACDNA),
                                  id = IMGTdbIO.readIMGTsql(tplist[0], field='HLATyping')[0]) 
            
            query = [SeqRecord(seq = Seq(re.sub("-", "", IMGTdbIO.readIMGTsql(tplist[0], field='AlignedGenomSeq')[0]),  IUPAC.ExtendedIUPACDNA),
                               id = str(ind)) for ind, seqV in enumerate(SeqObj[locus]['Sequence'])]
            
            SeqIO.write(reference, child.stdin, "fasta")
            SeqIO.write(query, child.stdin, "fasta")
            child.stdin.close()
            
            align = AlignIO.read(child.stdout, "clustal")
           
            #print(locus)
            if individual_seq[locus]['GLstring'][0] != "":
                
                query = 
                
                corrected_typing = correct_phase_typing(individual_ID, individual_seq[locus]['GLstring'], individual_seq[locus]['Sequence'], 
                                                        individual_seq[locus]['phase'], individual_seq[locus]['block'], locus)
 
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
    len(corrected_seq_table)
        