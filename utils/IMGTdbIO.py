#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Functions:
      1. Read IMGT/HLA aligned sequences

"""
from os import path
from collections import defaultdict
import re
import csv
import IMGTtools

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

def find_IMGT_alignment(filename, HLAtyping, header_line = 8):
    """
    Find the aligned
    """
    # Read alignment file
    alignment_list = read_IMGT_alignment(filename, header_line)
    
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

def read_IMGT_alignment(filename, seqType = 'gDNA', headOffset = 2, tailOffset = 4):
    """
    Read IMGT alignment file and construct a dictionary data structure
    headOffset: 2 for genomic and Prot alignment; 3 for CDS alignment
    tailOffset: 4 for genomic alignment; 3 for CDS and Prot alignment
    """
    #IMGT_db_fp = "../IMGTHLA/" 
    #filename = IMGT_db_fp + "/alignments/A_gen.txt"
    
    re_seqType = re.compile(r'\b'+seqType+'\\b') 
    # Read alignment file
    alignment_lines = open(filename).readlines() ### open(filename, 'r')
    
    # Line indices for the typing sequences through sequence type
    seqLineIndex = [Line + headOffset for Line, word in enumerate(alignment_lines) if re.search(re_seqType, word)]
    
    # building dictionary structure ; [key, value] = {'HLATyping': {seqType: Aligned sequeces}}
    IMGT_db = {}
    total_seq_num = seqLineIndex[1] - seqLineIndex[0] - tailOffset
    #ref_Aligned_seqs = []
    # Reference typing
    ref_HLAtyping = removeWhiteSpace(alignment_lines[seqLineIndex[0]])[0] # re.sub(" ","", alignment_lines[seqLineIndex[0]].rstrip().split("  ")[0])
    counter = 0
    for LineIndex in range(len(seqLineIndex)):
        
        # Reference sequence
        ref_Aligned_seqs = removeWhiteSpace(alignment_lines[seqLineIndex[LineIndex]])[1]# re.sub(" ", "", alignment_lines[seqLineIndex[LineIndex]].rstrip().split("  ")[2])
        
        type_index = 0
        while type_index < total_seq_num:
            
            HLAtyping, Aligned_seq_temp = removeWhiteSpace(alignment_lines[seqLineIndex[LineIndex] + type_index])#[0] # re.sub(" ","", alignment_lines[seqLineIndex[LineIndex] + type_index].rstrip().split("  ")[0])
            #try:
            #    Aligned_seq_temp = removeWhiteSpace(alignment_lines[seqLineIndex[LineIndex] + type_index])[1] # re.sub(" ", "", removeAllpattern(temp_seq, "")[])
            #except IndexError:
            #    Aligned_seq_temp = ''  # if the sequence doesn't exist in this region
                    
            # convert aligned "-" into corresponding nucleotide
            if HLAtyping != ref_HLAtyping and len(Aligned_seq_temp)>0:
                mat_index = findCharacter(Aligned_seq_temp, '-') # [ind for ind, x in enumerate(list(Aligned_seq_temp)) if x == '-']
                Aligned_seqs = list(Aligned_seq_temp)
                for x in mat_index: 
                    Aligned_seqs[x] = list(ref_Aligned_seqs)[x]
                Aligned_seqs = "".join(Aligned_seqs)
            else:
                Aligned_seqs = Aligned_seq_temp
                
            # convert "." into gap symbol "-"
            gap_index = findCharacter(Aligned_seqs, '.') # [ind for ind, x in enumerate(list(Aligned_seq_temp)) if x == '-']
            if len(gap_index)>0:
                Aligned_seqs_temp = list(Aligned_seqs)
                for x in gap_index: 
                    Aligned_seqs_temp[x] = '-'
                Aligned_seqs = "".join(Aligned_seqs_temp)
                
            # construct dict structure
            if HLAtyping != '':
                if not HLAtyping in IMGT_db: # new record
                
                    counter += 1
                    IMGT_db[HLAtyping] = {seqType:Aligned_seqs}
                    
                else: # existing record - add the sequence at the end
                    
                    IMGT_db[HLAtyping] = {seqType:IMGT_db[HLAtyping][seqType]+Aligned_seqs}
            
            type_index += 1
        #print(counter)
        
    return(IMGT_db)

def findCharacter(stringList, patternCharacter):
    """
    Find the specific character from the list and return their indices
    """
    return([ind for ind, x in enumerate(list(stringList)) if x == patternCharacter])

def removeAllpattern(stringList, patternCharacter):
    """
    Remove the specific character from the list and return
    """
    return([x for x in stringList if x != patternCharacter])

def removeWhiteSpace(stringList):
    """
    Remove the white space from the IMGT record
    """
    single_record = stringList.rstrip().split("  ")
    single_record_list = removeAllpattern(single_record, "")
    num_elem = len(single_record_list)
    lineRecord = []
    lineRecord.append(re.sub(" ","", single_record_list[0])) ## HLA typing
    if num_elem > 2:
        temp_record = []
        for ind in range(1, num_elem):
            temp_record += re.sub(" ","", single_record_list[ind])
        lineRecord.append(temp_record) ## sequence record
    else:
        try: 
            lineRecord.append(re.sub(" ", "", single_record_list[1]))
        except IndexError:
            lineRecord.append('') # if the sequence doesn't exist in this region
        
    return(lineRecord)


def parseExonSequences(seq_db, dbType = "CDS"):
    """
    Parse exon sequences from CDS or genomic sequences 
    seq_db: Dictionary structure
    dbType: "CDS" or "genomic" 
    """
    Exon_db = {}
    
    if dbType == "CDS":
        Exon_db = defaultdict(dict)
        for Typings, Seqs in seq_db.iteritems():
            boundaryIndex = findCharacter(Seqs["cDNA"], "|")
            numExons = len(boundaryIndex) + 1
            boundaryIndex.append(0)
            boundaryIndex.append(len(Seqs["cDNA"]))
            boundaryIndex.sort()
            temp_db = list()
            for index in range(numExons):
                temp_seqs = list(Seqs["cDNA"])
                Exon_seqs = []
                for x in range(boundaryIndex[index],boundaryIndex[index+1]):
                    if temp_seqs[x] != "|":
                        Exon_seqs.append(temp_seqs[x])
                temp_db.append({Typings:{'Exon'+str(index+1):"".join(Exon_seqs)}})
            for item in temp_db:
                for k, v in item.iteritems():
                    Exon_db[k].update(v)
    return(Exon_db)

def IMBTdb_2_dict(HLA_gene = "A", input_fp = "../IMGTHLA/", outfile):
    """
    Convert IMGT database into dictionary structure.
    """
    # genomic alignment file
    filename = input_fp + "/alignments/" + HLA_gene + "_gen.txt" 
    if path.exists(filename):
        gDNA_alignment = read_IMGT_alignment(filename, 'gDNA')
#    else:
#        gDNA_alignment = {}
    
    # CDS alignemtn file
    filename = input_fp + "/alignments/" + HLA_gene + "_nuc.txt"
    if path.exists(filename):    
        CDS_alignment = read_IMGT_alignment(filename, 'cDNA', 3, 4)
#    else:
#        CDS_alignment = {}
    if len(CDS_alignment)>0:
        ExonSequences = parseExonSequences(CDS_alignment)
    
    # protein alignment file
    filename = input_fp + "/alignments/" + HLA_gene + "_prot.txt"
    if path.exists(filename): 
        protein_alignment = read_IMGT_alignment(filename, 'Prot', 2, 3)
#    else:
#        protein_alignment = {}
    
 # merge into one single Dictionary structure
    combined_alignments = [gDNA_alignment, CDS_alignment, ExonSequences, protein_alignment]
    
    combined_dict = defaultdict(dict)
    for item in combined_alignments:
        for k, v in item.iteritems():
            combined_dict[k].update(v)
            
    return(combined_dict)