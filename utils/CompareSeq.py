#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 21:32:50 2017

@author: hhuang2
"""

from Bio import SeqIO #, AlignIO#, pairwise2,
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC #, generic_dna, generic_protein
#from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Align.Applications import MuscleCommandline
#from Bio.Alphabet import generic_dna
#from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os, shutil
from utils import IMGTdbIO
from collections import Counter
#import subprocess 
#import sys
#try:
#   from StringIO import StringIO # Python 2
#except ImportError:
#   from io import StringIO       # Python 3

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"


def compare_seqs(seq1, seq2, seq1_ID, seq2_ID, algn_file, saveFile = False, HLAtyping = None, DB_field ='*'):
   # muscle_cline = MuscleCommandline(clwstrict=True)
    '''
    DB_field: for Class I - 'UnalignedGenomSeq, SeqAnnotation' 
              for Class II - 'Exon2, SeqAnnotation' or 'Exon3, SeqAnnotation'
    Arguments:
        input:
            seq1 - first sequence
            seq2 - second sequence
            seq1_ID - the name of seq1
            seq2_ID - the name of seq2
            algn_file - output alignment file name, e.g. "../Output/testAlignment/TestAligned_wHLAseq.aln"
            saveFile - boolen; have to provide algn_file 
            HLAtyping - HLA typing for extracting the reference sequence
            DB_field - Reference sequence dababase column name
        
        output:
            alignment - aligned dictionary {seq1_ID: , seq2_ID:, RefSeq|HLAtyping: , AlignSymbol: }
            pos - mismatched position
            annotation - Mismatche position annotation and missing sequence part. {pos: {'PosAnn': , 'SeqFull': }}
        
        Also prints out the annotation results, e.g.
            ============
            Position 3111:
            Annotation: Intron7.69
            Sequence Completeness: Missing Partial 5'-UTR; Missing Partial 3'-UTR; 
            Donor:                  A
            Recipient:              G
            RefSeq|A*02:01:01:      G
            ============
    '''
    if HLAtyping is not None:
        
        if seq1 != seq2:
            SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
            SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')

            records = (SeqOne, SeqTwo)
            alignment, pos, CharSpaceNum_perLine = align_seqs(records, align_file = algn_file, saveFile)
            
        else:
            alignment = "Same sequence"
            pos = []
    else:
        
        if seq1 != seq2:
            # DB_field = 'UnalignedGenomSeq, SeqAnnotation'
            ReferenceSeq = IMGTdbIO.readIMGTsql(HLAtyping, field = DB_field)
            RefKey = 'RefSeq|'+HLAtyping
            RefAnnotation = ReferenceSeq[1].split(" ")
            
            SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
            SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')
            Ref = SeqRecord(Seq(ReferenceSeq[0]), id = 'RefSeq|'+HLAtyping, description = 'The second sequence')
            records = (SeqOne, SeqTwo, Ref)
            alignment, pos, CharSpaceNum_perLine = align_seqs(records, align_file , saveFile)
            
            annotation = posAnnotation(alignment, pos, RefKey, RefAnnotation)
                    
            print_MisMatchedPos(alignment, pos, CharSpaceNum_perLine, annotation)
            
        else:
            alignment = "Same sequence"
            pos = []
            annotation = {}
    
    return(alignment, pos, annotation)


def posAnnotation(alignment, pos, RefKey, RefAnnotation):
    '''
    '''
    annotation = {}
    
    AnnNum = [int(item) for item in list(set(RefAnnotation)) if item != '']
    for MMpos in pos:
        TempAnn = RefAnnotation[MMpos-alignment[RefKey][:MMpos].count('-')]
        if int(TempAnn) > 0: # Exon
            annotation[str(MMpos)] = {'PosAnn':'Exon' + TempAnn + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)}
        elif int(TempAnn) == 0: # 5'-UTR
            annotation[str(MMpos)] = {'PosAnn':'5\'-UTR' + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)}
        elif int(TempAnn) == min(AnnNum): # 3'-UTR
            annotation[str(MMpos)] = {'PosAnn':'3\'-UTR' + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)}
        else: # intron
            absTempAnn = re.sub('-', '', TempAnn)
            annotation[str(MMpos)] = {'PosAnn':'Intron' + absTempAnn + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)}
    ## Missing parts:
    Ann_Pos_frequency = Counter(RefAnnotation)
    if alignment['AlignSymbol'].find('>') != -1: # Partial length
        all_missing_pos = [ind for ind in range(len(alignment['AlignSymbol'])) if alignment['AlignSymbol'][ind] == '>']
        all_missing_ann = [RefAnnotation[ind-alignment[RefKey][:ind].count('-')] for ind in all_missing_pos]
            
        missing_ann_frequency = Counter(all_missing_ann)
        temp_missing_info = ''
        for MSpos, freq in missing_ann_frequency.items():
            if int(MSpos) > 0: # Exon
                if missing_ann_frequency[MSpos] == Ann_Pos_frequency[MSpos]: # missing full region
                    temp_missing_info = temp_missing_info + 'Missing full Exon' + MSpos + '; '
                else: # partial reigon
                    temp_missing_info = temp_missing_info + 'Missing Partial Exon' + MSpos + '; '
            elif int(MSpos) == 0: # 5'-UTR
                if missing_ann_frequency[MSpos] == Ann_Pos_frequency[MSpos]: # missing full region
                    temp_missing_info = temp_missing_info + 'Missing full 5\'-UTR' + '; '
                else: # partial reigon
                    temp_missing_info = temp_missing_info + 'Missing Partial 5\'-UTR' + '; '
                    
            elif int(MSpos) == min(AnnNum): # 3'-UTR
                if missing_ann_frequency[MSpos] == Ann_Pos_frequency[MSpos]: # missing full region
                    temp_missing_info = temp_missing_info + 'Missing full 3\'-UTR' + '; '
                else: # partial reigon
                    temp_missing_info = temp_missing_info + 'Missing Partial 3\'-UTR' + '; '
                
            else: # intron
                absMSpos = re.sub('-', '', MSpos)
                if missing_ann_frequency[MSpos] == Ann_Pos_frequency[MSpos]: # missing full region
                    temp_missing_info = temp_missing_info + 'Missing full Intron' + absMSpos + '; '
                else: # partial reigon
                    temp_missing_info = temp_missing_info + 'Missing Partial Intron'+ absMSpos + '; '
       
        annotation[str(MMpos)]['SeqFull'] = temp_missing_info
        
    else:
        annotation[str(MMpos)]['SeqFull'] = 'Completely full length sequences.'
        
    return(annotation)
        
def align_seqs(sequences, temp_fp = "../Output/Muscle_temp/", method = 'muscle', align_file = None, saveFile = False):
    '''
    '''
    if method == 'muscle':
        
        if not os.path.exists(temp_fp):
            os.makedirs(temp_fp)
        
        tempSeqFile =temp_fp+method+"_sequence_alignment_temp.fasta"
        tempAlnFile = temp_fp+method+"_aligned_sequences_temp.aln"
        
        SeqIO.write(sequences, tempSeqFile, "fasta")
        
        cline = MuscleCommandline(input=tempSeqFile, out=tempAlnFile, clwstrict=True)
        dump = cline()
        print(dump[1])
        
        reformated_algn, CharSpaceNum_perLine = reformat_clw(tempAlnFile, output = align_file, saveFile = saveFile)
        
        MM_positions = [pos for pos, char in enumerate(reformated_algn['AlignSymbol']) if char == 'X']
        
        #for ind in MM_positions:
        #    print("="*12)
        #    print("Position "+str(ind)+":")
        #    for key in reformated_algn.keys():
        #        if key != 'AlignSymbol':
        #            print(key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_algn[key][ind])
            #print("="*12)
        
    shutil.rmtree(temp_fp, ignore_errors=False, onerror=None)
    return(reformated_algn, MM_positions, CharSpaceNum_perLine)

def print_MisMatchedPos(reformated_algn, MM_positions, CharSpaceNum_perLine, Annotation = None):
    '''
    '''
    if Annotation == None:
        for ind in MM_positions:
            print("="*12)
            print("Position "+str(ind)+":")
            for key in reformated_algn.keys():
                if key != 'AlignSymbol':
                    print(key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_algn[key][ind])
            print("="*12)
    elif len(Annotation)>0: 
        for ind in MM_positions:
            print("="*12)
            print("Position "+str(ind)+":")
            print("Annotation: " + Annotation[str(ind)]['PosAnn'])
            print("Sequence Completeness: " + Annotation[str(ind)]['SeqFull'])
            for key in reformated_algn.keys():
                if key != 'AlignSymbol':
                    print(key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_algn[key][ind])
            print("="*12)       
 
       
def reformat_clw(aln_fp, output = None, saveFile = False):
    '''
    Reformat multi-lined clw alignment file into a single line format
    '''
    
    with open(aln_fp) as aln_File:
        fp_lines = aln_File.readlines()

    CharSpaceNum_perLine = fp_lines[3].rstrip().rindex(' ') + 1  # Also a sequence starting index
    #SeqcharNum_perLine = len(fp_lines[3].rstrip()) - CharSpaceNum_perLine
    
    reformated_alignment = {} 
    for line in fp_lines:
        if line != '\n' and line != 'CLUSTAL W (1.81) multiple sequence alignment\n': # actual alignment lines
            
            if line[0] == ' ': # alignment symbole line
                line_Info = line.rstrip().split(" ")
                if line_Info == ['']: # all mismatched; or this segment don't have sequence info
                    #Symbols = line.rstrip()[CharSpaceNum_perLine:]
                    #Symbols = re.sub("\*", ".", Symbols)
                    #Symbols = re.sub(" ", ">", Symbols)
                    Symbols = line[CharSpaceNum_perLine:]
                    Symbols = re.sub(" ", ">", Symbols)
                    Symbols = Symbols.rstrip()
                else: 
                    
                    Symbols = line[CharSpaceNum_perLine:]
                    Symbols = re.sub("\*", ".", Symbols)
                    Symbols = re.sub(" ", "X", Symbols)
                    Symbols = Symbols.rstrip()
                
                if 'AlignSymbol' not in reformated_alignment.keys():
                    reformated_alignment['AlignSymbol'] = Symbols
                else:    
                    reformated_alignment['AlignSymbol'] = reformated_alignment['AlignSymbol'] + Symbols
   
            else:
                Seq_ID = line[:CharSpaceNum_perLine].rstrip()
                if Seq_ID not in reformated_alignment.keys():
                    reformated_alignment[Seq_ID] = line[CharSpaceNum_perLine:].rstrip()
                else:
                    reformated_alignment[Seq_ID] = reformated_alignment[Seq_ID] + line[CharSpaceNum_perLine:].rstrip()
    
    MM_1stPos = reformated_alignment['AlignSymbol'].find('X')
    stopFlag = False
    if MM_1stPos != -1: 
        #offset = 0
        while not stopFlag:
            
            #MM_1stPos = MM_1stPos + offset
            if MM_1stPos != 0 and reformated_alignment['AlignSymbol'][MM_1stPos-1] == '>':
                trueNonExsit = False
                for seq_ID in reformated_alignment.keys():
                    if seq_ID != 'AlignSymbol' and 'RefSeq' not in seq_ID:
                        if reformated_alignment[seq_ID][MM_1stPos] == '-':
                            trueNonExsit = True
                        else:
                            trueNonExsit = False
                if trueNonExsit:
                    temp_symb = list(reformated_alignment['AlignSymbol'])
                    temp_symb[MM_1stPos] = '>'
                    reformated_alignment['AlignSymbol'] = ''.join(temp_symb)
                    if reformated_alignment['AlignSymbol'][MM_1stPos+1] == 'X':
                        MM_1stPos += 1
                    else:
                        stopFlag = True
                else:
                    stopFlag = True
            else:
                stopFlag = True
    
    
    MM_LastPos = reformated_alignment['AlignSymbol'].rfind('X')
    stopFlag = False
    if MM_LastPos != -1: 
        while not stopFlag:
            if MM_LastPos != len(reformated_alignment['AlignSymbol']) and reformated_alignment['AlignSymbol'][MM_LastPos+1] == '>':
                trueNonExsit = False
                for seq_ID in reformated_alignment.keys():
                    if seq_ID != 'AlignSymbol' and 'RefSeq' not in seq_ID:
                        if reformated_alignment[seq_ID][MM_LastPos] == '-':
                            trueNonExsit = True
                        else:
                            trueNonExsit = False
                if trueNonExsit:
                    temp_symb = list(reformated_alignment['AlignSymbol'])
                    temp_symb[MM_LastPos] = '>'
                    reformated_alignment['AlignSymbol'] = ''.join(temp_symb)
                    if reformated_alignment['AlignSymbol'][MM_LastPos-1] == 'X':
                        MM_LastPos -= 1
                    else:
                        stopFlag = True
                else:
                    stopFlag = True
            else:
                stopFlag = True
    
            
    if saveFile:
        with open(output, 'w') as f:
            for key in reformated_alignment.keys():
                if key != 'AlignSymbol':
                    f.write(key + ' '*(CharSpaceNum_perLine - len(key)) + reformated_alignment[key]+'\n')
            f.write('Align' + ' '*(CharSpaceNum_perLine - len('Align')) + reformated_alignment['AlignSymbol']+'\n')
    
    return(reformated_alignment, CharSpaceNum_perLine)
 