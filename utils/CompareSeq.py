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

def compare_DQB1_Targeted_Region(Sequences, params):
    '''
    Targeted region alignment and comparison for DQB1.
    DQB1*02 have two segments --  Exon2 and Exon3 
    the rest DQB1 have one segment inlcuding -- partial intron1 to Exon4
    -------------
    Arguments:
        input:
            Sequence: sequence query, Dict, {'seq1_ID': seq, 'seq2_ID': seq2, ....}
            params: {'algn_file': algn_file, 'saveFile': saveFile, 'HLAtyping': HLAtyping}
        
        output:
            Alignment - aligned dictionary {'Exon2': Exon2_alignment:, 'Exon3': Exon3_alignment}
            Annotation - Mismatche position annotation and missing sequence part. {'Exon2': Exon2_alignment:, 'Exon3': Exon3_alignment}
        
        Also prints out the annotation results, e.g.
            Sequence Completeness: Exon2, Exon3
            ============
            Position 314:
            Annotation: Exon2.75
            Recipient-PS2|Exon2:                C
            RefSeq|DRB1*15:03:01:01|Exon2:      C
            RefSeq|DRB1*15:03:01:02|Exon2:      C
            Donor-PS2|Exon2:                    T
            RefSeq|DRB1*15:01:01:01|Exon2:      T
            RefSeq|DRB1*15:01:01:02|Exon2:      T
            RefSeq|DRB1*15:01:01:03|Exon2:      T
            RefSeq|DRB1*15:01:01:04|Exon2:      T
            ============
    '''
    try:
        seqCount = len(Sequences)
        if seqCount == 2: 
            keys = list(Sequences.keys())
            seq1_ID = keys[0]
            seq1 = Sequences[seq1_ID]
            seq2_ID = keys[1]
            seq2 = Sequences[seq2_ID]
            
    except TypeError as Te:
        print(Te)
        
    try:
        algn_file = params['algn_file']
        saveFile = params['saveFile']
        HLAtyping = params['HLAtyping']
        SeroType = HLAtyping[0].split(':')[0]
    except TypeError as e:
        print(e)
        print('Please provide the following paramters into one Dictionary: algn_file, saveFile, HLAtyping, DB_field')
    
    #locus = HLAtyping[0].split('*')[0]
    # Exons = DB_field.split(',')
    Alignment = {}
    Annotation = {}
    if seqCount == 2:
        ReferenceSeq = {}
        AddtnlMessage = ''
        
        if SeroType == 'DQB1*02': # DQB1*02: has two segments
            
            seq1_Exons = seq1.split('-'*10)
            seq2_Exons = seq2.split('-'*10)
            
            if len(seq1_Exons) == 1:
                AddtnlMessage = AddtnlMessage + seq1_ID + ' is missing the Exon sequence(s).'
                seq1_Exons.append(seq1_Exons[0])
            
            if len(seq2_Exons) == 1:
                AddtnlMessage = AddtnlMessage + seq2_ID + ' is missing the Exon sequence(s).'
                seq2_Exons.append(seq2_Exons[0])
                
           # if len(seq1_Exons) == 2: # Exon2 and Exon3
            Exons = 'Exon2, Exon3'
            for typing in HLAtyping:
                
                RefKey = 'RefSeq|'+ typing
                tempRefseq = IMGTdbIO.readIMGTsql(typing, field = Exons)
                
                if tempRefseq[0] != '': # the original Field sequences exist
                    ## Exon 2
                    ReferenceSeq[RefKey+'|Exon2'] = tempRefseq[0] #[re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                    ReferenceSeq[RefKey+'|Exon3'] = tempRefseq[1]
    
                else:
                    AddtnlMessage = AddtnlMessage + "\n" + RefKey + ' doesn\'t have a confirmed Exon 2 sequence.'
    
            for exon_id in range(2):
                # Exons 2 and 3
                SeqOne = SeqRecord(Seq(seq1_Exons[exon_id]), id = seq1_ID+'|Exon'+str(exon_id+2), description = 'The first sequence Exon '+str(exon_id+2))
                SeqTwo = SeqRecord(Seq(seq2_Exons[exon_id]), id = seq2_ID+'|Exon'+str(exon_id+2), description = 'The second sequence Exon '+str(exon_id+2))
                
                records = (SeqOne, SeqTwo)
                #Ref = SeqRecord(Seq(ReferenceSeq[0]), id = 'RefSeq|'+HLAtyping, description = 'The second sequence')
                
                for itemID, itemSeq in ReferenceSeq.items():
                    if 'Exon'+str(exon_id+2) in itemID:
                        records += (SeqRecord(Seq(itemSeq), id = itemID, description = 'Reference IMGT/HLA sequence.'), )
                
                alignment, _, CharSpaceNum_perLine = align_seqs(records)
            
                alignment, PosAnnotation = correct_align_symbols(alignment, 'Exon'+str(exon_id+2))
                
                #annotation = Targeted_posAnnotation(alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount)
                PosAnnotation['SeqFull'] = Exons
                PosAnnotation['Note'] = AddtnlMessage
                PosAnnotation['SameSeqs'] = False
                if len([pos for pos, char in enumerate(alignment['AlignSymbol']) if char == 'X']) == 0:
                    PosAnnotation['Note'] += ' The sequences are seemingly the same at the Exon region.'
                    PosAnnotation['SameSeqs'] = True
                
                ## Check Intron region
                AlignSymbol = list(alignment['AlignSymbol'])
                #AnnNum = [int(item) for item in list(set(RefAnnotation)) if item != '']

                nonRefKeys = list(alignment.keys())
                for key  in alignment.keys():
                    if 'Ref'  in key or 'AlignSymbol'== key:
                        nonRefKeys.remove(key)

                # AlignedSymbol = list(alignment['AlignSymbol'])
                Intron_startID = 0
                Intron_endID = alignment['AlignSymbol'].find('.') - 1 # re.search('\\.+', alignment['AlignSymbol']).start()-1
                
                for ind in range(Intron_endID, Intron_startID-1, -1):
                    if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
                        
                        PosAnnotation[str(ind)] = 'Intron'+str(exon_id+1) + '.-' + str(Intron_endID-ind+1)
                        AlignSymbol[ind] = 'X'
                
                Intron_startID = alignment['AlignSymbol'].rfind('.') # re.search('\\.+', alignment['AlignSymbol']).end()
             
                for ind in range(len(alignment['AlignSymbol'])-1, Intron_startID-1, -1):
                    if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
            
                        PosAnnotation[str(ind)] = 'Intron'+str(exon_id+2) + '.' + str(Intron_endID-ind+1)
                        AlignSymbol[ind] = 'X'
                
                alignment['AlignSymbol'] = ''.join(AlignSymbol)
                
                
                if saveFile:
                    save_aln(alignment, CharSpaceNum_perLine, algn_file, PosAnnotation)
                
                print_MisMatchedPos(alignment, CharSpaceNum_perLine, PosAnnotation)
                
                Alignment['Exon'+str(exon_id+2)] = alignment
                Annotation['Exon'+str(exon_id+2)] = PosAnnotation
        
        else: # non-DQB1*03
            
            DBfields = 'UnalignedGenomSeq, SeqAnnotation'
            ReferenceSeq = {}
            potentialKeys = []
            for typing in HLAtyping:
                #print(typing)
                RefKey = 'RefSeq|'+ typing
                tempRefseq = IMGTdbIO.readIMGTsql(typing, field = DBfields)
                if tempRefseq[0] != '': # the original Field sequences exist
                    ReferenceSeq[RefKey] = [re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                    potentialKeys.append(RefKey)
            
            if len(potentialKeys) > 0: # if the genomic sequence exist, then align the sequences
                
                RefKey = potentialKeys[0]
                RefAnnotation = ReferenceSeq[RefKey][1].split(" ")
                Ref_frontMisCount = list(ReferenceSeq[RefKey][0])[:505].count('*')
                
                SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
                SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')
                
                records = (SeqOne, SeqTwo)
                #Ref = SeqRecord(Seq(ReferenceSeq[0]), id = 'RefSeq|'+HLAtyping, description = 'The second sequence')
                
                for itemID, itemSeq in ReferenceSeq.items():
                    records += (SeqRecord(Seq(itemSeq[0]), id = itemID, description = 'Reference IMGT/HLA sequence.'), )
                
                Alignment, pos, CharSpaceNum_perLine = align_seqs(records)
                
                Annotation = posAnnotation(Alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount)
                Annotation['Note'] = AddtnlMessage
                
                ## Check Intron region
                AlignSymbol = list(Alignment['AlignSymbol'])
                AnnNum = [int(item) for item in list(set(RefAnnotation)) if item != '']

                nonRefKeys = list(Alignment.keys())
                for key  in Alignment.keys():
                    if 'RefSeq'  in key or 'AlignSymbol'== key:
                        nonRefKeys.remove(key)
                
                # AlignedSymbol = list(alignment['AlignSymbol'])
                
                # Front intron regions
                Intron_startID = 0
                Intron_endID = Alignment['AlignSymbol'].find('.') - 1# re.search('\\.+', alignment['AlignSymbol']).start()-1
                RefSeq_alignment = Alignment[RefKey]
                for ind in range(Intron_endID, Intron_startID-1, -1):
                    if Alignment[nonRefKeys[0]][ind] != Alignment[nonRefKeys[1]][ind]:
                        RegionSymb = RefAnnotation[ind - RefSeq_alignment[:ind].count('-')]
                        if RegionSymb == '0':
                            RegionName = '5\'-UTR'
                        elif int(RegionSymb) == min(AnnNum):
                            RegionName = '3\'-UTR'
                        elif int(RegionSymb) > 0:
                            RegionName = 'Exon'+ RegionSymb
                        elif int(RegionSymb) < 0:
                            RegionName = 'Intron'+ str(abs(int(RegionSymb)))
                        Annotation[str(ind)] = RegionName + '.-' + str(Intron_endID-ind+1)
                        AlignSymbol[ind] = 'X'
                
                #alignment['AlignSymbol'] = ''.join(AlignSymbol)
                # Back intron regions
                Intron_startID = Alignment['AlignSymbol'].rfind('.') #re.search('\\.+', alignment['AlignSymbol']).end()
             
                for ind in range(len(Alignment['AlignSymbol'])-1, Intron_startID-1, -1):
                    if Alignment[nonRefKeys[0]][ind] != Alignment[nonRefKeys[1]][ind]:
                        RegionSymb = RefAnnotation[ind - RefSeq_alignment[:ind].count('-')]
                        if RegionSymb == '0':
                            RegionName = '5\'-UTR'
                        elif int(RegionSymb) == min(AnnNum):
                            RegionName = '3\'-UTR'
                        elif int(RegionSymb) > 0:
                            RegionName = 'Exon'+ RegionSymb
                        elif int(RegionSymb) < 0:
                            RegionName = 'Intron'+ str(abs(int(RegionSymb)))
                            
                        Annotation[str(ind)] = RegionName + '.' + str(Intron_endID-ind+1)
                        AlignSymbol[ind] = 'X'
                
                Alignment['AlignSymbol'] = ''.join(AlignSymbol)
            
                if saveFile:
                    save_aln(Alignment, CharSpaceNum_perLine, algn_file, Annotation)
                
                print_MisMatchedPos(Alignment, CharSpaceNum_perLine, Annotation)
                
            else: # if there is no genomic sequence
                
                DBfields = 'Exon2, Exon3, Exon4'
                ReferenceSeq = {}
                #potentialKeys = []
                for typing in HLAtyping:
                    #print(typing)
                    RefKey = 'RefSeq|'+ typing
                    tempRefseq = IMGTdbIO.readIMGTsql(typing, field = DBfields)

                    if tempRefseq[0] != '': # the original Field sequences exist
                        ## Exon 2
                        ReferenceSeq[RefKey+'|Exon2'] = tempRefseq[0] #[re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                        ReferenceSeq[RefKey+'|Exon3'] = tempRefseq[1]
                        ReferenceSeq[RefKey+'|Exon4'] = tempRefseq[1]

                    else:
                        AddtnlMessage = AddtnlMessage + "\n" + RefKey + ' doesn\'t have a confirmed Exon 2 sequence.'
                
                SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
                SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')
                
                # records = (SeqOne, SeqTwo)
                Intron_startID = 0 
                for exon_id in range(3):
                    # Exons 2 and 3 and 4
                    records = (SeqOne, SeqTwo)
                    for itemID, itemSeq in ReferenceSeq.items():
                        if 'Exon'+str(exon_id+2) in itemID:
                            records += (SeqRecord(Seq(itemSeq), id = itemID, description = 'Reference IMGT/HLA sequence.'), )
                    
                    alignment, _, CharSpaceNum_perLine = align_seqs(records)
                
                    alignment, PosAnnotation = correct_align_symbols(alignment, 'Exon'+str(exon_id+2))
                    
                    #annotation = Targeted_posAnnotation(alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount)
                    PosAnnotation['SeqFull'] = DBfields
                    PosAnnotation['Note'] = AddtnlMessage
                    PosAnnotation['SameSeqs'] = False
                    
                    ## Check Intron region
                    AlignSymbol = list(alignment['AlignSymbol'])
                    if exon_id < 2: # check front intron
                        nonRefKeys = list(alignment.keys())
                        for key  in alignment.keys():
                            if 'Ref'  in key or 'AlignSymbol'== key:
                                nonRefKeys.remove(key)
                        
                        
                        # AlignedSymbol = list(alignment['AlignSymbol'])
                        Intron_endID = alignment['AlignSymbol'].rfind('.') #re.search('\\.+', alignment['AlignSymbol']).start()-1
                        
                        for ind in range(Intron_endID, Intron_startID-1, -1):
                            if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
                                PosAnnotation[str(ind)] = 'Intron' + str(exon_id+1) + '.-' + str(Intron_endID-ind+1)
                                AlignSymbol[ind] = 'X'
                        
                        alignment['AlignSymbol'] = ''.join(AlignSymbol)
                        Intron_startID = alignment['AlignSymbol'].rfind('.')#re.search('\\.+', alignment['AlignSymbol']).end()
                    
                    else: # the last exon front and back introns
                        # AlignedSymbol = list(alignment['AlignSymbol'])
                        Intron_endID = alignment['AlignSymbol'].rfind('.')-1 #re.search('\\.+', alignment['AlignSymbol']).start()-1
                        
                        for ind in range(Intron_endID, Intron_startID-1, -1):
                            if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
                                PosAnnotation[str(ind)] = 'Intron' + str(exon_id+1) + '.-' + str(Intron_endID-ind+1)
                                AlignSymbol[ind] = 'X'
                        
                        #alignment['AlignSymbol'] = ''.join(AlignSymbol)
                        Intron_startID = alignment['AlignSymbol'].rfind('.') # re.search('\\.+', alignment['AlignSymbol']).end()
                        
                        for ind in range(len(alignment['AlignSymbol'])-1, Intron_startID-1, -1):
                            if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
                                PosAnnotation[str(ind)] = 'Intron' + str(exon_id+2) + '.' + str(Intron_endID-ind+1)
                                AlignSymbol[ind] = 'X'
                        
                        alignment['AlignSymbol'] = ''.join(AlignSymbol)
                        
                    ## Note and SamSeqs
                    if len([pos for pos, char in enumerate(alignment['AlignSymbol']) if char == 'X']) == 0:
                        PosAnnotation['Note'] += ' The sequences are seemingly the same at Exon'+ str(exon_id+2) +'.'
                        PosAnnotation['SameSeqs'] = True
                    
                    if saveFile:
                        save_aln(alignment, CharSpaceNum_perLine, algn_file, PosAnnotation)
                    
                    print_MisMatchedPos(alignment, CharSpaceNum_perLine, PosAnnotation)
                    
                    Alignment['Exon'+str(exon_id+2)] = alignment
                    Annotation['Exon'+str(exon_id+2)] = PosAnnotation    
                
    return(Alignment, Annotation)

def compare_Targeted_Region(Sequences, params):
    '''
    targeted alignment and comparison. Especially for Class II; comparing specific Exons 
    when there is no genomic reference sequence available.
    -------------
    Arguments:
        input:
            Sequence: sequence query, Dict, {'seq1_ID': seq, 'seq2_ID': seq2, ....}
            params: {'algn_file': algn_file, 'saveFile': saveFile, 'HLAtyping': HLAtyping}
        
        output:
            Alignment - aligned dictionary {'Exon2': Exon2_alignment:, 'Exon3': Exon3_alignment}
            Annotation - Mismatche position annotation and missing sequence part. {'Exon2': Exon2_alignment:, 'Exon3': Exon3_alignment}
        
        Also prints out the annotation results, e.g.
            Sequence Completeness: Exon2, Exon3
            ============
            Position 314:
            Annotation: Exon2.75
            Recipient-PS2|Exon2:                C
            RefSeq|DRB1*15:03:01:01|Exon2:      C
            RefSeq|DRB1*15:03:01:02|Exon2:      C
            Donor-PS2|Exon2:                    T
            RefSeq|DRB1*15:01:01:01|Exon2:      T
            RefSeq|DRB1*15:01:01:02|Exon2:      T
            RefSeq|DRB1*15:01:01:03|Exon2:      T
            RefSeq|DRB1*15:01:01:04|Exon2:      T
            ============
    '''
    try:
        seqCount = len(Sequences)
        if seqCount == 2: 
            keys = list(Sequences.keys())
            seq1_ID = keys[0]
            seq1 = Sequences[seq1_ID]
            seq2_ID = keys[1]
            seq2 = Sequences[seq2_ID]
            
    except TypeError as Te:
        print(Te)
        
    try:
        algn_file = params['algn_file']
        saveFile = params['saveFile']
        HLAtyping = params['HLAtyping']
        
    except TypeError as e:
        print(e)
        print('Please provide the following paramters into one Dictionary: algn_file, saveFile, HLAtyping, DB_field')
    
    #locus = HLAtyping[0].split('*')[0]
    # Exons = DB_field.split(',')
    Alignment = {}
    Annotation = {}
    if seqCount == 2:
        ReferenceSeq = {}
        AddtnlMessage = ''
        
        seq1_Exons = seq1.split('-'*10)
        seq2_Exons = seq2.split('-'*10)
        
        if len(seq1_Exons) == 1:
            AddtnlMessage = AddtnlMessage + seq1_ID + ' is missing the Exon sequence(s).'
            seq1_Exons.append(seq1_Exons[0])
        
        if len(seq2_Exons) == 1:
            AddtnlMessage = AddtnlMessage + seq2_ID + ' is missing the Exon sequence(s).'
            seq2_Exons.append(seq2_Exons[0])
            
       # if len(seq1_Exons) == 2: # Exon2 and Exon3
        Exons = 'Exon2, Exon3'
        for typing in HLAtyping:
            
            RefKey = 'RefSeq|'+ typing
            tempRefseq = IMGTdbIO.readIMGTsql(typing, field = Exons)
            
            if tempRefseq[0] != '': # the original Field sequences exist
                ## Exon 2
                ReferenceSeq[RefKey+'|Exon2'] = tempRefseq[0] #[re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                ReferenceSeq[RefKey+'|Exon3'] = tempRefseq[1]

            else:
                AddtnlMessage = AddtnlMessage + "\n" + RefKey + ' doesn\'t have a confirmed Exon 2 sequence.'

        for exon_id in range(2):
            # Exons 2 and 3
            SeqOne = SeqRecord(Seq(seq1_Exons[exon_id]), id = seq1_ID+'|Exon'+str(exon_id+2), description = 'The first sequence Exon '+str(exon_id+2))
            SeqTwo = SeqRecord(Seq(seq2_Exons[exon_id]), id = seq2_ID+'|Exon'+str(exon_id+2), description = 'The second sequence Exon '+str(exon_id+2))
            
            records = (SeqOne, SeqTwo)
            #Ref = SeqRecord(Seq(ReferenceSeq[0]), id = 'RefSeq|'+HLAtyping, description = 'The second sequence')
            
            for itemID, itemSeq in ReferenceSeq.items():
                if 'Exon'+str(exon_id+2) in itemID:
                    records += (SeqRecord(Seq(itemSeq), id = itemID, description = 'Reference IMGT/HLA sequence.'), )
            
            alignment, _, CharSpaceNum_perLine = align_seqs(records)
        
            alignment, PosAnnotation = correct_align_symbols(alignment, 'Exon'+str(exon_id+2))
            
            #annotation = Targeted_posAnnotation(alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount)
            PosAnnotation['SeqFull'] = Exons
            PosAnnotation['Note'] = AddtnlMessage
            PosAnnotation['SameSeqs'] = False
            if len([pos for pos, char in enumerate(alignment['AlignSymbol']) if char == 'X']) == 0:
                PosAnnotation['Note'] += ' The sequences are seemingly the same at the Exon region.'
                PosAnnotation['SameSeqs'] = True
            
            ## Check Intron region
            AlignSymbol = list(alignment['AlignSymbol'])
            #AnnNum = [int(item) for item in list(set(RefAnnotation)) if item != '']

            nonRefKeys = list(alignment.keys())
            for key  in alignment.keys():
                if 'Ref'  in key or 'AlignSymbol'== key:
                    nonRefKeys.remove(key)

            # AlignedSymbol = list(alignment['AlignSymbol'])
            Intron_startID = 0
            Intron_endID = alignment['AlignSymbol'].find('.') - 1 # re.search('\\.+', alignment['AlignSymbol']).start()-1
            
            for ind in range(Intron_endID, Intron_startID-1, -1):
                if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
                    
                    PosAnnotation[str(ind)] = 'Intron'+str(exon_id+1) + '.-' + str(Intron_endID-ind+1)
                    AlignSymbol[ind] = 'X'
            
            Intron_startID = alignment['AlignSymbol'].rfind('.') # re.search('\\.+', alignment['AlignSymbol']).end()
         
            for ind in range(len(alignment['AlignSymbol'])-1, Intron_startID-1, -1):
                if alignment[nonRefKeys[0]][ind] != alignment[nonRefKeys[1]][ind]:
        
                    PosAnnotation[str(ind)] = 'Intron'+str(exon_id+2) + '.' + str(Intron_endID-ind+1)
                    AlignSymbol[ind] = 'X'
            
            alignment['AlignSymbol'] = ''.join(AlignSymbol)
            
            
            
            if saveFile:
                save_aln(alignment, CharSpaceNum_perLine, algn_file, PosAnnotation)
            
            print_MisMatchedPos(alignment, CharSpaceNum_perLine, PosAnnotation)
            
            Alignment['Exon'+str(exon_id+2)] = alignment
            Annotation['Exon'+str(exon_id+2)] = PosAnnotation

            
    return(Alignment, Annotation)

def correct_align_symbols(alignment, ExonID):
    '''
    Correct the alignment symbols for class II. 
    True mismatches and missing sequences in reference sequence
    '''
    
    Ref_startID = float('inf') 
    Ref_endID = 0
    for itemID, itemSeq in alignment.items():
        if 'RefSeq' in itemID:
            tempStart = re.search('[ATCG]+', itemSeq).start()
            tempEnd = re.search('[ATCG]+', itemSeq).end() - 1
            if tempStart < Ref_startID:
                Ref_startID = tempStart
            if tempEnd > Ref_endID:
                Ref_endID = tempEnd
    
    AlignedSymbol = list(alignment['AlignSymbol'])
    PosAnnotation = {}
    for ind in range(len(AlignedSymbol)):
        if ind < Ref_startID:
            AlignedSymbol[ind] = '>'
        elif ind > Ref_endID:
            AlignedSymbol[ind] = '>'
        elif AlignedSymbol[ind] == 'X':
            PosAnnotation[str(ind)] = ExonID+'.' + str(ind - Ref_startID + 1)
            
    corrected_alignment = alignment
    corrected_alignment['AlignSymbol'] = ''.join(AlignedSymbol)
    
    return(corrected_alignment, PosAnnotation)

def compare_seqs(Sequence, params):#algn_file, saveFile = False, HLAtyping = None, DB_field ='*'):
   # muscle_cline = MuscleCommandline(clwstrict=True)
    '''
    DB_field: for Class I - 'UnalignedGenomSeq, SeqAnnotation' 
              
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
        
        input (updated):
            Sequence: sequence query, Dict, {'seq1_ID': seq, 'seq2_ID': seq2, ....}
            params: {'algn_file': algn_file, 'saveFile': saveFile, 'HLAtyping': HLAtyping, 'DB_field': DB_field}
        
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
    try:
        seqCount = len(Sequence)
        if seqCount == 2: 
            keys = list(Sequence.keys())
            seq1_ID = keys[0]
            seq1 = Sequence[seq1_ID]
            seq2_ID = keys[1]
            seq2 = Sequence[seq2_ID]
            
    except TypeError as Te:
        print(Te)
        
    try:
        algn_file = params['algn_file']
        saveFile = params['saveFile']
        HLAtyping = params['HLAtyping']
        DB_field = params['DB_field']
        DB_BackUpfield = params['DB_BackUpfield']
    except TypeError as e:
        print(e)
        print('Please provide the following paramters into one Dictionary: algn_file, saveFile, HLAtyping, DB_field')
    
    if HLAtyping is None:
        if seqCount == 2:
            annotation = {}
            if seq1 != seq2:
                SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
                SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')
    
                records = (SeqOne, SeqTwo)
                alignment, pos, CharSpaceNum_perLine = align_seqs(records, algn_file, saveFile)
                
                if saveFile:
                    save_aln(alignment, 20, algn_file)
            else:
                alignment = "Same sequence"
                pos = []
                
                print('The two sequences are exactly the same!')
                
        else: 
            records = ()
            counter = 0
            for key, seq in Sequence.items(): 
                records += (SeqRecord(Seq(seq), id = key, description = 'Sequence order: '+str(counter)), )
                counter += 1
            
            alignment, pos, CharSpaceNum_perLine = align_seqs(records)
            
            if saveFile:
                save_aln(alignment, 20, algn_file)
    else:
        locus = HLAtyping[0].split('*')[0]
        if seqCount == 2:  # two sequence +/- HLA ref alignment and annotation.
            if seq1 != seq2:
                # DB_field = 'UnalignedGenomSeq, SeqAnnotation'
                ReferenceSeq = {}
                
                AddtnlMessage = ''
                #potentialKeys = HLAtyping
                potentialKeys = []
                for typing in HLAtyping:
                    #print(typing)
                    RefKey = 'RefSeq|'+ typing
                    tempRefseq = IMGTdbIO.readIMGTsql(typing, field = DB_field)
                    if tempRefseq[0] != '': # the original Field sequences exist
                        ReferenceSeq[RefKey] = [re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                        potentialKeys.append(RefKey)
                        backUpField_flag = False
                        noRef = False
                    elif DB_BackUpfield != '': # If the original Field sequences don't exist, then compare the backup fields
                        # : Exon alignment
                        #tempRefseq = IMGTdbIO.readIMGTsql(typing, field = DB_BackUpfield)
                        #if tempRefseq[0] != '': # the original Field sequences exist
                        #    ReferenceSeq[RefKey] = [re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                        #    potentialKeys.append(RefKey)
                        backUpField_flag = True
                        AddtnlMessage = AddtnlMessage + '\nNo genomic refernce is available. Exon regions are assessed instead.'
                        noRef = False
                    else:
                        AddtnlMessage = AddtnlMessage + "\n" + RefKey + ' doesn\'t have a confirmed genomic sequence.'
                        backUpField_flag = False
                        noRef = True
                        #potentialKeys.remove(typing)
                        
                if len(potentialKeys) > 0: ## First choice field ref sequences exist
                    RefKey = potentialKeys[0]
                    RefAnnotation = ReferenceSeq[RefKey][1].split(" ")
                    Ref_frontMisCount = list(ReferenceSeq[RefKey][0])[:505].count('*')
                    
                    SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
                    SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')
                    
                    records = (SeqOne, SeqTwo)
                    #Ref = SeqRecord(Seq(ReferenceSeq[0]), id = 'RefSeq|'+HLAtyping, description = 'The second sequence')
                    
                    for itemID, itemSeq in ReferenceSeq.items():
                        records += (SeqRecord(Seq(itemSeq[0]), id = itemID, description = 'Reference IMGT/HLA sequence.'), )
                    
                    alignment, pos, CharSpaceNum_perLine = align_seqs(records)
                    
                    annotation = posAnnotation(alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount)
                    annotation['Note'] = AddtnlMessage
                    if saveFile:
                        save_aln(alignment, CharSpaceNum_perLine, algn_file, annotation)
                    
                    print_MisMatchedPos(alignment, CharSpaceNum_perLine, annotation)

                elif backUpField_flag: # back up filed provided in case the first choice field doesn't have ref sequences.

                    records = {seq1_ID: seq1, seq2_ID: seq2}
                    alignment = {} 
                    
                    #if locus in ['A', 'B', 'C']:
                    QueryAnnotation = SeqAnnotation(records, locus, Ref = HLAtyping[0])
                    
                    annotation = QueryAnnotation['annotation']
                    CharSpaceNum_perLine = QueryAnnotation['CharSpaceNum_perLine']
                    pos = QueryAnnotation['pos']
                    alignment = QueryAnnotation['alignment']
                    
                    #if saveFile:
                    #    save_aln(alignment, CharSpaceNum_perLine, algn_file, annotation)
                
                    #print_MisMatchedPos(alignment, CharSpaceNum_perLine, annotation)
                
                    Exons = DB_BackUpfield.split(', ')
                    
                    ExonSeqs = {}
                    
                    for exonID in range(len(Exons)):
                        for typing in HLAtyping:
                            RefKey = 'RefSeq|'+ typing
                            tempRefseq = IMGTdbIO.readIMGTsql(typing, field = DB_BackUpfield)
                            if tempRefseq[exonID] != '': # the original Field sequences exist
                                if '*' not in tempRefseq[exonID]:
                                    if Exons[exonID] in ExonSeqs.keys():
                                        ExonSeqs[Exons[exonID]][RefKey+'|'+Exons[exonID]] = tempRefseq[exonID]
                                        
                                        ExonNumber = re.search('[0-9]+', Exons[exonID]).group()
                                        QueryExonIndex = [ind for ind, Symbol in enumerate(QueryAnnotation['annotation']['PosAnnotation']) if Symbol == ExonNumber]
                                        
                                        for key in records.keys():
                                            tempReads = [QueryAnnotation['alignment'][key][ind] for ind in QueryExonIndex]
                                            ExonSeqs[Exons[exonID]][key+'|'+Exons[exonID]] = ''.join(tempReads)
                                        
                                    else:
                                        ExonSeqs[Exons[exonID]] = {RefKey+'|'+Exons[exonID]: tempRefseq[exonID]}
                                        
                                        ExonNumber = re.search('[0-9]+', Exons[exonID]).group()
                                        QueryExonIndex = [ind for ind, Symbol in enumerate(QueryAnnotation['annotation']['PosAnnotation']) if Symbol == ExonNumber]
                                        
                                        for key in records.keys():
                                            tempReads = [QueryAnnotation['alignment'][key][ind] for ind in QueryExonIndex]
                                            ExonSeqs[Exons[exonID]][key+'|'+Exons[exonID]] = ''.join(tempReads)
                    
                    # Align Exon region and save the annotation results                        
                    ExonAnnotation = AlignExons(ExonSeqs, locus, algn_file, saveFile)
                    annotation['Note'] = AddtnlMessage + 'No genomic Reference sequence available. Checking ARS...'
                    #try:
                    #    len(alignment)
                    #except NameError:
                    #    alignment = ['Check Exon Annotation']
                    
                    #try:
                    #    len(pos)
                    #except NameError:
                    #    pos = []
                    
                    #try:
                    #    len(annotation)
                    #except NameError:
                    #    annotation = {}
                        
                    annotation['ExonAnnotation'] = ExonAnnotation
                    
                    ### The rest region comaprison 
                    #annotation['RestRegionAnnotation'] = {}
                    All_regions = [int(ind) for ind in list(set(QueryAnnotation['annotation']['PosAnnotation']))]
                    Unchecked_regions = [ind for ind in All_regions if 'Exon'+str(ind) not in ExonAnnotation.keys()]
                    for posInd in range(len(alignment['AlignSymbol'])):
                        ## annotate mismatches between donor and recipient in other regions
                        regionSymbol = QueryAnnotation['annotation']['PosAnnotation'][posInd]
                        if int(regionSymbol) in Unchecked_regions and str(posInd) not in annotation.keys():
                            
                            if alignment[seq1_ID][posInd] != alignment[seq2_ID][posInd]: # missmatch
                                if regionSymbol == '0':
                                    RegionName = '5\'-UTR'
                                elif regionSymbol == min(All_regions):
                                    RegionName = '3\'-UTR'
                                elif regionSymbol > 0:
                                    RegionName = 'Exon'+ regionSymbol
                                elif regionSymbol < 0:
                                    RegionName = 'Intron'+ regionSymbol
                                annotation[str(posInd)] = RegionName + '.' + str(posInd - QueryAnnotation['annotation']['PosAnnotation'].index(regionSymbol) + 1)
                                #AlignSymbol[ind] = 'X'
                    
                    #### : save and print the mismatched positons outside the ARS.
                    
                    if saveFile:
                        save_aln(alignment, CharSpaceNum_perLine, algn_file, annotation)
                    
                    print_MisMatchedPos(alignment, CharSpaceNum_perLine, annotation)
                        
                if noRef: # no reference sequences provided
                    annotation = {}
                    if seq1 != seq2:
                        AddtnlMessage = ''
                        for ind in range(len(HLAtyping)):
                            if ind < len(HLAtyping)-1:
                                AddtnlMessage += HLAtyping[ind]+ ' and '
                            else:
                                AddtnlMessage += HLAtyping[ind]
                        #SeqOne = SeqRecord(Seq(seq1), id = seq1_ID, description = 'The first sequence')
                        #SeqTwo = SeqRecord(Seq(seq2), id = seq2_ID, description = 'The second sequence')
            
                        #records = (SeqOne, SeqTwo)
                        #alignment, pos, CharSpaceNum_perLine = align_seqs(records, algn_file, saveFile)
                        
                        records = {seq1_ID: seq1, seq2_ID: seq2}
                        
                        QueryAnnotation = SeqAnnotation(records, locus)
                        
                        annotation = QueryAnnotation['annotation']
                        CharSpaceNum_perLine = QueryAnnotation['CharSpaceNum_perLine']
                        pos = QueryAnnotation['pos']
                        alignment = QueryAnnotation['alignment']
                        #annotation['SeqFull'] = 'Reference genome sequence is needed to assess.'
                        annotation['Note'] = AddtnlMessage
                        if saveFile:
                            save_aln(alignment, 20, algn_file, annotation)
                            
                        print_MisMatchedPos(alignment, CharSpaceNum_perLine, annotation)
                        
                    else: # if the sequences are the same
                        alignment = "Same sequence"
                        pos = []
                        annotation['SeqFull'] = 'Reference genome sequence is needed to assess.'
                        annotation['Note'] = 'The two sequences are excatly the same.'
                        print('The two sequences are exactly the same!')
                        
            else: # if the sequences are the same
                alignment = "Same sequence"
                pos = []
                annotation = {}
                
        else: # Multiple sequences +/- HLA ref seq alignment and annotation.
            ReferenceSeq = {}
            AddtnlMessage = ''
            potentialKeys = HLAtyping
            for typing in HLAtyping:
                RefKey = 'RefSeq|'+ typing
                tempRefseq = IMGTdbIO.readIMGTsql(typing, field = DB_field)
                if tempRefseq[0] != '':
                    ReferenceSeq[RefKey] = [re.sub('\|', '', tempRefseq[0]), tempRefseq[1]]
                else:
                    AddtnlMessage = AddtnlMessage + "\n" + RefKey + ' doesn\'t have a confirmed genomic sequence.'
                    potentialKeys.remove(typing)

            RefKey = 'RefSeq|'+ potentialKeys[0]
           
            RefAnnotation = ReferenceSeq[RefKey][1].split(" ")
            Ref_frontMisCount = list(ReferenceSeq[RefKey][0])[:505].count('*')
            
            records = ()
            counter = 0
            for key, seq in Sequence.items(): 
                records += (SeqRecord(Seq(seq), id = key, description = 'Sequence order: '+str(counter)), )
                counter += 1
                
            for itemID, itemSeq in ReferenceSeq.items():
                records += (SeqRecord(Seq(itemSeq[0]), id = itemID, description = 'Reference IMGT/HLA sequence.'), )
            
            alignment, pos, CharSpaceNum_perLine = align_seqs(records)
            
            annotation = posAnnotation(alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount)
            annotation['Note'] = AddtnlMessage
            if saveFile:
                save_aln(alignment, CharSpaceNum_perLine, algn_file, annotation)
            
            print_MisMatchedPos(alignment, CharSpaceNum_perLine, annotation)
            
    return(alignment, pos, annotation)

def SeqAnnotation(query, locus, Ref = None):
    '''
    An unknown sequence with no known genomic reference sequence, aligned to a universe reference sequence, 
    and infer the Exon Intron regions.
    A: A*01:01:01:01
    B: B*07:02:01
    C: C*01:02:01
    DRB1: DRB1*01:02:01 Intron1+Exon2+Intron2; Intron2+Exon3+Intron3
    DPB1: DPB1*02:01:02 Intron1+Exon2+Intron2; Intron2+Exon3+Intron3
    DQB1: DQB1*05:01:01:01 Intron1+Exon2+Intron2; Intron2+Exon3+Intron3; or Intron1+Exon2+Intron2+Exon3+Intron3+Exon4+Exon4
    
    Arguments:
         input:
             query: dictionary {'seqID': sequence}
             locus: HLA locus
         output:
             QueryAnnotation: dictionary {'alignment', 'pos', 'CharSpaceNum_perLine', 'annotation'}
    '''
    QueryAnnotation = {}
    DB_field = 'UnalignedGenomSeq, SeqAnnotation'
    if locus == 'A':
        RefKey = 'A*01:01:01:01'
    elif locus == 'B':
        RefKey = 'B*07:02:01'
    elif locus == 'C':
        RefKey = 'C*01:02:01'
    elif locus == 'DRB1':
        RefKey = 'DRB1*01:02:01'
    elif locus == 'DPB1':
        RefKey = 'DPB1*02:01:02'
    elif locus == 'DQB1':
        RefKey = 'DQB1*05:01:01:01'
    

    RefSeq = IMGTdbIO.readIMGTsql(RefKey, field = DB_field)
    
    RefKey = 'RefSeq|' + RefKey
        
    records = ()
    counter = 0
    for key, seq in query.items(): 
        records += (SeqRecord(Seq(seq), id = key, description = 'Sequence order: '+str(counter)), )
        counter += 1
        
    records += (SeqRecord(Seq(re.sub('\|', '', RefSeq[0])), id = RefKey, description = 'Reference IMGT/HLA sequence.'), )
    
    alignment, pos, CharSpaceNum_perLine = align_seqs(records)
    
    RefAnnotation = RefSeq[1].split(' ')
    RefAnnotation.pop(RefAnnotation.index(''))

    RefAlingment = alignment.pop(RefKey,None)
    QueryAnnotation['alignment'] = alignment
    QueryAnnotation['pos'] = pos
    QueryAnnotation['CharSpaceNum_perLine'] = CharSpaceNum_perLine
    
    QueryAnnotation['annotation'] = {}
    if RefAlingment.find('-') == -1: # no insertion in reference, then directly read from the RefAnnotation
        QueryAnnotation['annotation']['PosAnnotation'] = RefAnnotation
    else:
        QueryAnnotation['annotation']['PosAnnotation'] = RefAnnotation
        insertCount = 0 
        for m in re.finditer('-', RefAlingment):
            insertInd = m.start() + insertCount
            QueryAnnotation['annotation']['PosAnnotation'].insert(insertInd, RefAnnotation[m.start()])
        
    Region_lengths = Counter(QueryAnnotation['annotation']['PosAnnotation'])
    AnnNum = [int(item) for item in list(set(RefAnnotation)) if item != '']
    seq_completeness = {}
    for key, seq in query.items():
        seq_completeness[key] = ''
        
        startIndex = re.search('[A-Z]+', QueryAnnotation['alignment'][key]).start()
        endIndex = re.search('[A-Z]+', QueryAnnotation['alignment'] [key]).end()
        
        FisrtRegionSymbol = QueryAnnotation['annotation']['PosAnnotation'][startIndex]
        queryFirstLength = Counter(QueryAnnotation['annotation']['PosAnnotation'][startIndex:])[FisrtRegionSymbol]
    
        LastRegionSymbol = QueryAnnotation['annotation']['PosAnnotation'][endIndex]
        queryLastLength = Counter(QueryAnnotation['annotation']['PosAnnotation'][:endIndex])[LastRegionSymbol]
        
        for MSpos, freq in Region_lengths.items():
            
            if abs(int(MSpos)) < abs(int(FisrtRegionSymbol)) or abs(int(MSpos)) > abs(int(LastRegionSymbol)): # missing full sequence
                if int(MSpos) > 0: # Exon
                    seq_completeness[key] += 'Missing full Exon' + MSpos + '; '
                elif int(MSpos) == 0: # 5'-UTR
                    seq_completeness[key] += 'Missing full 5\'-UTR' + '; '
                      
                elif int(MSpos) == min(AnnNum): # 3'-UTR
                    seq_completeness[key] += 'Missing full 3\'-UTR' + '; '
                
                else: # intron
                    absMSpos = re.sub('-', '', MSpos)
                    seq_completeness[key] += 'Missing full Intron' + absMSpos + '; '

            elif abs(int(MSpos)) == abs(int(FisrtRegionSymbol)): # partial
                if Region_lengths[MSpos] > queryFirstLength:
                    if int(MSpos) > 0: # Exon
                        seq_completeness[key] += 'Missing Partial Exon' + MSpos + '; '
                    elif int(MSpos) == 0: # 5'-UTR
                        seq_completeness[key] += 'Missing Partial 5\'-UTR' + '; '
                          
                    elif int(MSpos) == min(AnnNum): # 3'-UTR
                        seq_completeness[key] += 'Missing Partial 3\'-UTR' + '; '
                    else: # intron
                        absMSpos = re.sub('-', '', MSpos)
                        seq_completeness[key] += 'Missing Partial Intron' + absMSpos + '; '
            elif abs(int(MSpos)) == abs(int(LastRegionSymbol)): # partial
                if Region_lengths[MSpos] > queryLastLength:
                    if int(MSpos) > 0: # Exon
                        seq_completeness[key] += 'Missing Partial Exon' + MSpos + '; '
                    elif int(MSpos) == 0: # 5'-UTR
                        seq_completeness[key] += 'Missing Partial 5\'-UTR' + '; '
                          
                    elif int(MSpos) == min(AnnNum): # 3'-UTR
                        seq_completeness[key] += 'Missing Partial 3\'-UTR' + '; '
                    else: # intron
                        absMSpos = re.sub('-', '', MSpos)
                        seq_completeness[key] += 'Missing Partial Intron' + absMSpos + '; '
 
    QueryAnnotation['annotation']['SeqFull'] = seq_completeness
    QueryAnnotation['annotation']['Note'] = Ref + ' does not have a confirmed genomic sequence.'
    tempSymbol = list(QueryAnnotation['alignment']['AlignSymbol'])
    for MMpos in pos:
        
        tempRead = []
        for key, seq in query.items():
            tempRead.append(QueryAnnotation['alignment'][key][MMpos])
            
        if len(list(set(tempRead))) != 1: # some mis-match; real mismatched position among the query sequences
            
            if tempSymbol[MMpos] != 'X':
                tempSymbol[MMpos] = 'X'
                
            TempAnn = QueryAnnotation['annotation']['PosAnnotation'][MMpos]
            if int(TempAnn) > 0: # Exon
                QueryAnnotation['annotation'][str(MMpos)] = 'Exon' + TempAnn + '.' + str(MMpos - QueryAnnotation['annotation']['PosAnnotation'].index(TempAnn) + 1)
            elif int(TempAnn) == 0: # 5'-UTR
                QueryAnnotation['annotation'][str(MMpos)] = '5\'-UTR' + '.' + str(MMpos - QueryAnnotation['annotation']['PosAnnotation'].index(TempAnn) + 1)
            elif int(TempAnn) == min(AnnNum): # 3'-UTR
                QueryAnnotation['annotation'][str(MMpos)] = '3\'-UTR' + '.' + str(MMpos - QueryAnnotation['annotation']['PosAnnotation'].index(TempAnn) + 1)
            else: # intron
                absTempAnn = re.sub('-', '', TempAnn)
                QueryAnnotation['annotation'][str(MMpos)] = 'Intron' + absTempAnn + '.' + str(MMpos - QueryAnnotation['annotation']['PosAnnotation'].index(TempAnn) + 1)
        
        else: ## not a real mismatche among the queries
            if tempSymbol[MMpos] == 'X':
                tempSymbol[MMpos] = '.'
                
    QueryAnnotation['alignment']['AlignSymbol'] = "".join(tempSymbol)
    
    return(QueryAnnotation)
        
def posAnnotation(alignment, pos, RefKey, RefAnnotation, Ref_frontMisCount):
    '''
    output - 
       annotation: dictionary - {'SeqFull': , 'Note':, 'Pos':}
    
    '''
    annotation = {}
    
    AnnNum = [int(item) for item in list(set(RefAnnotation)) if item != '']
    for MMpos in pos:
        
        TempAnn = RefAnnotation[MMpos-alignment[RefKey][Ref_frontMisCount:MMpos].count('-')]
        if int(TempAnn) > 0: # Exon
            annotation[str(MMpos)] = 'Exon' + TempAnn + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)
        elif int(TempAnn) == 0: # 5'-UTR
            annotation[str(MMpos)] = '5\'-UTR' + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)
        elif int(TempAnn) == min(AnnNum): # 3'-UTR
            annotation[str(MMpos)] = '3\'-UTR' + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)
        else: # intron
            absTempAnn = re.sub('-', '', TempAnn)
            annotation[str(MMpos)] ='Intron' + absTempAnn + '.' + str(MMpos - RefAnnotation.index(TempAnn) + 1)
    ## Missing parts:
    Ann_Pos_frequency = Counter(RefAnnotation)
    if alignment['AlignSymbol'].find('>') != -1: # Partial length
        all_missing_pos = [ind for ind in range(len(alignment['AlignSymbol'])) if alignment['AlignSymbol'][ind] == '>']
        all_missing_ann = [RefAnnotation[ind-alignment[RefKey][Ref_frontMisCount:ind].count('-')] for ind in all_missing_pos]
            
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
                elif RefAnnotation.index(str(min(AnnNum))) in all_missing_pos:
                    temp_missing_info = temp_missing_info + 'Missing full 3\'-UTR' + '; ' 
                else: # partial reigon
                    temp_missing_info = temp_missing_info + 'Missing Partial 3\'-UTR' + '; '
            else: # intron
                absMSpos = re.sub('-', '', MSpos)
                if missing_ann_frequency[MSpos] == Ann_Pos_frequency[MSpos]: # missing full region
                    temp_missing_info = temp_missing_info + 'Missing full Intron' + absMSpos + '; '
                else: # partial reigon
                    temp_missing_info = temp_missing_info + 'Missing Partial Intron'+ absMSpos + '; '
       
        annotation['SeqFull'] = temp_missing_info
        
    else:
        annotation['SeqFull'] = 'Completely full length sequences.'
        
    return(annotation)

def AlignExons(ExonSequences, align_file = None, saveFile = False, method = 'muscle'):
    '''
    Align each exon reference sequence to the query, if there is no genomic sequences available.

    ExonAnnotation = AlignExons(ExonSeqs, CharSpaceNum_perLine, algn_file, saveFile)
    '''
    
    ExonAlignment = {}
    for ExonID, seqObj in ExonSequences.items():
        records = ()
        counter = 0
        ExonAlignment[ExonID] ={}
        for itemID, itemSeq in seqObj.items():
            records +=  (SeqRecord(Seq(itemSeq), id = itemID, description = 'Sequence order: '+str(counter)), )
        
        alignment, pos, CharSpaceNum_perLine = align_seqs(records, method = 'muscle')
        alignment['SeqFull'] = 'Exon Region only comparison'
        alignment['Note'] = 'No genomic Reference sequence available. Checking ARS...'
        
        if saveFile:
            save_aln(alignment, CharSpaceNum_perLine, align_file)
        
        ExonAlignment[ExonID]['alignment'] = alignment
        ExonAlignment[ExonID]['pos'] = pos
        ExonAlignment[ExonID]['CharSpaceNum_perLine'] = CharSpaceNum_perLine
       
    return(ExonAlignment)
        
def align_seqs(sequences, align_file = None, saveFile = False, method = 'muscle', temp_fp = "../Output/Muscle_temp/"):
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

def print_MisMatchedPos(reformated_algn, CharSpaceNum_perLine, Annotation = None):
    '''
    '''
    MM_positions = [pos for pos, char in enumerate(reformated_algn['AlignSymbol']) if char == 'X']
    if Annotation == None:
        for ind in MM_positions:
            print("="*12)
            print("Position "+str(ind)+":")
            for key in reformated_algn.keys():
                if key != 'AlignSymbol':
                    print(key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_algn[key][ind])
            print("="*12)
    elif len(Annotation)>0:
        if type(Annotation['SeqFull']) != dict:
            print("Sequence Completeness: " + Annotation['SeqFull'])
        else: 
            for key, item in Annotation['SeqFull'].items():
                print("Sequence Completeness:")
                print('>>' + key + ":\n"+ item)
                print("-"*12)
        for ind in MM_positions:
            print("="*12)
            print("Position "+str(ind)+":")
            print("Annotation: " + Annotation[str(ind)])
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
            if MM_LastPos != len(reformated_alignment['AlignSymbol'])-1 and reformated_alignment['AlignSymbol'][MM_LastPos+1] == '>':
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
        save_aln(reformated_alignment, CharSpaceNum_perLine, output)
        
    return(reformated_alignment, CharSpaceNum_perLine)
 
def save_aln(reformated_alignment, CharSpaceNum_perLine, output_fp, Annotation = None):

    MM_positions = [pos for pos, char in enumerate(reformated_alignment['AlignSymbol']) if char == 'X']
    
    with open(output_fp, 'a') as f:
        for key in reformated_alignment.keys():
            if key != 'AlignSymbol':
                f.write(key + ' '*(CharSpaceNum_perLine - len(key)) + reformated_alignment[key]+'\n')
        f.write('Align' + ' '*(CharSpaceNum_perLine - len('Align')) + reformated_alignment['AlignSymbol']+'\n')
        
        if Annotation == None:
            for ind in MM_positions:
                f.write("="*12+'\n')
                f.write("Position "+str(ind)+":" + "\n")
                for key in reformated_alignment.keys():
                    if key != 'AlignSymbol':
                        f.write('>>' + key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_alignment[key][ind] + '\n')
                f.write("="*12 + "\n")
        elif len(Annotation)>2: 
            if type(Annotation['SeqFull']) != dict:
                f.write("Sequence Completeness: " + Annotation['SeqFull'] + "\n")
            else: 
                for key, item in Annotation['SeqFull'].items():
                    f.write("Sequence Completeness: \n")
                    f.write('>>' + key + ":\n"+ item + '\n')
                    f.write("-"*12 + "\n")
                    
            f.write("Note: " + Annotation['Note'] + "\n")
            for ind in MM_positions:
                f.write("="*12+ "\n")
                f.write("Position "+str(ind)+":\n")
                f.write("Annotation: " + Annotation[str(ind)] + "\n")
                for key in reformated_alignment.keys():
                    if key != 'AlignSymbol':
                        f.write(key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_alignment[key][ind]+ '\n')
                f.write("="*12+ "\n\n") 
        elif len(Annotation) == 2 and list(Annotation.keys()) == ['SeqFull', 'Note']: # no mismatch
            f.write("Sequence Completeness: " + Annotation['SeqFull'] + "\n")
            f.write("Note: " + Annotation['Note'] + "\n")
            for ind in MM_positions:
                f.write("="*12+'\n')
                f.write("Position "+str(ind)+":" + "\n")
                for key in reformated_alignment.keys():
                    if key != 'AlignSymbol':
                        f.write(key + ":"+ ' '*(CharSpaceNum_perLine - len(key)) +reformated_alignment[key][ind] + '\n')
                f.write("="*12 + "\n")

def swapPS_comparison(Seq1, params1, Seq2, params2, caseID):
    '''
    If the two phase set mismatch numbers are too high, then try the swapped Typing alignment; 
    Might be the Phase set swapped case.
    '''

    algn_file = re.sub('.aln', '_PSswapped.aln', params1['algn_file'])
    saveFile = params1['saveFile']
    locus = params1['HLAtyping'][0].split('*')[0]
    swapped_alignment = {}
    
    DB_field = 'UnalignedGenomSeq, SeqAnnotation'

    if locus in ['A', 'C']:
        DB_BackUpfield = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7, Exon8'
    elif locus == 'B':
        DB_BackUpfield = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7'
    elif locus in ['DRB1', 'DPB1']:
        DB_BackUpfield = 'Exon2, Exon3'
    elif locus == 'DQB1':
        DB_BackUpfield = 'Exon2, Exon3, Exon4'
    
    for ind in range(2):
        
        if ind == 0:
            
            for item in list(Seq1.keys()):
                if 'Recipient' in item:
                    seq1_ID = item
                    seq1 = Seq1[item]
            
            for item in list(Seq2.keys()):
                if 'Donor' in item:
                    seq2_ID = item
                    seq2 = Seq2[item]
            if len(params2['HLAtyping'])>1:
                tplist = [params1['HLAtyping'][0], params2['HLAtyping'][1]]
            else:
                tplist = [params1['HLAtyping'][0], params2['HLAtyping'][0]]
                
            HLAtyping = []
            for tp in tplist:
                if tp.find('[') == -1:
                    if tp.find('/') != -1:
                        ambTPlist = tp.split('/')
                        HLAtyping.extend(ambTPlist)
                    else:
                        HLAtyping.append(tp)
                else:
                    possTPlist = re.sub('[\[\'\]]', '',tp) # remove possible characters
                    possTPlist = possTPlist.split(",")
                    
                    for item in possTPlist:
                        if item.find('/') != -1:
                            item_pos = item.replace(" ", "")
                            ambTPlist = item_pos.split('/')
                            HLAtyping.extend(ambTPlist)
                        else:    
                            #HLAtyping.extend(possTPlist)
                            HLAtyping.append(item.replace(" ", ""))
        else: 
            for item in list(Seq1.keys()):
                if 'Donor' in item:
                    seq1_ID = item
                    seq1 = Seq1[item]
            
            for item in list(Seq2.keys()):
                if 'Recipient' in item:
                    seq2_ID = item
                    seq2 = Seq2[item]
            if len(params1['HLAtyping']) > 1:
                tplist = [params1['HLAtyping'][1], params2['HLAtyping'][0]]
            else:
                tplist = [params1['HLAtyping'][0], params2['HLAtyping'][0]]
            
            HLAtyping = []
            for tp in tplist:
                if tp.find('[') == -1:
                    if tp.find('/') != -1:
                        ambTPlist = tp.split('/')
                        HLAtyping.extend(ambTPlist)
                    else:
                        HLAtyping.append(tp)
                else:
                    possTPlist = re.sub('[\[\'\]]', '',tp) # remove possible characters
                    possTPlist = possTPlist.split(",")
                    
                    for item in possTPlist:
                        if item.find('/') != -1:
                            item_pos = item.replace(" ", "")
                            ambTPlist = item_pos.split('/')
                            HLAtyping.extend(ambTPlist)
                        else:    
                            #HLAtyping.extend(possTPlist)
                            HLAtyping.append(item.replace(" ", ""))
                   # HLAtyping.append(tp)
                   
        Sequence= {seq1_ID: seq1, seq2_ID:seq2}
        
        if locus in ['A', 'B', 'C']:
            params = {'algn_file': algn_file, 'saveFile': saveFile, 'HLAtyping': HLAtyping, 'DB_field': DB_field, "DB_BackUpfield": DB_BackUpfield}
            print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
            print(params)
            
            alignment, pos, annotation = compare_seqs(Sequence, params)
            
            ## save results
            saveOBJ = {'seq': Sequence, 'params': params, 'alignment':alignment, 'MMannotation': annotation, 'MMpos': pos}
            
            Output_fname = re.sub('.aln','_annotation_swappedPS'+str(ind+1), algn_file)
            IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
            
            swapped_alignment['PS'+str(ind+1)] = saveOBJ
        elif locus in ['DRB1', 'DPB1']:
            
            params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
            print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
            print(params)
            
            alignment, annotation = compare_Targeted_Region(Sequence, params)
            
            ## save results # for multiple Exons
            for itemID, itemDict in alignment.items():
                
                saveOBJ = {'seq': Sequence, 'params': params, 'alignment':itemDict, 'MMannotation': annotation[itemID], 'SameSeqs': annotation[itemID]['SameSeqs']}
                
                if annotation[itemID]['SameSeqs']: # same seqs
                    Output_fname = re.sub('.aln','_annotation_swappedPS'+str(ind+1)+'_'+itemID+'_SameSeqs', algn_file)
                else:
                    Output_fname = re.sub('.aln','_annotation_swappedPS'+str(ind+1)+'_'+itemID+'_MisMatchSeqs', algn_file) 
                IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
                
                swapped_alignment['PS'+str(ind+1)+'_'+itemID] = saveOBJ
                
        elif locus in ['DQB1']:
            
            params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
            print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
            print(params)
            
            alignment, annotation = compare_DQB1_Targeted_Region(Sequence, params)
            
            if any("Exon" in s for s in alignment.keys()):
                ## save results # for multiple Exons
                for itemID, itemDict in alignment.items():
                    
                    saveOBJ = {'seq': Sequence, 'params': params, 'alignment':itemDict, 'MMannotation': annotation[itemID], 'SameSeqs': annotation[itemID]['SameSeqs']}
                    
                    if annotation[itemID]['SameSeqs']: # same seqs
                        Output_fname = re.sub('.aln','_annotation_swappedPS'+str(ind+1)+'_'+itemID+'_SameSeqs', algn_file)
                    else:
                        Output_fname = re.sub('.aln','_annotation_swappedPS'+str(ind+1)+'_'+itemID+'_MisMatchSeqs', algn_file) 
                    IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
                    swapped_alignment['PS1'+str(ind+1)+'_'+itemID] = saveOBJ
            else: 
                ## save results -- for one single sequence
                saveOBJ = {'seq': Sequence, 'params': params, 'alignment':alignment, 'MMannotation': annotation}
            
                Output_fname = re.sub('.aln','_annotation_swappedPS'+str(ind+1), algn_file)
                IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
                swapped_alignment['PS1'+str(ind)] = saveOBJ

    return(swapped_alignment)

def RegionCount(MMAnnotation, locus):
    '''
    Count the numbers of mismatches in ARS, non-ARS exon and Intron region
    '''
    RegionMMCount = {'ARS': 0, 'Non_ARS_exon': 0, 'Intron': 0}
    if locus in ['A', 'B', 'C']:
        ARS = ['Exon2', 'Exon3']
    elif locus in ['DRB1', 'DQB1', 'DPB1']:
        ARS = ['Exon2']
    
    for key, item in MMAnnotation.items():
        if key.isdigit():
            if item.split('.')[0] in ARS:
                RegionMMCount['ARS'] += 1
            elif 'Exon' in item.split('.')[0]:
                RegionMMCount['Non_ARS_exon'] += 1
            elif 'Intron' in item.split('.')[0]:
                RegionMMCount['Intron'] += 1
            #elif 'UTR' in
     
    return(RegionMMCount)          
                
                
