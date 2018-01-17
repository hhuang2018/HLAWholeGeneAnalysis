#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 10:30:43 2017

@author: hhuang2
"""

# import glob
import sqlite3 as sql
# from utils import phase_block_check as ps
from utils import IMGTdbIO, CompareSeq
import os
import re


locus = 'A'

pkl_fp = '../Output/SG41_52/SG41_52_DRpairs/SG41_52_HLA_'+ locus +'_paired.pkl'

DRpair_seqInfo = IMGTdbIO.load_pickle2dict(pkl_fp)

case_count = len(DRpair_seqInfo)
print('Locus '+locus+ ' has '+ str(case_count) + ' paired cases.')

DB_fp = "../Output/SG41_52/SG41_52_DRpairs/SG41_52_HLA_"+ locus +"_paired.db"
conn = sql.connect(DB_fp)
cursor = conn.cursor()
cursor.execute('''CREATE TABLE IF NOT EXISTS DR_pair_comparison 
               (BMT_caseID text, QC text, 
               PS1_HLATyping text, 
               PS1_GLstringM text, PS1_SeqM text, 
               PS2_HLATyping text,
               PS2_GLstringM text, PS2_SeqM text,
               Audit text, Active text, Comment text)''')

conn.commit()
for caseID, SeqInfo in DRpair_seqInfo.items():
    Active = SeqInfo['Active']
    Audit = SeqInfo['Audit']
    Comment = SeqInfo['Comment']
    QC = SeqInfo['QC']
    
    PS1_HLATyping = SeqInfo['PS1']['HLATyping']
    PS2_HLATyping = SeqInfo['PS2']['HLATyping']
    
    if QC == 'PASS':
        PS1_GLstringM = 'Y'
        PS2_GLstringM = 'Y'
    else:
        PS12_qc = QC.split(";")
        if PS12_qc[0] == 'PASS':
            PS1_GLstringM = 'Y'
        else:
            PS1_GLstringM = 'N'
            
        if PS12_qc[1] == 'PASS':
            PS2_GLstringM = 'Y'
        else:
            PS2_GLstringM = 'N'
    
    if SeqInfo['PS1']['Recipient'] == SeqInfo['PS1']['Donor']:
        PS1_SeqM = 'Y'
    else: 
        PS1_SeqM = 'N'
    DRpair_seqInfo[caseID]['PS1']['isSeqMatch'] = PS1_SeqM
        
    if SeqInfo['PS2']['Recipient'] == SeqInfo['PS2']['Donor']:
        PS2_SeqM = 'Y'
    else: 
        PS2_SeqM = 'N'
        
    DRpair_seqInfo[caseID]['PS2']['isSeqMatch'] = PS2_SeqM
    
    record = (caseID, QC, PS1_HLATyping, PS1_GLstringM, PS1_SeqM, PS2_HLATyping, PS2_GLstringM, PS2_SeqM, Active, Audit, Comment,)
    cursor.execute('INSERT INTO DR_pair_comparison VALUES (?,?,?,?,?,?,?,?,?,?,?)', record)

conn.commit()   
conn.close()

fname = '../Output/SG41_52/SG41_52_DRpairs/SG41_52_HLA_' + locus + '_wComparison'
IMGTdbIO.save_dict2pickle(DRpair_seqInfo, fname)

########### check GL-string match, sequence matching
DRpair_seqInfo = {}
Loci = ['A','B', 'C', 'DRB1', 'DQB1', 'DPB1']

All_caseIDs = []
for locus  in Loci:
    fname = '../Output/SG41_52/SG41_52_DRpairs/SG41_52_HLA_' + locus + '_wComparison.pkl'
    DRpair_seqInfo[locus] = IMGTdbIO.load_pickle2dict(fname)
    All_caseIDs += list(DRpair_seqInfo[locus].keys())
All_caseIDs = list(set(All_caseIDs)) # 3412 total

Matching_cases_stats = {"ClassI_paired":[], "fiveLoci_paired":[], "All_paired":[], 
                        "A_both_SeqMatch":[], "A_both_Seqmm":[], "A_one_Seqmm":[], 'A_GLmm_SeqM':[], 'A_GlM_Seqmm':[],
                        "B_both_SeqMatch":[], "B_both_Seqmm":[], "B_one_Seqmm":[], 'B_GLmm_SeqM':[], 'B_GlM_Seqmm':[],
                        "C_both_SeqMatch":[], "C_both_Seqmm":[], "C_one_Seqmm":[], 'C_GLmm_SeqM':[], 'C_GlM_Seqmm':[],
                        "DRB1_both_SeqMatch":[], "DRB1_both_Seqmm":[], "DRB1_one_Seqmm":[], 'DRB1_GLmm_SeqM':[], 'DRB1_GlM_Seqmm':[],
                        "DPB1_both_SeqMatch":[], "DPB1_both_Seqmm":[], "DPB1_one_Seqmm":[], 'DPB1_GLmm_SeqM':[], 'DPB1_GlM_Seqmm':[],
                        "DQB1_both_SeqMatch":[], "DQB1_both_Seqmm":[], "DQB1_one_Seqmm":[], 'DQB1_GLmm_SeqM':[], 'DQB1_GlM_Seqmm':[]}
for caseID in All_caseIDs:
    
    classI_flag = True
    fiveLoci_flag = True
    All_flag = True
    
    for locus in Loci:
        isExist = caseID in DRpair_seqInfo[locus].keys()
        All_flag = All_flag and isExist
        if locus in ['A', 'B', 'C']:
            classI_flag = classI_flag and isExist
        if locus != 'DPB1':
            fiveLoci_flag = fiveLoci_flag and isExist
        
        if isExist:
            if DRpair_seqInfo[locus][caseID]['PS1']['isSeqMatch'] == 'Y' and DRpair_seqInfo[locus][caseID]['PS2']['isSeqMatch'] == 'Y':
                Matching_cases_stats[locus+'_both_SeqMatch'].append(caseID)
                if DRpair_seqInfo[locus][caseID]['QC'] != 'PASS':
                    Matching_cases_stats[locus+'_GLmm_SeqM'].append(caseID)
            elif DRpair_seqInfo[locus][caseID]['PS1']['isSeqMatch'] != 'Y' and DRpair_seqInfo[locus][caseID]['PS2']['isSeqMatch'] != 'Y':
                Matching_cases_stats[locus+'_both_Seqmm'].append(caseID)
                if DRpair_seqInfo[locus][caseID]['QC'] == 'PASS':
                    Matching_cases_stats[locus+'_GlM_Seqmm'].append(caseID)
            else: # one mismatch
                Matching_cases_stats[locus+'_one_Seqmm'].append(caseID)
                if DRpair_seqInfo[locus][caseID]['QC'] == 'PASS':
                    Matching_cases_stats[locus+'_GlM_Seqmm'].append(caseID)
    if All_flag:
        Matching_cases_stats["All_paired"].append(caseID)
    
    if classI_flag:
        Matching_cases_stats["ClassI_paired"].append(caseID)
        
    if fiveLoci_flag:
        Matching_cases_stats["fiveLoci_paired"].append(caseID)
    
fname = '../Output/SG41_52/SG41_52_DRpair_Stats/SG41_52_pairedCases_Stats'
IMGTdbIO.save_dict2pickle(Matching_cases_stats, fname)

print("Paired at all 6 loci cases: " + str(len(Matching_cases_stats['All_paired'])))
print("Paired at all 5 loci cases: "+ str(len(Matching_cases_stats['fiveLoci_paired'])))
print("Paired at Class I loci cases: "+ str(len(Matching_cases_stats['ClassI_paired'])))
for locus in Loci:
    print("Both sequences matched at locus " + locus + ": " + str(len(Matching_cases_stats[locus+'_both_SeqMatch'])))

for locus in Loci:
    print("Both sequences MisMatched at locus " + locus + ": " + str(len(Matching_cases_stats[locus+'_both_Seqmm'])))

for locus in Loci:
    print("One sequence MisMatched at locus " + locus + ": " + str(len(Matching_cases_stats[locus+'_one_Seqmm'])))
    
for locus in Loci:
    print("One or both sequences matched but GL-string mismatched at locus " + locus + ": " + str(len(Matching_cases_stats[locus+'_GLmm_SeqM'])))

for locus in Loci:
    print("One or both sequences MisMatched but GL-string matched at locus " + locus + ": " + str(len(Matching_cases_stats[locus+'_GlM_Seqmm'])))

############################
# Check unmatched cases
############################

fname = '../Output/SG41_52/SG41_52_DRpair_Stats/SG41_52_pairedCases_Stats.pkl'
Matching_cases_stats = IMGTdbIO.load_pickle2dict(fname)

locus = 'A'
# both mismatched sequences
MM_caseID = Matching_cases_stats[locus+'_both_Seqmm']

DB_fp = '../Output/SG41_52/SG41_52_DRpairs/SG41_52_HLA_'+ locus +'_paired.db'
con = sql.connect(DB_fp)
con.row_factory = sql.Row
cur = con.cursor()

bothMM_output = "../Output/SG41_52/SG41_52_bothMisMatched_locus_" + locus + "_1218_TargetedAlignment/"
if not os.path.exists(bothMM_output):
    os.makedirs(bothMM_output) 

### 
DB_field = 'UnalignedGenomSeq, SeqAnnotation'

if locus in ['A', 'C']:
    DB_BackUpfield = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7, Exon8'
elif locus == 'B':
    DB_BackUpfield = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7'
elif locus in ['DRB1', 'DPB1']:
    DB_BackUpfield = 'Exon2, Exon3'
elif locus == 'DQB1':
    DB_BackUpfield = 'Exon2, Exon3, Exon4'

## Class I
for caseID in MM_caseID:
    algn_file = bothMM_output + 'CaseID_' + caseID + '_' + locus + '_aligned.aln'
    
    if not os.path.exists(algn_file):
        t = (caseID,)
        cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
        case_records = cur.fetchall()
        for ind in range(2):
            seq1_ID = 'Recipient-PS'+str(ind+1)
            seq2_ID = 'Donor-PS'+str(ind+1)
            seq1 = case_records[ind][seq1_ID.split('-')[0]]
            seq2 = case_records[ind][seq2_ID.split('-')[0]]
            HLAtyping_list = case_records[ind]['HLATyping']
            tplist = HLAtyping_list.split("+")
            HLAtyping = []
            for tp in tplist:
                if tp.find('/') != -1:
                    ambTPlist = tp.split('/')
                    HLAtyping.extend(ambTPlist)
                else:
                    HLAtyping.append(tp)
            
            Sequence= {seq1_ID: seq1, seq2_ID:seq2}
    
            params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping, 'DB_field': DB_field, "DB_BackUpfield": DB_BackUpfield}
            print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
            print(params)
            alignment, pos, annotation = CompareSeq.compare_seqs(Sequence, params)
            
            ## save results
            saveOBJ = {'seq': Sequence, 'params': params, 'alignment':alignment, 'MMannotation': annotation, 'MMpos': pos}
            
            Output_fname = bothMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)
            IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
con.close()

## Class II
locus = 'DPB1'
for caseID in MM_caseID:
    algn_file = bothMM_output + 'CaseID_'+ caseID + '_' + locus + '_aligned.aln'
    
    if not os.path.exists(algn_file):
        t = (caseID,)
        cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
        case_records = cur.fetchall()
        for ind in range(2):
            seq1_ID = 'Recipient-PS'+str(ind+1)
            seq2_ID = 'Donor-PS'+str(ind+1)
            seq1 = case_records[ind][seq1_ID.split('-')[0]]
            seq2 = case_records[ind][seq2_ID.split('-')[0]]
            HLAtyping_list = case_records[ind]['HLATyping']
            tplist = HLAtyping_list.split("+")
            HLAtyping = []
            for tp in tplist:
                if tp.find('/') != -1:
                    ambTPlist = tp.split('/')
                    HLAtyping.extend(ambTPlist)
                elif tp.find(',')!= -1:
                    tp = tp.replace('[\'', '')
                    tp = tp.replace('\']', '')
                    tp = tp.replace('\'', '')
                    ambTPlist = tp.split(', ')
                    HLAtyping.extend(ambTPlist)
                else:
                    HLAtyping.append(tp)
            
            Sequence= {seq1_ID: seq1, seq2_ID:seq2}
    
            params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
            print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
            print(params)
            alignment, annotation = CompareSeq.compare_Targeted_Region(Sequence, params)
            
            ## save results # for multiple Exons
            for itemID, itemDict in alignment.items():
                
                saveOBJ = {'seq': Sequence, 'params': params, 'alignment':itemDict, 'MMannotation': annotation[itemID], 'SameSeqs': annotation[itemID]['SameSeqs']}
                
                if annotation[itemID]['SameSeqs']: # same seqs
                    Output_fname = bothMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_SameSeqs'
                else:
                    Output_fname = bothMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_MisMatchSeqs'
                IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)

            
con.close()

#### DQB1
locus = 'DQB1'
for caseID in MM_caseID:
    algn_file = bothMM_output + 'CaseID_'+ caseID + '_' + locus + '_aligned.aln'
    
    if not os.path.exists(algn_file):
        t = (caseID,)
        cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
        case_records = cur.fetchall()
        for ind in range(2):
            
            seq1_ID = 'Recipient-PS'+str(ind+1)
            seq2_ID = 'Donor-PS'+str(ind+1)
            seq1 = case_records[ind][seq1_ID.split('-')[0]]
            seq2 = case_records[ind][seq2_ID.split('-')[0]]
            HLAtyping_list = case_records[ind]['HLATyping']

            #HLAtyping_list_1 = case_records[ind]['HLATyping']
            #tplist_1 = HLAtyping_list_1.split("+")
            #HLAtyping_list_2 = case_records[ind-1]['HLATyping']
            #tplist_2 = HLAtyping_list_2.split("+")
            #tplist = [tplist_1[0], tplist_2[1]]
            tplist = HLAtyping_list.split("+")
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
    
            params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
            print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
            print(params)
            alignment, annotation = CompareSeq.compare_DQB1_Targeted_Region(Sequence, params)
            
            if any("Exon" in s for s in alignment.keys()):
                ## save results # for multiple Exons
                for itemID, itemDict in alignment.items():
                    
                    saveOBJ = {'seq': Sequence, 'params': params, 'alignment':itemDict, 'MMannotation': annotation[itemID], 'SameSeqs': annotation[itemID]['SameSeqs']}
                    
                    if annotation[itemID]['SameSeqs']: # same seqs
                        Output_fname = bothMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_SameSeqs'
                    else:
                        Output_fname = bothMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_MisMatchSeqs'
                    IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
            else: 
                ## save results -- for one single sequence
                saveOBJ = {'seq': Sequence, 'params': params, 'alignment':alignment, 'MMannotation': annotation}
            
                Output_fname = bothMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)
                IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
        
        
con.close()

############################
# Check unmatched sequences
############################
# single mismatched sequences
fname = '../Output/SG41_52/SG41_52_DRpair_Stats/SG41_52_pairedCases_Stats.pkl'
Matching_cases_stats = IMGTdbIO.load_pickle2dict(fname)

locus = 'DQB1'
singleMM_caseID = Matching_cases_stats[locus+'_one_Seqmm']

DB_fp = '../Output/SG41_52/SG41_52_DRpairs/SG41_52_HLA_'+ locus +'_paired.db'
con = sql.connect(DB_fp)
con.row_factory = sql.Row
cur = con.cursor()

singleMM_output = "../Output/SG41_52/SG41_52_singleMisMatched_" + locus + "_1220_TargetedAlignment/"
if not os.path.exists(singleMM_output):
    os.makedirs(singleMM_output) 

## Class I
DB_field = 'UnalignedGenomSeq, SeqAnnotation'
### 
if locus in ['A', 'C']:
    DB_BackUpfield = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7, Exon8'
elif locus == 'B':
    DB_BackUpfield = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7'
elif locus in ['DRB1', 'DPB1']:
    DB_BackUpfield = 'Exon2, Exon3'
elif locus == 'DQB1':
    DB_BackUpfield = 'Exon2, Exon3, Exon4'
    
for caseID in singleMM_caseID:
    algn_file = singleMM_output + caseID + '_' + locus + '_aligned.aln'
    
    if not os.path.exists(algn_file):
        t = (caseID,)
        cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
        case_records = cur.fetchall()

        #if case_records[0]['QC'].split("; ")[0] == 'PASS':
        #    ind = 1
        #else:
        #    ind = 0
        
        cur.execute('SELECT * FROM DR_pair_comparison WHERE BMT_caseID = ?', t)
        case_comparison_records = cur.fetchall()
        if case_comparison_records[0]['PS1_SeqM'] == 'N':
            ind = 0
        elif case_comparison_records[0]['PS2_SeqM'] == 'N':
            ind = 1
        else: 
            ind = -1
        
        seq1_ID = 'Recipient-PS'+str(ind+1)
        seq2_ID = 'Donor-PS'+str(ind+1)
        seq1 = case_records[ind][seq1_ID.split('-')[0]]
        seq2 = case_records[ind][seq2_ID.split('-')[0]]
        HLAtyping_list = case_records[ind]['HLATyping']
        tplist = HLAtyping_list.split("+")
        HLAtyping = []
        for tp in tplist:
            if tp.find('/') != -1:
                ambTPlist = tp.split('/')
                HLAtyping.extend(ambTPlist)
            elif tp.find(',') != -1:
                ambTPlist = re.sub('[\[\'\]]', '',tp)
                ambTPlist = ambTPlist.split(', ')
                HLAtyping.extend(ambTPlist)
            else:
                HLAtyping.append(tp)
                
        Sequence= {seq1_ID: seq1, seq2_ID:seq2}
        #params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping, 'DB_field': DB_field}
        params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping, 'DB_field': DB_field, "DB_BackUpfield": DB_BackUpfield}
        print('CaseID: '+ caseID + ' ; Locus: ' + locus)
        alignment, pos, annotation = CompareSeq.compare_seqs(Sequence, params)

        ## save results
        saveOBJ = {'seq': Sequence, 'params': params, 'alignment':alignment, 'MMannotation': annotation, 'MMpos': pos}
        
        Output_fname = singleMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation'
        IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
con.close()

## Class II: DRB1 and DPB1
locus = 'DPB1'
for caseID in singleMM_caseID:
    algn_file = singleMM_output + 'CaseID_'+ caseID + '_' + locus + '_aligned.aln'
    
    if not os.path.exists(algn_file):
        t = (caseID,)
        cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
        case_records = cur.fetchall()
        #if case_records[0]['QC'].split("; ")[0] == 'PASS':
        #    ind = 1
        #else:
        #    ind = 0
        cur.execute('SELECT * FROM DR_pair_comparison WHERE BMT_caseID = ?', t)
        case_comparison_records = cur.fetchall()
        if case_comparison_records[0]['PS1_SeqM'] == 'N':
            ind = 0
        elif case_comparison_records[0]['PS2_SeqM'] == 'N':
            ind = 1
        else: 
            ind = -1
            
        seq1_ID = 'Recipient-PS'+str(ind+1)
        seq2_ID = 'Donor-PS'+str(ind+1)
        seq1 = case_records[ind][seq1_ID.split('-')[0]]
        seq2 = case_records[ind][seq2_ID.split('-')[0]]
        HLAtyping_list = case_records[ind]['HLATyping']
        tplist = HLAtyping_list.split("+")
        HLAtyping = []
        for tp in tplist:
            if tp.find('/') != -1:
                ambTPlist = tp.split('/')
                HLAtyping.extend(ambTPlist)
            else:
                HLAtyping.append(tp)
        
        Sequence= {seq1_ID: seq1, seq2_ID:seq2}

        params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
        print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
        print(params)
        alignment, annotation = CompareSeq.compare_Targeted_Region(Sequence, params)
        
        ## save results
        for itemID, itemDict in alignment.items():
            
            saveOBJ = {'seq': Sequence, 'params': params, 'alignment':itemDict, 'MMannotation': annotation[itemID], 'SameSeqs': annotation[itemID]['SameSeqs']}
            
            if annotation[itemID]['SameSeqs']: # same seqs
                Output_fname = singleMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_SameSeqs'
            else:
                Output_fname = singleMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_MisMatchSeqs'
            IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
            
con.close()


#### DQB1
locus = 'DQB1'
for caseID in singleMM_caseID:
    algn_file = singleMM_output + 'CaseID_'+ caseID + '_' + locus + '_aligned.aln'
    
    if not os.path.exists(algn_file):
        t = (caseID,)
        cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
        case_records = cur.fetchall()
        #if case_records[0]['QC'].split("; ")[0] == 'PASS':
        #    ind = 1
        #else:
        #    ind = 0
        cur.execute('SELECT * FROM DR_pair_comparison WHERE BMT_caseID = ?', t)
        case_comparison_records = cur.fetchall()
        if case_comparison_records[0]['PS1_SeqM'] == 'N':
            ind = 0
        elif case_comparison_records[0]['PS2_SeqM'] == 'N':
            ind = 1
        else: 
            ind = -1
            
        seq1_ID = 'Recipient-PS'+str(ind+1)
        seq2_ID = 'Donor-PS'+str(ind+1)
        seq1 = case_records[ind][seq1_ID.split('-')[0]]
        seq2 = case_records[ind][seq2_ID.split('-')[0]]
        HLAtyping_list = case_records[ind]['HLATyping']
        tplist = HLAtyping_list.split("+")
        HLAtyping = []
        for tp in tplist:
            if tp.find('/') != -1:
                ambTPlist = tp.split('/')
                tempABMlist = []
                for amTP in ambTPlist:
                    possTPlist = re.sub('[\[\'\]]', '',amTP) # remove possible characters
                    possTPlist = possTPlist.split(", ")
                    tempABMlist.extend(possTPlist)
                HLAtyping.extend(tempABMlist)
                    
            elif tp.find(',') != -1:
                ambTPlist = re.sub('[\[\'\]]', '',tp)
                ambTPlist = ambTPlist.split(', ')
                HLAtyping.extend(ambTPlist)
            else:
                possTPlist = re.sub('[\[\'\]]', '',tp) # remove possible characters
                possTPlist = possTPlist.split(",")
                for item in possTPlist:
                    #HLAtyping.extend(possTPlist)
                    HLAtyping.append(item.replace(" ", ""))

        Sequence= {seq1_ID: seq1, seq2_ID:seq2}

        params = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
        print('CaseID: '+ caseID + ' ; Locus: ' + locus + '; \nParamters:\n')
        print(params)
        alignment, annotation = CompareSeq.compare_DQB1_Targeted_Region(Sequence, params)
        
        if any("Exon" in s for s in alignment.keys()):
            ## save results - for multiple Exon sequences
            for itemID, itemDict in alignment.items():
                
                saveOBJ = {'seq': Sequence, 'params': params, 'alignment':itemDict, 'MMannotation': annotation[itemID], 'SameSeqs': annotation[itemID]['SameSeqs']}
                
                if annotation[itemID]['SameSeqs']: # same seqs
                    Output_fname = singleMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_SameSeqs'
                else:
                    Output_fname = singleMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS'+str(ind+1)+'_'+itemID+'_MisMatchSeqs'
                IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
        else: 
        ## save results -- for one single sequence
            saveOBJ = {'seq': Sequence, 'params': params, 'alignment':alignment, 'MMannotation': annotation}
        
            Output_fname = singleMM_output + 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation'
            IMGTdbIO.save_dict2pickle(saveOBJ, Output_fname)
        
con.close()
