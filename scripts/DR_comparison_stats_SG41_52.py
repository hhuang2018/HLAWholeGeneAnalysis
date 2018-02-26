#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:42:06 2017

@author: hhuang2
"""

import glob
import sqlite3 as sql
# from utils import phase_block_check as ps
from utils import IMGTdbIO, CompareSeq
import os
import re


fname = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/SG41_52_pairedCases_Stats.pkl'
Matching_cases_stats = IMGTdbIO.load_pickle2dict(fname)

## 'All_paired'
groupType = 'fiveLoci_paired' # groupType = 'ClassI_paired' # groupType = 'All_paired'
group_caseIDs = Matching_cases_stats[groupType]
All_loci = ['A', 'B', 'C', 'DRB1', 'DQB1']#, 'DPB1']
ClassI_loci = ['A', 'B', 'C']
ClassII_loci = ['DRB1', 'DQB1']

CaseStats = {}
LocusStats = {}
#MatchStats = {}
for caseID in group_caseIDs:
    # 
    for locus in ClassI_loci:
        ARSregion = ['Exon2', 'Exon3']
        bothMM_output = "../Output/SG41_52/2018/IMGTv3310/SG41_52_bothMisMatched_locus_" + locus + "_0125_TargetedAlignment/" # "_1218_TargetedAlignment/"
        
        singleMM_output = "../Output/SG41_52/2018/IMGTv3310/SG41_52_singleMisMatched_" + locus + "_0125_TargetedAlignment/"
            
        ### Cases where both sequences don't match
        if caseID in Matching_cases_stats[locus+'_both_Seqmm']:
            mm_file_PS1 = bothMM_output+ 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS1.pkl'
            mm_file_PS2 = bothMM_output+ 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS2.pkl'
            
            mm_locus_stats_PS1 = IMGTdbIO.load_pickle2dict(mm_file_PS1)
            mm_locus_stats_PS1 = CompareSeq.rmRefAln(mm_locus_stats_PS1)
            mm_locus_stats_PS2 = IMGTdbIO.load_pickle2dict(mm_file_PS2)
            mm_locus_stats_PS2 = CompareSeq.rmRefAln(mm_locus_stats_PS2)
            
            #if len(mm_locus_stats_PS1['MMpos']) > 20 or len(mm_locus_stats_PS2['MMpos']) > 20:
            if CompareSeq.isARSmm(mm_locus_stats_PS1['MMannotation'].values(), ARSregion) and CompareSeq.isARSmm(mm_locus_stats_PS2['MMannotation'].values(), ARSregion):
                # probably phase set swap.
                seq_ps1 = mm_locus_stats_PS1['seq']
                seq_ps2 = mm_locus_stats_PS2['seq']
                params_ps1 = mm_locus_stats_PS1['params']
                params_ps2 = mm_locus_stats_PS2['params']
                #tp = params_ps2['HLAtyping']
                #tp = [tp[1], tp[0]]
                #params_ps2['HLAtyping'] = tp
                swapped_alignment = CompareSeq.swapPS_comparison(seq_ps1, params_ps1, seq_ps2, params_ps2, caseID)
                
                #if max([len(swapped_alignment['PS1']['MMpos']), len(swapped_alignment['PS2']['MMpos'])]) < max([len(mm_locus_stats_PS1['MMpos']), len(mm_locus_stats_PS2['MMpos'])]):
                if not CompareSeq.isARSmm(swapped_alignment['PS1']['MMannotation'].values(), ARSregion) or not CompareSeq.isARSmm(swapped_alignment['PS2']['MMannotation'].values(), ARSregion):
                    
                    # if swapped case is better, then use the swapped case
                    mm_locus_stats_PS1 = swapped_alignment['PS1']
                    mm_locus_stats_PS2 = swapped_alignment['PS2']
            
            params_ps1 = mm_locus_stats_PS1['params']
            params_ps2 = mm_locus_stats_PS2['params']
            # caseStats
            if caseID in CaseStats.keys():
                CaseStats[caseID][locus] = {'PS1': CompareSeq.RegionCount(mm_locus_stats_PS1['MMannotation'], locus, True), 'PS2': CompareSeq.RegionCount(mm_locus_stats_PS2['MMannotation'], locus, True)}
                CaseStats[caseID][locus]['HLAtyping'] = params_ps1['HLAtyping'] + params_ps2['HLAtyping']

            else:
                CaseStats[caseID] = {locus:{'PS1': CompareSeq.RegionCount(mm_locus_stats_PS1['MMannotation'], locus, True), 'PS2': CompareSeq.RegionCount(mm_locus_stats_PS2['MMannotation'], locus, True), 
                         'HLAtyping':params_ps1['HLAtyping'] + params_ps2['HLAtyping']}}
            
            if len(mm_locus_stats_PS1['params']['HLAtyping']) == 1: 
                typing = mm_locus_stats_PS1['params']['HLAtyping'][0]
                if typing in LocusStats.keys():
                    if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                        if 'ARS' in LocusStats[typing].keys():
                            LocusStats[typing]['ARS'].append(caseID) 
                        else:
                            LocusStats[typing] = {'ARS': [caseID]}
                        #LocusStats[typing]['ARS'].append(caseID) 
                    if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                        if 'Non_ARS_exon' in LocusStats[typing].keys():
                            LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        else:
                            LocusStats[typing] = {'Non_ARS_exon': [caseID]}
                        #LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                    if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                        if 'Intron' in LocusStats[typing].keys():
                            LocusStats[typing]['Intron'].append(caseID) 
                        else:
                            LocusStats[typing] = {'Intron': [caseID]}
                        #LocusStats[typing]['Intron'].append(caseID)
                else:
                    LocusStats[typing] = {'ARS': [], 'Non_ARS_exon': [], 'Intron': []}
                    if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                        LocusStats[typing]['ARS'].append(caseID) 
                    if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                        LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                    if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                        LocusStats[typing]['Intron'].append(caseID)
                        
                for key, item in mm_locus_stats_PS1['MMannotation'].items():
                    if key.isdigit():
                        if item in LocusStats[typing].keys():
                            LocusStats[typing][item].append(caseID)
                        else:
                            LocusStats[typing] = {item: [caseID]}
                                    
            if len(mm_locus_stats_PS2['params']['HLAtyping']) == 1: 
                typing = mm_locus_stats_PS2['params']['HLAtyping'][0]
                if typing in LocusStats.keys():
                    if CaseStats[caseID][locus]['PS2']['ARS'] > 0:
                        if 'ARS' in LocusStats[typing].keys():
                            LocusStats[typing]['ARS'].append(caseID) 
                        else:
                            LocusStats[typing] = {'ARS': [caseID]}
                        #LocusStats[typing]['ARS'].append(caseID) 
                    if CaseStats[caseID][locus]['PS2']['Non_ARS_exon'] >0:
                        if 'Non_ARS_exon' in LocusStats[typing].keys():
                            LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        else:
                            LocusStats[typing] = {'Non_ARS_exon': [caseID]}
                        #LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                    if CaseStats[caseID][locus]['PS2']['Intron'] > 0:
                        if 'Intron' in LocusStats[typing].keys():
                            LocusStats[typing]['Intron'].append(caseID) 
                        else:
                            LocusStats[typing] = {'Intron': [caseID]}
                        #LocusStats[typing]['Intron'].append(caseID)
                else:
                    LocusStats[typing] = {'ARS': [], 'Non_ARS_exon': [], 'Intron': []}
                    if CaseStats[caseID][locus]['PS2']['ARS'] > 0:
                        LocusStats[typing]['ARS'].append(caseID) 
                    if CaseStats[caseID][locus]['PS2']['Non_ARS_exon'] >0:
                        LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                    if CaseStats[caseID][locus]['PS2']['Intron'] > 0:
                        LocusStats[typing]['Intron'].append(caseID)
                
                for key, item in mm_locus_stats_PS2['MMannotation'].items():
                    if key.isdigit():
                        if item in LocusStats[typing].keys():
                            LocusStats[typing][item].append(caseID)
                        else:
                            LocusStats[typing] = {item: [caseID]}
                            
        elif caseID in Matching_cases_stats[locus+'_one_Seqmm']:
        ### Cases where only one sequence doesn't match
            mm_file = singleMM_output+ 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation.pkl'
            
            mm_locus_stats = IMGTdbIO.load_pickle2dict(mm_file)
            mm_locus_stats = CompareSeq.rmRefAln(mm_locus_stats)
            
            params_singmm = mm_locus_stats['params']
            # caseStats
            if caseID in CaseStats.keys():
                CaseStats[caseID][locus] = {'PS1': CompareSeq.RegionCount(mm_locus_stats['MMannotation'], locus, True)}
                CaseStats[caseID][locus]['HLAtyping'] = params_singmm['HLAtyping'] 
                
            else:
                CaseStats[caseID] = {locus:{'PS1': CompareSeq.RegionCount(mm_locus_stats['MMannotation'], locus, True), 
                         'HLAtyping':params_singmm['HLAtyping']}}
            
            if len(params_singmm['HLAtyping']) == 1: 
                typing = params_singmm['HLAtyping'][0]
                if typing in LocusStats.keys():
                    if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                        if 'ARS' in LocusStats[typing].keys():
                            LocusStats[typing]['ARS'].append(caseID) 
                        else:
                            LocusStats[typing] = {'ARS': [caseID]}
                    if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                        if 'Non_ARS_exon' in LocusStats[typing].keys():
                            LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        else:
                            LocusStats[typing] = {'Non_ARS_exon': [caseID]}
                         
                    if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                        if 'Intron' in LocusStats[typing].keys():
                            LocusStats[typing]['Intron'].append(caseID) 
                        else:
                            LocusStats[typing] = {'Intron': [caseID]}
                else:
                    LocusStats[typing] = {'ARS': [], 'Non_ARS_exon': [], 'Intron': []}
                    if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                        LocusStats[typing]['ARS'].append(caseID) 
                    if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                        LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                    if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                        LocusStats[typing]['Intron'].append(caseID)
               
                for key, item in mm_locus_stats['MMannotation'].items():
                    if key.isdigit():
                        if item in LocusStats[typing].keys():
                            LocusStats[typing][item].append(caseID)
                        else:
                            LocusStats[typing] = {item: [caseID]}

ClassI_stats = {'CaseStats': CaseStats, 'LocusStats': LocusStats}
IMGTdbIO.save_dict2pickle(ClassI_stats, '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/ClassI_Stats_0125_'+groupType) #1220_'+groupType)

# Class II
#Group_fname = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/ClassI_Stats_0125_' + groupType + '.pkl'
#Stats_Dict = IMGTdbIO.load_pickle2dict(Group_fname)

#CaseStats = Stats_Dict['CaseStats']
#LocusStats = Stats_Dict['LocusStats']
for caseID in group_caseIDs:
    # 
    for locus in ClassII_loci:
        bothMM_output = "../Output/SG41_52/2018/IMGTv3310/SG41_52_bothMisMatched_locus_" + locus + "_0125_TargetedAlignment/"
        
        singleMM_output = "../Output/SG41_52/2018/IMGTv3310/SG41_52_singleMisMatched_" + locus + "_0125_TargetedAlignment/"
            
        ### Cases where both sequences don't match
        if caseID in Matching_cases_stats[locus+'_both_Seqmm']:
            mm_file_PS1 = glob.glob(bothMM_output+ 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS1*.pkl')

            for file_id in mm_file_PS1:
                mm_locus_stats_PS1 = IMGTdbIO.load_pickle2dict(file_id)
                mm_locus_stats_PS1 = CompareSeq.rmRefAln(mm_locus_stats_PS1)
                #if mm_locus_stats_PS1['SameSeqs']: # if the exons are the same
                if 'Exon' in file_id:
                    for key, item in mm_locus_stats_PS1['MMannotation'].items():
                        if key.isdigit():
                            tempItem = item.split('.')
                            if int(tempItem[1]) <0 and int(tempItem[0][-1])== int(file_id.split('_')[-2][-1]):
                                tempItem[0] = 'Intron'+ str(int(file_id.split('_')[-2][-1])-1)
                            mm_locus_stats_PS1['MMannotation'][key] = '.'.join(tempItem)
                
                if caseID in CaseStats.keys():
                    if locus not in CaseStats[caseID].keys():
                        CaseStats[caseID][locus] = {}
                    if 'PS1' not in CaseStats[caseID][locus].keys():
                        CaseStats[caseID][locus]['PS1'] = CompareSeq.RegionCount(mm_locus_stats_PS1['MMannotation'], locus)
                    else:
                        tempStats = CompareSeq.RegionCount(mm_locus_stats_PS1['MMannotation'], locus)
                        for key, item in tempStats.items():
                            CaseStats[caseID][locus]['PS1'][key] += item
                        
                    if 'HLAtyping' not in CaseStats[caseID][locus].keys():
                        CaseStats[caseID][locus]['HLAtyping'] = mm_locus_stats_PS1['params']['HLAtyping']
    
                else:
                    CaseStats[caseID] = {locus:{'PS1': CompareSeq.RegionCount(mm_locus_stats_PS1['MMannotation'], locus), 
                             'HLAtyping':mm_locus_stats_PS1['params']['HLAtyping']}}
       
                if len(mm_locus_stats_PS1['params']['HLAtyping']) == 1: 
                    typing = mm_locus_stats_PS1['params']['HLAtyping'][0]
                    if typing in LocusStats.keys():
                        if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                            if 'ARS' in LocusStats[typing].keys():
                                LocusStats[typing]['ARS'].append(caseID) 
                            else:
                                LocusStats[typing] = {'ARS': [caseID]}
                            #LocusStats[typing]['ARS'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                            if 'Non_ARS_exon' in LocusStats[typing].keys():
                                LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                            else:
                                LocusStats[typing] = {'Non_ARS_exon': [caseID]}
                            #LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                            if 'Intron' in LocusStats[typing].keys():
                                LocusStats[typing]['Intron'].append(caseID) 
                            else:
                                LocusStats[typing] = {'Intron': [caseID]}
                            #LocusStats[typing]['Intron'].append(caseID)
                    else:
                        LocusStats[typing] = {'ARS': [], 'Non_ARS_exon': [], 'Intron': []}
                        if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                            LocusStats[typing]['ARS'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                            LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                            LocusStats[typing]['Intron'].append(caseID)
                            
                    for key, item in mm_locus_stats_PS1['MMannotation'].items():
                        if key.isdigit():
                            if item in LocusStats[typing].keys():
                                LocusStats[typing][item].append(caseID)
                            else:
                                LocusStats[typing] = {item: [caseID]}
            
            mm_file_PS2 = glob.glob(bothMM_output+ 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation_PS2*.pkl')

            for file_id in mm_file_PS2:
                mm_locus_stats_PS2 = IMGTdbIO.load_pickle2dict(file_id)
                mm_locus_stats_PS2 = CompareSeq.rmRefAln(mm_locus_stats_PS2)
                if 'Exon' in file_id:
                    for key, item in mm_locus_stats_PS2['MMannotation'].items():
                        if key.isdigit():
                            tempItem = item.split('.')
                            if int(tempItem[1]) <0 and int(tempItem[0][-1])== int(file_id.split('_')[-2][-1]):
                                tempItem[0] = 'Intron'+ str(int(file_id.split('_')[-2][-1])-1)
                            mm_locus_stats_PS2['MMannotation'][key] = '.'.join(tempItem)
                
                if caseID in CaseStats.keys():
                    if locus not in CaseStats[caseID].keys():
                        CaseStats[caseID][locus] = {}
                    if 'PS2' not in CaseStats[caseID][locus].keys():
                        CaseStats[caseID][locus]['PS2'] = CompareSeq.RegionCount(mm_locus_stats_PS2['MMannotation'], locus)
                    else:
                        tempStats = CompareSeq.RegionCount(mm_locus_stats_PS2['MMannotation'], locus)
                        for key, item in tempStats.items():
                            CaseStats[caseID][locus]['PS2'][key] += item
                        
                    CaseStats[caseID][locus]['HLAtyping'] += mm_locus_stats_PS2['params']['HLAtyping']
    
                else:
                    CaseStats[caseID] = {locus:{'PS2': CompareSeq.RegionCount(mm_locus_stats_PS2['MMannotation'], locus), 
                             'HLAtyping':mm_locus_stats_PS2['params']['HLAtyping']}}
       
                if len(mm_locus_stats_PS2['params']['HLAtyping']) == 1: 
                    typing = mm_locus_stats_PS2['params']['HLAtyping'][0]
                    if typing in LocusStats.keys():
                        if CaseStats[caseID][locus]['PS2']['ARS'] > 0:
                            if 'ARS' in LocusStats[typing].keys():
                                LocusStats[typing]['ARS'].append(caseID) 
                            else:
                                LocusStats[typing] = {'ARS': [caseID]}
                            #LocusStats[typing]['ARS'].append(caseID) 
                        if CaseStats[caseID][locus]['PS2']['Non_ARS_exon'] >0:
                            if 'Non_ARS_exon' in LocusStats[typing].keys():
                                LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                            else:
                                LocusStats[typing] = {'Non_ARS_exon': [caseID]}
                            #LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        if CaseStats[caseID][locus]['PS2']['Intron'] > 0:
                            if 'Intron' in LocusStats[typing].keys():
                                LocusStats[typing]['Intron'].append(caseID) 
                            else:
                                LocusStats[typing] = {'Intron': [caseID]}
                            #LocusStats[typing]['Intron'].append(caseID)
                    else:
                        LocusStats[typing] = {'ARS': [], 'Non_ARS_exon': [], 'Intron': []}
                        if CaseStats[caseID][locus]['PS2']['ARS'] > 0:
                            LocusStats[typing]['ARS'].append(caseID) 
                        if CaseStats[caseID][locus]['PS2']['Non_ARS_exon'] >0:
                            LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        if CaseStats[caseID][locus]['PS2']['Intron'] > 0:
                            LocusStats[typing]['Intron'].append(caseID)
                            
                    for key, item in mm_locus_stats_PS2['MMannotation'].items():
                        if key.isdigit():
                            if item in LocusStats[typing].keys():
                                LocusStats[typing][item].append(caseID)
                            else:
                                LocusStats[typing] = {item: [caseID]}   

        # Single mismatch cases      
        elif caseID in Matching_cases_stats[locus+'_one_Seqmm']:
        ### Cases where only one sequence doesn't match
            mm_file_PS = glob.glob(singleMM_output+ 'CaseID_'+ caseID + '_Locus_' + locus + '_annotation*.pkl')

            for file_id in mm_file_PS:
                mm_locus_stats_PS = IMGTdbIO.load_pickle2dict(file_id)
                mm_locus_stats_PS = CompareSeq.rmRefAln(mm_locus_stats_PS)
                if 'Exon' in file_id:
                    for key, item in mm_locus_stats_PS['MMannotation'].items():
                        if key.isdigit():
                            tempItem = item.split('.')
                            if int(tempItem[1]) <0 and int(tempItem[0][-1])== int(file_id.split('_')[-2][-1]):
                                tempItem[0] = 'Intron'+ str(int(file_id.split('_')[-2][-1])-1)
                            mm_locus_stats_PS['MMannotation'][key] = '.'.join(tempItem)
                
                if caseID in CaseStats.keys():
                    if locus not in CaseStats[caseID].keys():
                        CaseStats[caseID][locus] = {}
                    if 'PS1' not in CaseStats[caseID][locus].keys():
                        CaseStats[caseID][locus]['PS1'] = CompareSeq.RegionCount(mm_locus_stats_PS['MMannotation'], locus)
                    else:
                        tempStats = CompareSeq.RegionCount(mm_locus_stats_PS['MMannotation'], locus)
                        for key, item in tempStats.items():
                            CaseStats[caseID][locus]['PS1'][key] += item
                        
                    if 'HLAtyping' not in CaseStats[caseID][locus].keys():
                        CaseStats[caseID][locus]['HLAtyping'] = mm_locus_stats_PS['params']['HLAtyping']
    
                else:
                    CaseStats[caseID] = {locus:{'PS1': CompareSeq.RegionCount(mm_locus_stats_PS['MMannotation'], locus), 
                             'HLAtyping':mm_locus_stats_PS['params']['HLAtyping']}}
       
                if len(mm_locus_stats_PS['params']['HLAtyping']) == 1: 
                    typing = mm_locus_stats_PS['params']['HLAtyping'][0]
                    if typing in LocusStats.keys():
                        if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                            if 'ARS' in LocusStats[typing].keys():
                                LocusStats[typing]['ARS'].append(caseID) 
                            else:
                                LocusStats[typing] = {'ARS': [caseID]}
                            #LocusStats[typing]['ARS'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                            if 'Non_ARS_exon' in LocusStats[typing].keys():
                                LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                            else:
                                LocusStats[typing] = {'Non_ARS_exon': [caseID]}
                            #LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                            if 'Intron' in LocusStats[typing].keys():
                                LocusStats[typing]['Intron'].append(caseID) 
                            else:
                                LocusStats[typing] = {'Intron': [caseID]}
                            #LocusStats[typing]['Intron'].append(caseID)
                    else:
                        LocusStats[typing] = {'ARS': [], 'Non_ARS_exon': [], 'Intron': []}
                        if CaseStats[caseID][locus]['PS1']['ARS'] > 0:
                            LocusStats[typing]['ARS'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Non_ARS_exon'] >0:
                            LocusStats[typing]['Non_ARS_exon'].append(caseID) 
                        if CaseStats[caseID][locus]['PS1']['Intron'] > 0:
                            LocusStats[typing]['Intron'].append(caseID)
                            
                    for key, item in mm_locus_stats_PS['MMannotation'].items():
                        if key.isdigit():
                            if item in LocusStats[typing].keys():
                                LocusStats[typing][item].append(caseID)
                            else:
                                LocusStats[typing] = {item: [caseID]}    
        
ClassII_stats = {'CaseStats': CaseStats, 'LocusStats': LocusStats}
IMGTdbIO.save_dict2pickle(ClassII_stats, '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/ClassII_Stats_0125_'+groupType)


#ClassI_stats = IMGTdbIO.load_pickle2dict('../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/ClassI_Stats_0125_fiveLoci_paired.pkl')
#ClassII_stats = IMGTdbIO.load_pickle2dict('../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/ClassII_Stats_0125_fiveLoci_paired.pkl')

fiveLociPaired_stats = {'CaseStats': CaseStats, 'LocusStats': LocusStats}
IMGTdbIO.save_dict2pickle(fiveLociPaired_stats, '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/fiveLoci_paired_Stats_0125_'+groupType)



####### Swapped cases for DQB1
## Swapped case check
            if 'PS1' in CaseStats[caseID][locus].keys() and 'PS2' in CaseStats[caseID][locus].keys():
                if CaseStats[caseID][locus]['PS1']['ARS'] > 5 and CaseStats[caseID][locus]['PS2']['ARS'] > 5:
                    DB_fp = '../Output/SG39_DRpairs/SG39_HLA_'+ locus +'_paired.db'
                    con = sql.connect(DB_fp)
                    con.row_factory = sql.Row
                    cur = con.cursor()
                    t = (caseID,)
                    cur.execute('SELECT * FROM OriginalSeqs WHERE BMT_caseID = ?', t)
                    case_records = cur.fetchall()
                    Sequence = {}
                    Params = {}
                    for ind in range(2):
                        
                        seq1_ID = 'Recipient-PS'+str(ind+1)
                        seq2_ID = 'Donor-PS'+str(ind+1)
                        seq1 = case_records[ind][seq1_ID.split('-')[0]]
                        seq2 = case_records[ind][seq2_ID.split('-')[0]]
                        HLAtyping_list = case_records[ind]['HLATyping']
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
                        
                        Sequence[str(ind)]= {seq1_ID: seq1, seq2_ID:seq2}
                        if ind == 0:
                            algn_file = mm_locus_stats_PS1['params']['algn_file']
                        else:
                            algn_file = mm_locus_stats_PS2['params']['algn_file']
                            
                        Params[str(ind)] = {'algn_file': algn_file, 'saveFile': True, 'HLAtyping': HLAtyping}
                        
                    swapped_alignment = CompareSeq.swapPS_comparison(Sequence['0'], Params['0'], Sequence['1'], Params['1'], caseID)
                    
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
                        