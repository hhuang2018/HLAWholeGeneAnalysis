#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 22:12:52 2017

@author: hhuang2
"""

from utils import IMGTdbIO#, CompareSeq
from collections import Counter

groupType = 'All_paired' # groupType = 'ClassI_paired' # groupType = 'All_paired'

All_loci = ['A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1']
Five_loci = ['A', 'B', 'C', 'DRB1', 'DQB1']
ClassI_loci = ['A', 'B', 'C']
ClassII_loci = ['DRB1', 'DQB1', 'DPB1']

Group_fname = '../Output/Stats/ClassI_Stats_1003_' + groupType + '.pkl'
Stats_Dict = IMGTdbIO.load_pickle2dict(Group_fname)

CaseStats = Stats_Dict['CaseStats']
LocusStats = Stats_Dict['LocusStats']

## TODO1: paired cases HLA typing stats
fname = '../Output/SG39_DRpairs/SG39_pairedCases_Stats.pkl'
Matching_cases_stats = IMGTdbIO.load_pickle2dict(fname)

AlleleStats = {}
for locus in All_loci:
    AlleleStats[locus] = {}
    DRpaired_file = '../Output/SG39_DRpairs/SG39_HLA_'+locus+'_paired.pkl'
    
    DRpaired_table = IMGTdbIO.load_pickle2dict(DRpaired_file)
    num_total = len(DRpaired_table)
    
    for key, item in DRpaired_table.items():
        for ps in ['PS1', 'PS2']:
            typinglist = item[ps]['HLATyping'].split('+')
            if len(typinglist) > 1:
                if 'mismatch' not in AlleleStats[locus].keys():
                    AlleleStats[locus]['mismatch'] = {'count': 1, 'cases': [key]}
                else:
                    AlleleStats[locus]['mismatch']['count'] += 1
                    AlleleStats[locus]['mismatch']['cases'].append(key)
            for tp in typinglist:
                if tp not in AlleleStats[locus].keys():
                    AlleleStats[locus][tp] = {'count': 1, 'cases': [key]}
                else:
                    if key not in AlleleStats[locus][tp]['cases']:
                        AlleleStats[locus][tp]['count'] += 1
                        AlleleStats[locus][tp]['cases'].append(key)
    BothMismatch = AlleleStats[locus]['mismatch']['count'] - len(list(set(AlleleStats[locus]['mismatch']['cases'])))
    SingleMismatch = len(list(set(AlleleStats[locus]['mismatch']['cases']))) - BothMismatch
    
    print('Total number of paired cases at Locus '+locus+': '+ str(num_total))
    print('# One mismached cases:'+ str(SingleMismatch))
    print('# Both mismached cases: ' + str(BothMismatch))
    print('-'*30)

#################################
AlleleStats_fiveLociPaired = {} 
AlleleStats_AllLociPaired = {}
All_caseIDs = []
for locus in All_loci:
    AlleleStats_AllLociPaired[locus] = {}
    AlleleStats_fiveLociPaired[locus] = {}
    DRpaired_file = '../Output/SG39_DRpairs/SG39_HLA_'+locus+'_paired.pkl'
    
    DRpaired_table = IMGTdbIO.load_pickle2dict(DRpaired_file)
    #num_total = len(DRpaired_table)
    AllLociPaired_counter = len(Matching_cases_stats['All_paired'])
    fiveLociPaired_counter = len(Matching_cases_stats['fiveLoci_paired'])
    All_caseIDs.extend(list(DRpaired_table.keys()))
    for key, item in DRpaired_table.items():
        if key in Matching_cases_stats['All_paired']:
            #AllLociPaired_counter += 1
            for ps in ['PS1', 'PS2']:
                typinglist = item[ps]['HLATyping'].split('+')
                if len(typinglist) > 1:
                    if 'mismatch' not in AlleleStats_AllLociPaired[locus].keys():
                        AlleleStats_AllLociPaired[locus]['mismatch'] = {'count': 1, 'cases': [key]}
                    else:
                        AlleleStats_AllLociPaired[locus]['mismatch']['count'] += 1
                        AlleleStats_AllLociPaired[locus]['mismatch']['cases'].append(key)
                for tp in typinglist:
                    if tp not in AlleleStats_AllLociPaired[locus].keys():
                        AlleleStats_AllLociPaired[locus][tp] = {'count': 1, 'cases': [key]}
                    else:
                        if key not in AlleleStats_AllLociPaired[locus][tp]['cases']:
                            AlleleStats_AllLociPaired[locus][tp]['count'] += 1
                            AlleleStats_AllLociPaired[locus][tp]['cases'].append(key)
        if key in Matching_cases_stats['fiveLoci_paired']:
            #fiveLociPaired_counter += 1
            for ps in ['PS1', 'PS2']:
                typinglist = item[ps]['HLATyping'].split('+')
                if len(typinglist) > 1:
                    if 'mismatch' not in AlleleStats_fiveLociPaired[locus].keys():
                        AlleleStats_fiveLociPaired[locus]['mismatch'] = {'count': 1, 'cases': [key]}
                    else:
                        AlleleStats_fiveLociPaired[locus]['mismatch']['count'] += 1
                        AlleleStats_fiveLociPaired[locus]['mismatch']['cases'].append(key)
                for tp in typinglist:
                    if tp not in AlleleStats_fiveLociPaired[locus].keys():
                        AlleleStats_fiveLociPaired[locus][tp] = {'count': 1, 'cases': [key]}
                    else:
                        if key not in AlleleStats_fiveLociPaired[locus][tp]['cases']:
                            AlleleStats_fiveLociPaired[locus][tp]['count'] += 1
                            AlleleStats_fiveLociPaired[locus][tp]['cases'].append(key)
            
            
    BothMismatch_AllLociPaired = AlleleStats_AllLociPaired[locus]['mismatch']['count'] - len(list(set(AlleleStats_AllLociPaired[locus]['mismatch']['cases'])))
    SingleMismatch_AllLociPaired = len(list(set(AlleleStats_AllLociPaired[locus]['mismatch']['cases']))) - BothMismatch_AllLociPaired
    print('>> All 6 loci paired cases:')
    print('Total number of paired cases at Locus '+locus+': '+ str(AllLociPaired_counter))
    print('# One mismached cases:'+ str(SingleMismatch_AllLociPaired))
    print('# Both mismached cases: ' + str(BothMismatch_AllLociPaired))
    print('-'*30)
    
    BothMismatch_fiveLociPaired = AlleleStats_fiveLociPaired[locus]['mismatch']['count'] - len(list(set(AlleleStats_fiveLociPaired[locus]['mismatch']['cases'])))
    SingleMismatch_fiveLociPaired = len(list(set(AlleleStats_fiveLociPaired[locus]['mismatch']['cases']))) - BothMismatch_fiveLociPaired
    print('>> Five loci paired cases:')
    print('Total number of paired cases at Locus '+locus+': '+ str(fiveLociPaired_counter))
    print('# One mismached cases:'+ str(SingleMismatch_fiveLociPaired))
    print('# Both mismached cases: ' + str(BothMismatch_fiveLociPaired))
    print('-'*30)

AllLoci_case_freq = {}
fiveLoci_case_freq = {}
for locus in All_loci:
    AllLoci_case_freq[locus] = Counter(AlleleStats_AllLociPaired[locus]['mismatch']['cases'])
    fiveLoci_case_freq[locus] = Counter(AlleleStats_fiveLociPaired[locus]['mismatch']['cases'])

All_caseIDs = list(set(All_caseIDs))
MatchingType_AllLoci = {}
MatchingType_fiveLoci = {}
for key in All_caseIDs:
    if key in Matching_cases_stats['All_paired']:
        Match_allele_count_AllLoci = 12
        for locus in All_loci:
            if key in AllLoci_case_freq[locus].keys():
                Match_allele_count_AllLoci -= AllLoci_case_freq[locus][key]
        if str(Match_allele_count_AllLoci) in MatchingType_AllLoci.keys():
            MatchingType_AllLoci[str(Match_allele_count_AllLoci)]['count'] += 1
            MatchingType_AllLoci[str(Match_allele_count_AllLoci)]['cases'].append(key)
        else:
            MatchingType_AllLoci[str(Match_allele_count_AllLoci)] = {'count': 1, 'cases': [key]}
            
    if key in Matching_cases_stats['fiveLoci_paired']:
        Match_allele_count_fiveLoci = 10
        for locus in Five_loci:
            if key in fiveLoci_case_freq[locus].keys():
                Match_allele_count_fiveLoci -= fiveLoci_case_freq[locus][key]
        if str(Match_allele_count_fiveLoci) in MatchingType_fiveLoci.keys():
            MatchingType_fiveLoci[str(Match_allele_count_fiveLoci)]['count'] += 1
            MatchingType_fiveLoci[str(Match_allele_count_fiveLoci)]['cases'].append(key)
        else:
            MatchingType_fiveLoci[str(Match_allele_count_fiveLoci)] = {'count': 1, 'cases': [key]}        
print('-*'*30)    
print('All 6 loci paired cases:')
for key, item in MatchingType_AllLoci.items():
    print(key+'/12 matched cases: '+str(MatchingType_AllLoci[key]['count']))      
print('-*'*30)          
print('Five loci paired cases:')
for key, item in MatchingType_fiveLoci.items():
    print(key+'/10 matched cases: '+str(MatchingType_fiveLoci[key]['count']))      
print('-*'*30)     
                
## TODO2: Count of each mismatch types -- by Alelle; by mismatch region; by Specific mismatch
groupType = 'All_paired' # groupType = 'ClassI_paired' # groupType = 'All_paired'

Group_fname = '../Output/Stats/ClassII_Stats_1003_' + groupType + '.pkl'
Stats_Dict = IMGTdbIO.load_pickle2dict(Group_fname)

CaseStats = Stats_Dict['CaseStats']
LocusStats = Stats_Dict['LocusStats']
MisMatch_region = {}
for key, item in LocusStats.items():
    locus = key.split('*')[0]
    for itemKey, itemItem in item.items():
        keyWord = key + '_' + itemKey
        if keyWord not in MisMatch_region.keys():
            MisMatch_region[keyWord] = 1
        else:
            MisMatch_region[keyWord] += 1