#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 10:38:30 2017

@author: hhuang2
"""

import glob
import csv

from utils import IMGTdbIO, CompareSeq


groupType = 'fiveLoci_paired' # groupType = 'ClassI_paired' # groupType = 'All_paired' ; 'fiveLoci_paired'

All_loci = ['A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1']
Five_loci = ['A', 'B', 'C', 'DRB1', 'DQB1']
ClassI_loci = ['A', 'B', 'C']
ClassII_loci = ['DRB1', 'DQB1', 'DPB1']

Group_fname = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/fiveLoci_paired_Stats_0125_' + groupType + '.pkl'
Stats_Dict = IMGTdbIO.load_pickle2dict(Group_fname)

CaseStats = Stats_Dict['CaseStats']
LocusStats = Stats_Dict['LocusStats']

db_fp= '../Database/'
#key = '83687'
#CaseStats[key] 
#key in group_caseIDs
## : paired cases HLA typing stats
fname = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/SG41_52_pairedCases_Stats.pkl'
Matching_cases_stats = IMGTdbIO.load_pickle2dict(fname)

group_caseIDs = Matching_cases_stats[groupType]

CaseMatchTable = {}
for locus in Five_loci:
    
    singleMismatch_fp = '../Output/SG41_52/2018/IMGTv3310/SG41_52_singleMisMatched_'+locus+'_0125_TargetedAlignment/'
    bothMisMatch_fp = '../Output/SG41_52/2018/IMGTv3310/SG41_52_bothMisMatched_locus_'+locus+'_0125_TargetedAlignment/'
    
    if locus in ClassI_loci:
        ARS_exon = ['Exon2', 'Exon3']
    elif locus in ClassII_loci:
        ARS_exon = ['Exon2']
    #CaseMatchTable[locus] = {}
    DRpaired_file = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpairs/SG41_52_HLA_'+locus+'_wComparison.pkl'
    
    DRpaired_table = IMGTdbIO.load_pickle2dict(DRpaired_file)

    num_total = len(DRpaired_table)
    
    #for key, item in DRpaired_table.items():
    for key in group_caseIDs:
        item = DRpaired_table[key]
        
        if key in Matching_cases_stats[locus+'_both_SeqMatch']: # both sequences are matched
            if item['Audit'] == 'Y' and item['Active'] == 'Y':
                if key not in CaseMatchTable.keys():
                    CaseMatchTable[key] = {}
                
                if locus not in CaseMatchTable[key].keys():
                    CaseMatchTable[key][locus] = {}
                    
                for ps in ['PS1', 'PS2']:
                    typinglist = item[ps]['HLATyping'].split('+')
                    twoField = IMGTdbIO.Full2TwoField(typinglist)
                    if item[ps]['isSeqMatch'] == 'Y':
                        CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': True,
                                      'ARS': [], 'non_ARS_exon': [], 'Intron': [], 'UTR': []}
                    else:
                        CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': False,
                                      'ARS': [], 'non_ARS_exon': [], 'Intron': [], 'UTR': []}
                        print('Case ID '+key+ ' Locus '+locus+' '+ ps + ' Sequences don\'t match. Please double check!')
                    
        elif key in Matching_cases_stats[locus+'_one_Seqmm']: # Only one sequence is matched
            if item['Audit'] == 'Y' and item['Active'] == 'Y':
                if key not in CaseMatchTable.keys():
                    CaseMatchTable[key] = {}
                if locus not in CaseMatchTable[key].keys():
                    CaseMatchTable[key][locus] = {}
                for ps in ['PS1', 'PS2']:
                    typinglist = item[ps]['HLATyping'].split('+')
                    twoField = IMGTdbIO.Full2TwoField(typinglist)
                    if item[ps]['isSeqMatch'] == 'Y':
                        CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': True,
                                      'ARS': [], 'non_ARS_exon': [], 'Intron': [], 'UTR': []}
                    else:
                        '''if locus in ClassI_loci:
                            IntronKeys = []
                            ARSExonKeys =[]
                            NonARSExonKeys = []
                            UTRKeys =[]
                            for annKey, annItem in CaseStats[key][locus][ps]['MMannotation'].items():
                                if annKey.isdigit(): 
                                    #annReads = IMGTdbIO.annotationFormat(annKey, annItem, seqAlgn_stats['alignment']) #annItem+'[D:'+seqAlgn_stats['alignment']['Donor-'+ps]
                                    annReads = annItem + '>'+ 'single_mm_Case-' +key+'-'+locus 
                                    if 'Intron' in annItem.split('.')[0]:
                                        IntronKeys.append(annReads)
                                    elif 'UTR' in annItem.split('.')[0]:
                                        UTRKeys.append(annReads)
                                    elif 'Exon' in annItem.split('.')[0]:
                                        if annItem.split('.')[0] in ARS_exon:
                                            ARSExonKeys.append(annReads)
                                        else:
                                            NonARSExonKeys.append(annReads) 
                            CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': False,
                                          'ARS': ARSExonKeys, 
                                          'non_ARS_exon': NonARSExonKeys, 
                                          'Intron': IntronKeys, 
                                          'UTR': UTRKeys}
                        else:'''
                        CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': False,
                                      'ARS': [], 'non_ARS_exon': [], 'Intron': [], 'UTR': []}
                       
                        
                        AlgnFiles = glob.glob(singleMismatch_fp+'CaseID_'+key+'_Locus_'+locus+'*.pkl')
                        
                        for algnFile in AlgnFiles:
                            seqAlgn_stats = IMGTdbIO.load_pickle2dict(algnFile)
                            seqAlgn_stats = CompareSeq.rmRefAln(seqAlgn_stats)
                            
                            for annKey, annItem in seqAlgn_stats['MMannotation'].items():
                                if annKey.isdigit(): 
                                    print(annKey)
                                    annReads = IMGTdbIO.annotationFormat(annKey, annItem, seqAlgn_stats['alignment'])
                                    if 'Intron' in annItem.split('.')[0]:
                                        CaseMatchTable[key][locus][ps]['Intron'].append(annReads)
                                    elif 'UTR' in annItem.split('.')[0]:
                                        CaseMatchTable[key][locus][ps]['UTR'].append(annReads)
                                    elif 'Exon' in annItem.split('.')[0]:
                                        if annItem.split('.')[0] in ARS_exon:
                                            CaseMatchTable[key][locus][ps]['ARS'].append(annReads)
                                        else:
                                            CaseMatchTable[key][locus][ps]['non_ARS_exon'].append(annReads)

        elif key in Matching_cases_stats[locus+'_both_Seqmm']: # both sequence are mismatched
            # remove Non-audited cases 71110, 72169, 76358
            if item['Audit'] == 'Y' and item['Active'] == 'Y':
                if key not in CaseMatchTable.keys():
                    CaseMatchTable[key] = {}
                if locus not in CaseMatchTable[key].keys():
                    CaseMatchTable[key][locus] = {}
                for ps in ['PS1', 'PS2']:
                    #if ps == 'PS1':
                    ##    typinglist = [item['PS1']['HLATyping'].split('+')[0],item['PS2']['HLATyping'].split('+')[1]]
                    #else:
                    #    typinglist = [item['PS1']['HLATyping'].split('+')[1],item['PS2']['HLATyping'].split('+')[0]]
                    
                    typinglist = item[ps]['HLATyping'].split('+')
                    twoField = IMGTdbIO.Full2TwoField(typinglist)
                    if item[ps]['isSeqMatch'] == 'Y':
                        CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': True,
                                      'ARS': [], 'non_ARS_exon': [], 'Intron': [], 'UTR': []}
                    else:
                        if locus in ClassI_loci: 
                            IntronKeys = []
                            ARSExonKeys =[]
                            NonARSExonKeys = []
                            UTRKeys =[]
                            for annKey, annItem in CaseStats[key][locus][ps]['MMannotation'].items():
                                if annKey.isdigit(): 
                                    #annReads = IMGTdbIO.annotationFormat(annKey, annItem, seqAlgn_stats['alignment']) #annItem+'[D:'+seqAlgn_stats['alignment']['Donor-'+ps]
                                    annReads = annItem + '>'+ 'swapped_Case-' +key+'-'+locus 
                                    if 'Intron' in annItem.split('.')[0]:
                                        IntronKeys.append(annReads)
                                    elif 'UTR' in annItem.split('.')[0]:
                                        UTRKeys.append(annReads)
                                    elif 'Exon' in annItem.split('.')[0]:
                                        if annItem.split('.')[0] in ARS_exon:
                                            ARSExonKeys.append(annReads)
                                        else:
                                            NonARSExonKeys.append(annReads) 
                            CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': False,
                                          'ARS': ARSExonKeys, 
                                          'non_ARS_exon': NonARSExonKeys, 
                                          'Intron': IntronKeys, 
                                          'UTR': UTRKeys}
                        else:
                            
                        
                            # for swapped case:
                            #AlgnFiles = glob.glob(bothMisMatch_fp+'CaseID_'+key+'_Locus_'+locus+'*'+ps+'*.pkl')
                            #swapped_tpyinglist = seqAlgn_stats['params']['HLAtyping']
                            #swapped_twoField = IMGTdbIO.Full2TwoField(swapped_tpyinglist)
                            #CaseMatchTable[key][locus][ps]['HLAtyping'] = swapped_tpyinglist
                            #CaseMatchTable[key][locus][ps]['TwoField'] = swapped_twoField
                            CaseMatchTable[key][locus][ps] = {'HLAtyping': typinglist, 'TwoField': twoField, 'SameSeq': False,
                                      'ARS': [], 'non_ARS_exon': [], 'Intron': [], 'UTR': []}
                            
                            AlgnFiles = glob.glob(bothMisMatch_fp+'CaseID_'+key+'_Locus_'+locus+'_annotation_'+ps+'*.pkl')
                            for algnFile in AlgnFiles:
                                seqAlgn_stats = IMGTdbIO.load_pickle2dict(algnFile)
                                seqAlgn_stats = CompareSeq.rmRefAln(seqAlgn_stats)
                                
                                #CaseStats[key][locus]
                                for annKey, annItem in seqAlgn_stats['MMannotation'].items():
                                    if annKey.isdigit(): 
                                        annReads = IMGTdbIO.annotationFormat(annKey, annItem, seqAlgn_stats['alignment']) #annItem+'[D:'+seqAlgn_stats['alignment']['Donor-'+ps]
                                        if 'Intron' in annItem.split('.')[0]:
                                            CaseMatchTable[key][locus][ps]['Intron'].append(annReads)
                                        elif 'UTR' in annItem.split('.')[0]:
                                            CaseMatchTable[key][locus][ps]['UTR'].append(annReads)
                                        elif 'Exon' in annItem.split('.')[0]:
                                            if annItem.split('.')[0] in ARS_exon:
                                                CaseMatchTable[key][locus][ps]['ARS'].append(annReads)
                                            else:
                                                CaseMatchTable[key][locus][ps]['non_ARS_exon'].append(annReads) 
                           
    IMGTdbIO.save_dict2pickle(CaseMatchTable, '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/'+groupType+'_case_MatchRecord_Locus_'+locus)
        
IMGTdbIO.save_dict2pickle(CaseMatchTable, '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/'+groupType+'_case_MatchRecord')


################
# Count the cases
################
CaseMatchTable = IMGTdbIO.load_pickle2dict('../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/'+groupType+'_case_MatchRecord.pkl')

#caseID = '44107'
#CaseMatchTable[caseID]

## : paired cases HLA typing stats
fname = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/SG41_52_pairedCases_Stats.pkl'
Matching_cases_stats = IMGTdbIO.load_pickle2dict(fname)
## Matching_cases_stats['fiveLoci_paired']

TwoField_allele_stats = {}
AllField_allele_stats = {}
HLA_allele_count = {}
HLA_twoField_allele_count = {}

for caseID, Records in CaseMatchTable.items():
    if caseID in Matching_cases_stats['fiveLoci_paired']: # All_paired; fiveLoci_paired
        if IMGTdbIO.checkSubList(list(CaseMatchTable[caseID].keys()), Five_loci): # All_loci, Five_loci
            for locus in Five_loci:
                if locus in Records.keys():
                    if locus not in TwoField_allele_stats.keys():
                        TwoField_allele_stats[locus] = {}
                    if locus not in AllField_allele_stats.keys():
                        AllField_allele_stats[locus] = {}    
                    
                    if locus not in HLA_allele_count.keys():
                        HLA_allele_count[locus] = {}
                    if locus not in HLA_twoField_allele_count.keys():
                        HLA_twoField_allele_count[locus] = {}
                        
                    for ps in ['PS1', 'PS2']: 
                        ##### All field gl-string
                        typings = Records[locus][ps]['HLAtyping']
                        temp_typings = []
                        for tp in typings: 
                            if tp.find('[') != -1:
                                temp_typings.append(tp[(tp.find('[')+2):(tp.find(',')-1)])
                            else:
                                temp_typings.append(tp)
                        
                        typings = temp_typings
                        
                        '''    HLAtyping = []
                        for tp in typings:
                            if tp.find(',')!= -1:
                                tp1 = tp.replace('[\'', '')
                                tp1 = tp1.replace('\']', '')
                                tp1 = tp1.replace('\'', '')
                                ambTPlist = tp1.split(', ')
                                for tp11 in ambTPlist: 
                                    if tp11.find('/')!=-1:
                                        tp2 = tp11.split('/')
                                        HLAtyping.extend(tp2)
                                    else:
                                        HLAtyping.append(tp11)
                                #HLAtyping.extend(ambTPlist)
                            elif tp.find('/') != -1:
                                ambTPlist = tp.split('/')
                                HLAtyping.extend(ambTPlist)
                                
                            else:
                                HLAtyping.append(tp)
                        typings = HLAtyping
                        '''
                        
                        for tp in typings:
                            
                            if tp not in HLA_allele_count[locus].keys():
                                HLA_allele_count[locus][tp] = [caseID]
                            else:
                                HLA_allele_count[locus][tp].append(caseID)
                        
                        ##### Two field gl-string
                        #twoFieldtypings = Records[locus][ps]['TwoField']
                        twoFieldtypings = [':'.join((tp.split(':')[0], tp.split(':')[1])) for tp in typings]
                        twoFieldtypings = list(set(twoFieldtypings))
                            
                        for Ttp in twoFieldtypings:
                            if Ttp not in HLA_twoField_allele_count[locus].keys():
                                HLA_twoField_allele_count[locus][Ttp] = [caseID]
                            else:
                                HLA_twoField_allele_count[locus][Ttp].append(caseID)
                        
                        if not Records[locus][ps]['SameSeq']:
                            ##### All field gl-string
                            typings.sort()
                            if len(typings) > 1:
                                FullField = ''
                                
                                for tp in typings:
                                    if len(FullField) == 0:
                                        FullField += tp
                                    else:
                                        FullField += '_'+tp
                            else: 
                                FullField = typings[0]
                                
                            if FullField not in AllField_allele_stats[locus].keys():
                                AllField_allele_stats[locus][FullField] = {'ARS': [], 'ARS_nonSynonymous': [],
                                                     'non_ARS_exon': [], 'non_ARS_exon_nonSyn': [], 'Intron': []}
                            
                            ##### Two field field gl-string
                            twoFieldtypings = Records[locus][ps]['TwoField']
                            twoFieldtypings.sort()
                            if len(twoFieldtypings) > 1:
                                TwoField = ''
                                for Ttp in twoFieldtypings:
                                    if len(TwoField) ==0:
                                        TwoField += Ttp
                                    else:
                                        TwoField += '_' + Ttp
                            else:
                                TwoField = twoFieldtypings[0]
                            if TwoField not in TwoField_allele_stats[locus].keys():
                                TwoField_allele_stats[locus][TwoField] = {'ARS': [], 'ARS_nonSynonymous': [], 
                                                     'non_ARS_exon': [], 'non_ARS_exon_nonSyn': [], 'Intron': []}
        
                            if len(Records[locus][ps]['ARS']) > 0:
                                AllField_allele_stats[locus][FullField]['ARS'].append(caseID)
                                TwoField_allele_stats[locus][TwoField]['ARS'].append(caseID)
                                
                                #ExonList = [exn.split('.')[0] for exn in Records[locus][ps]['ARS']]
                                #ExonList = list(set(ExonList))
                                #Exons = ''
                                #for i in range(len(ExonList)):
                                #    if i == 0:
                                #        Exons += ExonList[i]
                                #    else:
                                #        Exons += ', ' + ExonList[i]
                                        
                                if not IMGTdbIO.checkSynonymMutation(FullField, 'Exon1, Exon2, Exon3', db_fp): # non-synonymous
                                    AllField_allele_stats[locus][FullField]['ARS_nonSynonymous'].append(caseID)
                                    TwoField_allele_stats[locus][TwoField]['ARS_nonSynonymous'].append(caseID)
                                    
                            if len(Records[locus][ps]['non_ARS_exon']) > 0:
                                AllField_allele_stats[locus][FullField]['non_ARS_exon'].append(caseID)
                                TwoField_allele_stats[locus][TwoField]['non_ARS_exon'].append(caseID)
                                
                                ExonList = [exn.split('.')[0] for exn in Records[locus][ps]['non_ARS_exon']]
                                ExonList = list(set(ExonList))
                                Exons = ''
                                for i in range(len(ExonList)):
                                    if i == 0:
                                        Exons += ExonList[i]
                                    else:
                                        Exons += ', ' + ExonList[i]
                                        
                                if not IMGTdbIO.checkSynonymMutation(FullField, Exons, db_fp): # non-synonymous
                                    AllField_allele_stats[locus][FullField]['non_ARS_exon_nonSyn'].append(caseID)
                                    TwoField_allele_stats[locus][TwoField]['non_ARS_exon_nonSyn'].append(caseID)

                            if len(Records[locus][ps]['Intron']) > 0:
                                AllField_allele_stats[locus][FullField]['Intron'].append(caseID)
                                TwoField_allele_stats[locus][TwoField]['Intron'].append(caseID)

out_fp = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/fiveLoci_paired_cases_0328/' # /fiveLoci_paired_cases/ ; allLoci_paired_cases
Loci_num = 'fiveLoci_'
for locus in Five_loci:
    with open(out_fp +'AllField/MisMatchCases/HLA_'+locus+'_'+Loci_num+'PairedCases_wNonSyn_0328.csv', 'w+') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in AllField_allele_stats[locus].items():
           writer.writerow([key, value['ARS'], value['ARS_nonSynonymous'], value['non_ARS_exon'], value['non_ARS_exon_nonSyn'], value['Intron']])
    
    with open(out_fp + 'TwoField/MisMatchCases/HLA_'+locus+'_'+Loci_num+'PairedCases_wNonSyn_0328.csv', 'w+') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in TwoField_allele_stats[locus].items():
           writer.writerow([key, value['ARS'], value['ARS_nonSynonymous'], value['non_ARS_exon'], value['non_ARS_exon_nonSyn'], value['Intron']])

######## counts
HLA_AF_allele_name = {}
HLA_AF_allele_count = {}
HLA_AF_allele_name_count = {}
HLA_TF_allele_name = {}  
HLA_TF_allele_count = {}  
HLA_TF_allele_name_count = {}       
for locus in Five_loci:
    if locus not in HLA_AF_allele_name.keys():
        HLA_AF_allele_name[locus] = ()
    if locus not in HLA_AF_allele_count.keys():
        HLA_AF_allele_count[locus] = []
    if locus not in HLA_AF_allele_name_count.keys():
        HLA_AF_allele_name_count[locus] = {}
    if locus not in HLA_TF_allele_name.keys():
        HLA_TF_allele_name[locus] = ()
    if locus not in HLA_TF_allele_count.keys():
        HLA_TF_allele_count[locus] = []
    if locus not in HLA_TF_allele_name_count.keys():
        HLA_TF_allele_name_count[locus] = {}
    
    for Allele, cases in HLA_allele_count[locus].items():
        HLA_AF_allele_name[locus] += (Allele,)
        HLA_AF_allele_count[locus].append(len(list(set(cases))))
        HLA_AF_allele_name_count[locus][Allele] = len(list(set(cases)))
    
    for TAllele, Tcases in HLA_twoField_allele_count[locus].items():
        HLA_TF_allele_name[locus] += (TAllele,)
        HLA_TF_allele_count[locus].append(len(list(set(Tcases))))
        HLA_TF_allele_name_count[locus][TAllele] = len(list(set(Tcases)))
        
    print('-'*40)
    print('Locus '+locus+ ' allele count: '+str(len(HLA_allele_count[locus])))
    print('Locus '+locus+ ' Two-Field allele count: '+str(len(HLA_twoField_allele_count[locus])))
    
### Write to csv
out_fp = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/fiveLoci_paired_cases_0328/'
for locus in Five_loci:
    with open(out_fp + 'AllField/AlleleCount/HLA_'+locus+'_'+Loci_num+'Cases_Corrected_0328.csv', 'w+') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in HLA_AF_allele_name_count[locus].items():
           writer.writerow([key, value])
    
    with open(out_fp + 'TwoField/AlleleCount/HLA_'+locus+'_'+Loci_num+'Cases_Corrected_0328.csv', 'w+') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in HLA_TF_allele_name_count[locus].items():
           writer.writerow([key, value])

### Mismatch counts

HLA_AF_allele_MisMatch_name_count = {}

HLA_TF_allele_MisMatch_name_count = {}       
for locus in Five_loci:
    
    if locus not in HLA_AF_allele_MisMatch_name_count.keys():
        HLA_AF_allele_MisMatch_name_count[locus] = {}
    
    if locus not in HLA_TF_allele_MisMatch_name_count.keys():
        HLA_TF_allele_MisMatch_name_count[locus] = {}
    
    for Allele, cases in AllField_allele_stats[locus].items():
        
        HLA_AF_allele_MisMatch_name_count[locus][Allele] = {'ARS': len(list(set(cases['ARS']))), 
                                        'ARS_nonSynonymous': len(list(set(cases['ARS_nonSynonymous']))),
                                'non_ARS_exon': len(list(set(cases['non_ARS_exon']))),
                                'non_ARS_exon_nonSyn':len(list(set(cases['non_ARS_exon_nonSyn']))),
                                'Intron': len(list(set(cases['Intron'])))}
    
    for TAllele, Tcases in TwoField_allele_stats[locus].items():
        
        HLA_TF_allele_MisMatch_name_count[locus][TAllele] = {'ARS': len(list(set(Tcases['ARS']))), 
                                         'ARS_nonSynonymous': len(list(set(Tcases['ARS_nonSynonymous']))),
                                'non_ARS_exon': len(list(set(Tcases['non_ARS_exon']))),
                                'non_ARS_exon_nonSyn':len(list(set(Tcases['non_ARS_exon_nonSyn']))),
                                'Intron': len(list(set(Tcases['Intron'])))}
     
    print('=*'*40)
    print('Locus '+locus+ ' Mismatched allele count: '+str(len(HLA_AF_allele_MisMatch_name_count[locus])))
    print('Locus '+locus+ ' Mismatched Two-Field allele count: '+str(len(HLA_TF_allele_MisMatch_name_count[locus])))

## save
out_fp = out_fp = '../Output/SG41_52/2018/IMGTv3310/SG41_52_DRpair_Stats/fiveLoci_paired_cases_0328/' # /fiveLoci_paired_cases/ ; allLoci_paired_cases
for locus in Five_loci:
    with open(out_fp + 'AllField/MisMatchCount/HLA_AllField_'+locus+'_MismatchedCounts_'+Loci_num+'Cases_Corrected_wNonSyn_0328.csv', 'w+') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in HLA_AF_allele_MisMatch_name_count[locus].items():
            writer.writerow([key, str(value['ARS']), str(value['ARS_nonSynonymous']), str(value['non_ARS_exon']), str(value['non_ARS_exon_nonSyn']), str(value['Intron'])])
    
    with open(out_fp +'TwoField/MisMatchCount/HLA_TwoField_'+locus+'_MismatchedCounts_'+Loci_num+'Cases_Corrected_wNonSyn_0328.csv', 'w+') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in HLA_TF_allele_MisMatch_name_count[locus].items():
            writer.writerow([key, str(value['ARS']), str(value['ARS_nonSynonymous']), str(value['non_ARS_exon']), str(value['non_ARS_exon_nonSyn']), str(value['Intron'])])



####################
## Check mismatched cases, if they are PS swapped case
## AllField_allele_stats[locus][Typing]['ARS']
## CaseMatchTable[caseID][Typing]
# A*68:01:01:02_A*68:01:02:02   -- 22 cases
# A*68:01:01:02_A*68:01:02:02   -- 11 cases
# A*02:01:01:01/A*02:01:01:02L_A*26:01:01 -- 3 cases
# A*24:02:01:01_A*24:03:01
# A*03:01:01:01_A*24:03:01
# A*01:01:01:01_A*03:01:01:01
# A*11:01:01:01_A*29:02:01:01
# A*68:01:02:02_A*68:02:01:01
# A*01:01:01:01_A*23:01:01
# A*02:01:01:01/A*02:01:01:02L_A*02:01:05
# A*02:01:01:01/A*02:01:01:02L_A*23:01:01
# A*02:01:01:01/A*02:01:01:02L_A*30:01:01
# A*01:01:01:01_A*02:01:01:01/A*02:01:01:02L
# A*03:01:01:01_A*31:01:02:01
# A*02:01:01:01/A*02:01:01:02L_A*31:01:02:01
# A*02:01:01:01/A*02:01:01:02L_A*02:06:01
# A*03:01:01:01_A*29:02:01:02
# A*01:01:01:01_A*30:01:01
# A*03:01:01:01_A*24:02:01:01
# A*02:01:01:01/A*02:01:01:02L_A*03:01:01:01
# A*02:01:01:01/A*02:01:01:02L_A*02:05:01
# A*02:01:01:01/A*02:01:01:02L_A*11:01:01:01

ARS_mismatched_alleles_A = ['A*68:01:01:02_A*68:01:02:02', 'A*68:01:01:02_A*68:01:02:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*26:01:01',
                            'A*24:02:01:01_A*24:03:01',
                            'A*03:01:01:01_A*24:03:01',
                            'A*01:01:01:01_A*03:01:01:01',
                            'A*11:01:01:01_A*29:02:01:01',
                            'A*68:01:02:02_A*68:02:01:01',
                            'A*01:01:01:01_A*23:01:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*02:01:05',
                            'A*02:01:01:01/A*02:01:01:02L_A*23:01:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*30:01:01',
                            'A*01:01:01:01_A*02:01:01:01/A*02:01:01:02L',
                            'A*03:01:01:01_A*31:01:02:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*31:01:02:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*02:06:01',
                            'A*03:01:01:01_A*29:02:01:02',
                            'A*01:01:01:01_A*30:01:01',
                            'A*03:01:01:01_A*24:02:01:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*03:01:01:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*02:05:01',
                            'A*02:01:01:01/A*02:01:01:02L_A*11:01:01:01',
                            'A*68:01:01:02_A*68:01:02:01']
ARS_mismatched_alleles_B = ['B*08:01:01_B*08:01:02',
                            'B*39:01:01:02L/B*39:01:01:03_B*39:06:02',
                            'B*35:01:01:01/B*35:01:01:02_B*35:03:01',
                            'B*18:01:01:02_B*58:01:01',
                            'B*38:02:01_B*38:02:02',
                            'B*15:01:01:01_B*55:01:01',
                            'B*51:01:01:01_B*51:08:01']
ARS_mismatched_alleles_C = ['C*04:01:01:01_C*06:02:01:01', 'C*03:03:01_C*03:03:02']
ARS_mismatched_alleles_DRB1 = ['DRB1*01:01:01_DRB1*01:02:01', 
                               'DRB1*03:01:01:01/DRB1*03:01:01:02_DRB1*03:01:02',
                               'DRB1*11:01:01_DRB1*11:01:02',
                               'DRB1*15:01:01:01/DRB1*15:01:01:02/DRB1*15:01:01:03/DRB1*15:01:01:04_DRB1*15:03:01:01/DRB1*15:03:01:02']
ARS_mismatched_alleles_DQB1 =['DQB1*03:01:01:01_DQB1*03:02:01',
                              'DQB1*02:02:01_DQB1*03:03:02:01',
                              'DQB1*03:01:01:02/DQB1*03:01:01:03_DQB1*05:01:01:02/DQB1*05:01:01:03',
                              'DQB1*05:01:01:02/DQB1*05:01:01:03_DQB1*06:02:01',
                              'DQB1*06:02:01_DQB1*06:03:01',
                              'DQB1*05:02:01_DQB1*06:02:01',
                              'DQB1*06:04:01_DQB1*06:09:01',
                              'DQB1*03:01:01:01/DQB1*03:01:01:02/DQB1*03:01:01:03_DQB1*03:02:01',
                              'DQB1*03:01:01:02/DQB1*03:01:01:03_DQB1*03:03:02:01',
                              'DQB1*02:02:01_DQB1*03:03:02:01/DQB1*03:03:02:02/DQB1*03:03:02:03/DQB1*03:03:02:04',
                              'DQB1*03:02:01_DQB1*03:05:01',
                              'DQB1*03:01:01:02/DQB1*03:01:01:03_DQB1*06:03:01',
                              'DQB1*05:01:01:02/DQB1*05:01:01:03_DQB1*06:09:01',
                              'DQB1*03:01:01:01_DQB1*06:02:01',
                              'DQB1*05:01:01_DQB1*05:02:01',
                              'DQB1*03:19_DQB1*06:02:01',
                            'DQB1*02:01:01_DQB1*06:02:01',
                            'DQB1*06:02:01_DQB1*06:09:01',
                            'DQB1*03:02:01_DQB1*04:02:01',
                            'DQB1*03:01:01_DQB1*03:02:01',
                            'DQB1*03:01:01:02/DQB1*03:01:01:03_DQB1*03:02:01',
                            'DQB1*03:01:01_DQB1*05:01:01:02/DQB1*05:01:01:03',
                            'DQB1*06:02:01_DQB1*06:04:01']
bothMM_cases = []
singleMM_cases = []
locus = 'DQB1'
for allele in ARS_mismatched_alleles_DQB1:
    for cases in AllField_allele_stats[locus][allele]['ARS']:
        if cases in Matching_cases_stats[locus+'_both_Seqmm']:
            bothMM_cases.append(cases)
        if cases in Matching_cases_stats[locus+'_one_Seqmm']:
            singleMM_cases.append(cases)


## plots
import datetime
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

with PdfPages('../Output/Stats/HLA_counts_FiveLoci_paired_cases.pdf') as pdf:
    plt.figure(figsize=(9, 3))
    for locus in All_loci:
        AF_sorted_allele = ()
        AF_count = []
        for allele in sorted(HLA_AF_allele_name_count[locus], key = HLA_AF_allele_name_count[locus].get, reverse = True):
            AF_sorted_allele += (allele,) 
            AF_count.append( HLA_AF_allele_name_count[locus][allele])
        
        TF_sorted_allele = ()
        TF_count = []
        for allele in sorted(HLA_TF_allele_name_count[locus], key = HLA_TF_allele_name_count[locus].get, reverse = True):
            TF_sorted_allele += (allele,) 
            TF_count.append( HLA_TF_allele_name_count[locus][allele])
        
        y_pos = np.arange(len(AF_sorted_allele))
 
        plt.bar(y_pos, AF_count, align='center', alpha=0.5)
        
        plt.xticks(y_pos, AF_sorted_allele, rotation = 90)
        plt.ylabel('Case count')
        plt.title('HLA-'+locus+" Alleles")
        
        pdf.savefig()
        plt.close()

################## 
## Matching counts based on ARS
#################
All_caseIDs = list(set(CaseMatchTable))
MatchingType_AllLoci = {}
MatchingType_fiveLoci = {}
fiveLoci_paired_cases = [] # after removing non-audit cases # 2756
AllLoci_paired_cases = []
for key in All_caseIDs:
    if key in Matching_cases_stats['fiveLoci_paired']:
        
        if IMGTdbIO.checkSubList(list(CaseMatchTable[key].keys()), Five_loci):
            fiveLoci_paired_cases.append(key)
            Match_allele_count_fiveLoci = 10
            for locus in Five_loci:
                
                for ps in ['PS1', 'PS2']:
                    if not CaseMatchTable[key][locus][ps]['SameSeq']:# check ARS
                        if len(CaseMatchTable[key][locus][ps]['ARS'])>0:
                            Match_allele_count_fiveLoci -= 1
            if str(Match_allele_count_fiveLoci) in MatchingType_fiveLoci.keys():
                MatchingType_fiveLoci[str(Match_allele_count_fiveLoci)]['count'] += 1
                MatchingType_fiveLoci[str(Match_allele_count_fiveLoci)]['cases'].append(key)
            else:
                MatchingType_fiveLoci[str(Match_allele_count_fiveLoci)] = {'count': 1, 'cases': [key]}

    if key in Matching_cases_stats['All_paired']:
        if IMGTdbIO.checkSubList(list(CaseMatchTable[key].keys()), All_loci):
            AllLoci_paired_cases.append(key)
            Match_allele_count_AllLoci = 12
            
            for locus in All_loci:
                
                for ps in ['PS1', 'PS2']:
                    if not CaseMatchTable[key][locus][ps]['SameSeq']:# check ARS
                        if len(CaseMatchTable[key][locus][ps]['ARS'])>0:
                            Match_allele_count_AllLoci -= 1
            if str(Match_allele_count_AllLoci) in MatchingType_AllLoci.keys():
                MatchingType_AllLoci[str(Match_allele_count_AllLoci)]['count'] += 1
                MatchingType_AllLoci[str(Match_allele_count_AllLoci)]['cases'].append(key)
            else:
                MatchingType_AllLoci[str(Match_allele_count_AllLoci)] = {'count': 1, 'cases': [key]}     

print('-*'*30)    
print('All 6 loci paired cases:')
for key, item in MatchingType_AllLoci.items():
    print(key+'/12 ARS matched cases: '+str(MatchingType_AllLoci[key]['count']))      
print('-*'*30)          
print('Five loci paired cases:')
for key, item in MatchingType_fiveLoci.items():
    print(key+'/10 ARS matched cases: '+str(MatchingType_fiveLoci[key]['count']))      
print('-*'*30)     
 
'''
save: 
All_caseIDs = list(set(CaseMatchTable))
MatchingType_AllLoci = {}
MatchingType_fiveLoci = {}
fiveLoci_paired_cases = [] # after removing non-audit cases # 2756
AllLoci_paired_cases = []
'''
saveOBJ = {'All_caseIDs': All_caseIDs, 'MatchingType_AllLoci': MatchingType_AllLoci, 
           'MatchingType_fiveLoci':MatchingType_fiveLoci, 'CaseMatchTable': CaseMatchTable,
           'fiveLoci_paired_cases': fiveLoci_paired_cases,
           'AllLoci_paired_cases':AllLoci_paired_cases}

IMGTdbIO.save_dict2pickle(saveOBJ, '../Output/Stats/Final_all_case_Mathcin_Info')


############
import datetime
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Create the PdfPages object to which we will save the pages:
# The with statement makes sure that the PdfPages object is closed properly at
# the end of the block, even if an Exception occurs.
with PdfPages('multipage_pdf.pdf') as pdf:
    plt.figure(figsize=(3, 3))
    plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
    plt.title('Page One')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    plt.rc('text', usetex=True)
    plt.figure(figsize=(8, 6))
    x = np.arange(0, 5, 0.1)
    plt.plot(x, np.sin(x), 'b-')
    plt.title('Page Two')
    pdf.savefig()
    plt.close()

    plt.rc('text', usetex=False)
    fig = plt.figure(figsize=(4, 5))
    plt.plot(x, x*x, 'ko')
    plt.title('Page Three')
    pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
    plt.close()

    # We can also set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d['Title'] = 'Multipage PDF Example'
    d['Author'] = u'Jouni K. Sepp\xe4nen'
    d['Subject'] = 'How to create a multipage pdf file and set its metadata'
    d['Keywords'] = 'PdfPages multipage keywords author title subject'
    d['CreationDate'] = datetime.datetime(2009, 11, 13)
    d['ModDate'] = datetime.datetime.today()