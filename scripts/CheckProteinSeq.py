#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 13:24:16 2017

@author: hhuang2
"""
from utils import IMGTdbIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

typing1 = 'A*02:01:01:01'
Refseq1 = IMGTdbIO.readIMGTsql(typing1, field = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7, Exon8')

typing2 = 'A*02:05:01'
Refseq2 = IMGTdbIO.readIMGTsql(typing2, field = 'Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7, Exon8')

coding_dna = Seq(Refseq1[1], generic_dna)
coding_dna.translate()

################## 
## Check Synonymous
#################
out_fp = '../Output/Stats/allLoci_paired_cases/'
Loci_num = 'AllLoci'
for locus in All_loci:
    for Allele, Item in AllField_allele_stats[locus].items():
        if Item['ARS']
    AllField_allele_stats[locus]
    TwoField_allele_stats[locus]        
            
            
            
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
 