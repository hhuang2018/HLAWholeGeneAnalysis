#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Functions:
      1. Read IMGT/HLA aligned sequences file
      2. Parse Exon/Intron sequences from the IMGT database
      3. Convert the IMGT database into the Dictionary sturcutre
      4. Save the IMGT Dictionary into a Sqlite3 file (or a Pickle file)
      5. Load a IMGT Sqlite3 (or Pickle) file into a Dictionary structure
      6. Build IMGT sequence dattabase (SQLite3)
"""
from os import path
from collections import defaultdict
import re
import pickle # import csv
import utils.IMGTtools as imgt
import sqlite3
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"


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
    ref_HLAtyping = imgt.removeWhiteSpace(alignment_lines[seqLineIndex[0]])[0] # re.sub(" ","", alignment_lines[seqLineIndex[0]].rstrip().split("  ")[0])
    counter = 0
    for LineIndex in range(len(seqLineIndex)):
        
        # Reference sequence
        ref_Aligned_seqs = imgt.removeWhiteSpace(alignment_lines[seqLineIndex[LineIndex]])[1]# re.sub(" ", "", alignment_lines[seqLineIndex[LineIndex]].rstrip().split("  ")[2])
        
        type_index = 0
        while type_index < total_seq_num:
            
            #HLAtyping, Aligned_seq_temp = imgt.removeWhiteSpace(alignment_lines[seqLineIndex[LineIndex] + type_index])#[0] # re.sub(" ","", alignment_lines[seqLineIndex[LineIndex] + type_index].rstrip().split("  ")[0])
            try:
                HLAtyping, Aligned_seq_temp = imgt.removeWhiteSpace(alignment_lines[seqLineIndex[LineIndex] + type_index])# [1] # re.sub(" ", "", removeAllpattern(temp_seq, "")[])
            except IndexError:
                Aligned_seq_temp = ''  # if the sequence doesn't exist in this region
                    
            # convert aligned "-" into corresponding nucleotide
            if HLAtyping != ref_HLAtyping and len(Aligned_seq_temp)>0:
                mat_index = imgt.findCharacter(Aligned_seq_temp, '-') # [ind for ind, x in enumerate(list(Aligned_seq_temp)) if x == '-']
                Aligned_seqs = list(Aligned_seq_temp)
                for x in mat_index: 
                    Aligned_seqs[x] = list(ref_Aligned_seqs)[x]
                Aligned_seqs = "".join(Aligned_seqs)
            else:
                Aligned_seqs = Aligned_seq_temp
                
            # convert "." into gap symbol "-"
            gap_index = imgt.findCharacter(Aligned_seqs, '.') # [ind for ind, x in enumerate(list(Aligned_seq_temp)) if x == '-']
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

def parseExonSequences(seq_db, dbType = "CDS"):
    """
    Parse exon sequences from CDS or genomic sequences 
    seq_db: Dictionary structure -- aligned
    dbType: "CDS" or "genomic" 
    """
    Exon_db = {}
    
    if dbType == "CDS":
        Exon_db = defaultdict(dict)
        for Typings, Seqs in seq_db.iteritems():
            boundaryIndex = imgt.findCharacter(Seqs["cDNA"], "|")
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

def parseExonSequences_py3(seq_db, dbType = "CDS"):
    """
    Parse exon sequences from CDS or genomic sequences 
    seq_db: Dictionary structure -- aligned
    dbType: "CDS" or "genomic" 
    """
    Exon_db = {}
    
    if dbType == "CDS":
        Exon_db = defaultdict(dict)
        for Typings, Seqs in seq_db.items():
            
            temp_Seqs = re.sub("-", "", Seqs["cDNA"])
            temp_segments = temp_Seqs.split("|")
            
            numExons = len(temp_segments)
            temp_db = []
            for index in range(numExons):
                temp_db.append({Typings:{'Exon'+str(index+1):"".join(temp_segments[index])}})
           
            for item in temp_db:
                for k, v in item.items():
                    Exon_db[k].update(v)
   
    return(Exon_db)

def readIMGT_alleleList(version, input_fp):
    '''
    Read Allele list and return as a dictionary
    HLA gl-string: HLA ID 
    '''
    alleleList_lines = open(input_fp+"Allelelist."+version+".txt").readlines() ### open(filename, 'r')
    
    alleleList = {}
    for LineIndex in range(len(alleleList_lines)):
        alleleList[alleleList_lines[LineIndex].rstrip().split(" ")[1]] = alleleList_lines[LineIndex].rstrip().split(" ")[0]
        
    return(alleleList)
        
def extractValues(fieldName, keyName):
    """
    Extract the value from a dictionary given a key name
    """
    try: 
        filedValue = fieldName[keyName]
    except KeyError:
        filedValue = ""
    return(filedValue)
    
def IMGTdb_2_dict(HLA_locus = "A", version = "3250", input_fp = "../IMGTHLA/"):
    """
    Convert IMGT aligned genomic and CDS sequence files into dictionary structure.
    """
    
    alleleList = readIMGT_alleleList(version, input_fp)
    if version == "3250":
        if HLA_locus in ["A", "B", "C"]: # class I
            locus = HLA_locus
        elif HLA_locus in ["DRB1", "DQB1", "DPB1"]: # class II
            locus = re.sub("1", "", HLA_locus)
    else:
        locus = HLA_locus
    # genomic alignment file    
    filename = input_fp + "/alignments/" + locus + "_gen.txt" 
    if path.exists(filename):
        gDNA_alignment = read_IMGT_alignment(filename, 'gDNA')
  
    # CDS alignemtn file
    filename = input_fp + "/alignments/" + locus + "_nuc.txt"
    if path.exists(filename):    
        CDS_alignment = read_IMGT_alignment(filename, 'cDNA', 3, 4)
    else:
        print("File %s does not exist.", filename)
        CDS_alignment = {}
    
    if len(CDS_alignment)>0:
        try:  ## Python 2.7
            ExonSequences = parseExonSequences(CDS_alignment)
        except AttributeError:## Python 3.6
            ExonSequences = parseExonSequences_py3(CDS_alignment)
            
    # protein alignment file
    filename = input_fp + "/alignments/" + locus + "_prot.txt"
    if path.exists(filename): 
        protein_alignment = read_IMGT_alignment(filename, 'Prot', 2, 3)
#    else:
#        protein_alignment = {}
    
 # merge into one single Dictionary structure
    combined_alignments = [gDNA_alignment, CDS_alignment, ExonSequences, protein_alignment]

    combined_dict = defaultdict(dict)
    try: # python 2.7
        for item in combined_alignments:
            for k, v in item.iteritems():
                combined_dict[k].update(v)
    except AttributeError: ## Python 3.6
        for item in combined_alignments:
            for k, v in item.items():
                combined_dict[k].update(v)
        
    return(combined_dict, alleleList)

## buil IMGT sqlite3 database 
def buildIMGTsql(Locus, version = "3250", output_fp = "../Database/"):
    """
    Build a Sqllite3 database for each locus: A, B, C, DRB1, DQB1, DPB1
    Include - HLA-ID
              HLA gl-string, 
              genome sequence, CDS sequence, 
              aligned genome sequence, aligned CDS sequence, 
              protein sequence, algined protein sequence,
              Exons and introns sequences (optional)
    """
    # output_fp = "Database/"
    file_name = "IMGT-" + version + "_HLA-" + Locus + ".db"

    seq_db, alleleList = IMGTdb_2_dict(Locus, version)
    
#   alleleList = readIMGT_alleleList(version, input_fp)
    if path.exists(output_fp+file_name):
       # if the the database already exists, then check if it needs to be updated
       conn = sqlite3.connect(output_fp+file_name)
       c = conn.cursor()
    else: # build a new database
       conn = sqlite3.connect(output_fp+file_name)
       c = conn.cursor()
       c.execute('''CREATE TABLE Sequences
                 (AlleleID text, HLATyping text, AlignedCDS text, AlignedGenomSeq text, UnalignedGenomSeq text, SeqAnnotation text, Exon1 text, Exon2 text, Exon3 text,
                 Exon4 text, Exon5 text, Exon6 text, Exon7 text, Exon8 text, Protein text)''')

 
    # Create table
    #c.execute('''CREATE TABLE Sequences
    #          (HLATyping text, AlignedCDS text, AlignedGenomSeq text, Exon1 text, Exon2 text, Exon3 text,
    #           Exon4 text, Exon5 text, Exon6 text, Exon7 text, Exon8 text, Protein text)''')

    # Insert a row of data
    # If the HLA typing record exists, then update or add its corresponding sequence
    # otherwise, insert a new record
    # 
    # Do this instead
    # t = (symbol,)
    # c.execute('select * from stocks where symbol=?', t)
    #
    # Larger example
    # for t in [('2006-03-28', 'BUY', 'IBM', 1000, 45.00),
    #          ('2006-04-05', 'BUY', 'MSOFT', 1000, 72.00),
    #          ('2006-04-06', 'SELL', 'IBM', 500, 53.00),
    #          ]:
    #    c.execute('insert into stocks values (?,?,?,?,?)', t)
    
    counter = 0
    for Typing, Value in seq_db.items():
        
        try:
            if Typing.index(Locus) >= 0:
                proceed_flag = True
        except ValueError:
            proceed_flag = False
            
        if proceed_flag and Typing[0] == Locus[0]:
            ID = alleleList[Typing]
            AlignedCDS = extractValues(Value, "cDNA")
            AlignedGenomSeq = extractValues(Value, "gDNA")
            
            Exon1 = extractValues(Value, "Exon1")
            Exon2 = extractValues(Value, "Exon2")
            Exon3 = extractValues(Value, "Exon3")
            Exon4 = extractValues(Value, "Exon4")
            Exon5 = extractValues(Value, "Exon5")
            Exon6 = extractValues(Value, "Exon6")
            Exon7 = extractValues(Value, "Exon7")
            Exon8 = extractValues(Value, "Exon8")
             
            if AlignedGenomSeq !='':
                UnalignedGenomSeq = re.sub("-", "", AlignedGenomSeq)
                
                annotatedSeq = UnalignedGenomSeq.split("|")
                temp_SeqAnnotation = []
                for ind in range(len(annotatedSeq)):
                    if ind%2 == 0: ## intron - negative numbers
                        temp_SeqAnnotation.append((str(int(-ind/2))+' ')*len(annotatedSeq[ind]))
                        #print(len(annotatedSeq[ind]))
                    else: # Exon - positive numbers
                        temp_SeqAnnotation.append((str(int((ind+1)/2))+' ')*len(annotatedSeq[ind]))
                        #print(len(annotatedSeq[ind]))
                        if extractValues(Value, "Exon"+str(int((ind+1)/2))) != annotatedSeq[ind]: # double check the exons sequences
                            print(Typing + " Exon "+ str(int((ind+1)/2))+ " sequence doesn't match to the Genomic sequence. \n Please Double check!!")
                SeqAnnotation = "".join(temp_SeqAnnotation)
            else:
                UnalignedGenomSeq = ''
                SeqAnnotation = ''
            Protein = extractValues(Value, "Prot")
            #ExonIndex = extractValues(Value, "ExonIndex")
            record = (ID, Typing, AlignedCDS, AlignedGenomSeq, UnalignedGenomSeq, SeqAnnotation, Exon1, Exon2, Exon3, Exon4, Exon5, Exon6, Exon7, Exon8, Protein)
            
            c.execute('INSERT INTO Sequences VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', record)
            counter += 1
            # Save (commit) the changes
    conn.commit()
    
    print("\nTotal allele count at locus HLA-"+ Locus + " is " + str(counter))
            
    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()
    
def readIMGTsql(HLAtyping, db_fp = "Database/", field = '*', version = "3310",unaligned = False):
    """
    Load a Sqllite3 database [Loci: A, B, C, DRB1, DQB1, DPB1]
    Include - HLA gl-string, 
              genome sequence, CDS sequence, 
              aligned genome sequence, aligned CDS sequence, 
              protein sequence, algined protein sequence,
              Exons and introns sequences (optional)
    """ 
    #HLAtyping = "A*01:01:01:02N"
    filename = db_fp + "IMGT-" + version + "_HLA-" + HLAtyping.split("*")[0]+".db"
    #field = 'Exon2, Exon3'
    if path.exists(filename):
        con = sqlite3.connect(filename)
        cur = con.cursor()
        t = (HLAtyping,)
        cur.execute('SELECT ' + field + ' FROM Sequences WHERE HLATyping = ?', t)
        sequences_temp = cur.fetchone()
        
        if sequences_temp == None:
            
            t = (HLAtyping+":01%",)
            cur.execute('SELECT ' + field + ' FROM Sequences WHERE HLATyping LIKE ?', t)
            sequences_temp = cur.fetchone()
            if sequences_temp == None:
                t = (HLAtyping+":02%",)
                cur.execute('SELECT ' + field + ' FROM Sequences WHERE HLATyping LIKE ?', t)
                sequences_temp = cur.fetchone()
                if sequences_temp == None:
                    print("No record of " + HLAtyping + "\nPlease check the HLA typing.")
                    sequences = ['', '', '', '', '', '', '', '']
                else:
                    if unaligned:
                        sequences = [re.sub("-", "", seq) for seq in sequences_temp]
                    else:
                        sequences = sequences_temp
            else:
                if unaligned:
                    sequences = [re.sub("-", "", seq) for seq in sequences_temp]
                else:
                    sequences = sequences_temp
           
        else:
            if unaligned:
                sequences = [re.sub("-", "", seq) for seq in sequences_temp]
            else:
                sequences = sequences_temp
            
        cur.close()
        return(sequences)
    else:
        print("No database is available. Please build an SQL database first. \n" +
              "To build a new SQL database, use the following command:\n>>> buildIMGTsql(\""+ HLAtyping.split("*")[0] +"\")")
 
    
################## pickle output 
def save_dict2pickle(dict_obj, fname):
    """
    Save the IMGT database into a pickle file
    """
    with open(fname + '.pkl', 'wb') as fileHandle:
        pickle.dump(dict_obj, fileHandle, pickle.HIGHEST_PROTOCOL)

def load_pickle2dict(fname):
    """
    Load the the IMGT database pickle file into a Dictionary structure
    """
    with open(fname, 'rb') as fileHandle:
        return pickle.load(fileHandle)

###################
def Full2TwoField(TypingList):
    '''
    Convert multi-field gl-string into 2-field typing
    '''
    
    HLAtyping = []
    for tp in TypingList:
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
                                
    twoFieldList = [':'.join((tp.split(':')[0], tp.split(':')[1])) for tp in HLAtyping]
    twoFieldList = list(set(twoFieldList)) # remove duplicates
    
    return twoFieldList

def annotationFormat(MMPos, Annotation, Alignment):
    '''
    Annotation: {'pos': 'annotation'}
    Alignment:{'AlignmSymbol': , 'RefSeq':, 'Donor-PS': , 'Recipient-PS': }
    Annotation + Donor Read + Recipient Read
    '''
    for key, item in Alignment.items():
        if 'Donor' in key:
            Donor = item[int(MMPos)]
        if 'Recipient' in key:
            Recipient = item[int(MMPos)]
    
    reformatted_annotation = Annotation + ':D'+Donor+';R'+Recipient
    return(reformatted_annotation)
    
##################
def checkSubList(List1, List2):
    '''
    Used to check if a case has all loci that listed in List2
    '''
    Result = True
    for item in List2:
        if item in List1:
            Result = Result and True
        else: 
            Result = Result and False
    return(Result)
    
##################
def checkSynonymMutation(Allele, Exons, db_fp = 'Database/'):
    '''
    Check if the nucleotide change is synonymous or not; IGMT/HLA alleles
    '''
    typingList = Allele.split('_')
    HLAtyping = []
    for tp in typingList:
        if tp.find('/') != -1:
            ambTPlist = tp.split('/')
            HLAtyping.extend(ambTPlist)
        else:
            possTPlist = re.sub('[\[\'\]]', '',tp) # remove possible characters
            possTPlist = possTPlist.split(",")
            for item in possTPlist:
                #HLAtyping.extend(possTPlist)
                HLAtyping.append(item.replace(" ", ""))
    
    ExonSeq = {}
    TransSeq = {}
    
    if 'Exon1' not in Exons:
        ExonIDs = [int(s) for s in re.findall(r'\d+', Exons)]
        ExonIDs.sort(reverse = True)
        maxID = ExonIDs[0]
        Exons = ','.join(['Exon'+str(s+1) for s in range(maxID)])
    
    for tp in HLAtyping:
        Refseq = readIMGTsql(tp, db_fp, field = Exons)
        DNAseq = ''
        for seq in Refseq:
            DNAseq += seq
        ExonSeq[tp] = Seq(DNAseq, generic_dna) 
        try: 
            TransSeq[tp] = ExonSeq[tp].translate()
            Synonymous = True # Default: same protein sequence
        except CodonTable.TranslationError:
            TransSeq[tp] = '*********'
            Synonymous = False
            return(Synonymous)

    for i in range(len(HLAtyping)):
        for j in range(len(HLAtyping)):
            if i != j:
                 if TransSeq[HLAtyping[i]] != TransSeq[HLAtyping[j]]:
                     Synonymous = Synonymous and False
    return (Synonymous)