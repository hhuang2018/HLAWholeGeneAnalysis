#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Functions:
      1. Read IMGT/HLA aligned sequences

"""
import re

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

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
