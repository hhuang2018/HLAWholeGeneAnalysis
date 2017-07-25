#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Functions:
      1. Read IMGT/HLA aligned sequences

"""
from os import path
from collections import defaultdict
import re
import csv

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