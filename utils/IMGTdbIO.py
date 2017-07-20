#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Functions:
      1. Read IMGT/HLA aligned sequences
      2.
"""
#from os import walk, environ
from Bio import SeqIO, AlignIO
import re

__author__ = "Hu Huang"
__copyright__ = "Copyright 2017, Hu Huang"
__credits__ = ["Add names"]
__license__ = "GPL"
__version__ = "dev"
__maintainer__ = "Hu Huang"
__email__ = "hwangtiger@gmail.com"

def read_IMGT_alignment(filename):
    """
    Read IMGT alignment file and convert into specific format
    """
    alignment_file = open(filename, 'r')
    alignment_file.readline()
    
    alignment_file.close()
     