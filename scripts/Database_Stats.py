#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 12:28:22 2018

@author: hhuang2
"""
import sqlite3 as sql
from utils import IMGTdbIO

version = '3310'
IMGTdbIO.buildIMGTsql('DQB1', version = version, output_fp = "../Database/")


locus = 'DQB1'
db_fp = '../Database/'

filename = db_fp + "IMGT-" + version + "_HLA-" + locus+".db"

con = sql.connect(filename)
cur = con.cursor()
field1 = 'HLATyping'
cur.execute('SELECT ' + field1 + ' FROM Sequences')
Typings_temp = cur.fetchall()
count = 0
field2 = 'AlignedGenomSeq'
for tp in Typings_temp:
    cur.execute('SELECT ' + field2 + ' FROM Sequences WHERE HLATyping = ?', tp)
    sequences_temp = cur.fetchone()
    if sequences_temp[0] != '':
        count += 1

    
con.close()