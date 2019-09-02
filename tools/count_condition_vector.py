#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.6.17
"""count_condition

Usage:
    count_coundition <start> <end> <binsize> <basename> <tab_file>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .

Author: Mengzhu LIU
Last Update:2019.6.17

"""

import os
import pysam
import re
from time import time
from docopt import docopt
import pandas as pd
from Bio import SearchIO
from Bio.Seq import Seq

def count(start=None,end=None,binsize=None,basename=None,tab_file=None):
    start = int(start)
    end = int(end)
    binsize = int(binsize)
    data = pd.read_table(tab_file,sep = '\t',index_col=False,low_memory=False)
    count_tab = open(basename+"_count.tab", "w")
    count_tab.write("Start"+"\t"+"End"+"\t"+"Pos_total"+"\t"+"Neg_total"+"\n")
    for i in range(start,end,binsize):
        binstart = i
        binend = i+binsize
        # condition0 = data['Rname'] == chrom
        condition01 = data['Strand'] == '+'
        condition02 = data['Strand'] == '-'
        condition1 = data['Junction'] >= binstart
        condition2 = data['Junction'] < binend
        # condition3 = data['MH_end'].str.len() > 0
        #Pos_total
        datb = data[condition01][condition1][condition2]
        count3 = datb['Qname'].count()
        #Neg_total
        datb = data[condition02][condition1][condition2]
        count4 = datb['Qname'].count()
        count_tab.write(str(binstart)+"\t"+str(binend)+"\t"+str(count3)+"\t"+str(count4)+"\n")       
    count_tab.close()
    
def percentage(basename=None,tab_file=None):
    data = pd.read_table(basename+"_count.tab",sep = '\t',index_col=False,low_memory=False)
    total_pos = sum(data['Pos_total'])
    total_neg = sum(data['Neg_total'])
    total = total_pos+total_neg
    print(total,total_pos,total_neg )
    data['Pos%'] = data['Pos_total']/total*100
    data['Neg%'] = -data['Neg_total']/total*100
    data.to_csv(basename+"_count.tab",header=True,sep='\t',index=False) 
       
def main():
    args = docopt(__doc__,version='count_condition 1.0')
    
    kwargs = {'start':args['<start>'], 'end':args['<end>'], 'binsize':args['<binsize>'], 'basename':args['<basename>'], 'tab_file':args['<tab_file>']}
    print('[superCasQ] start: ' + str(kwargs['start']))
    print('[superCasQ] end: ' + str(kwargs['end']))
    print('[superCasQ] binsize: ' + str(kwargs['binsize']))
    print('[superCasQ] basename: ' + str(kwargs['basename']))
    print('[superCasQ] tab_file: ' + str(kwargs['tab_file']))

    count(**kwargs)
    percentage(str(kwargs['basename']))
    
if __name__ == '__main__':
    main()
    