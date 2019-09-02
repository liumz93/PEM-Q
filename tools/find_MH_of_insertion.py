#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.8.21

"""find_MH_of_insertion

Usage:
    find_MH_of_insertion <vector_fasta> <vector_sam> <transloc_tab>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script help you find microhomolog sequence of insetions from vector.

Input file: vector/vector.fa
            vector alignment sam file
            transloc all tab file
Output file: insertion_mh file

Author: Mengzhu LIU
Last Update:2019.8.21

"""
import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO

def find_mh_of_insertion(vector_fasta, vector_sam, transloc_tab):
    
    transloc_file = pd.read_csv(transloc_tab,sep = '\t')
    vector_sam_file = pysam.AlignmentFile(vector_sam, 'rb')
    vector_seq = SeqIO.read(vector_fasta, "fasta").seq
    
    insertion_mh_tab = open("Inserion_mh.tab", "w")
    insertion_mh_tab.write('Qname'+"\t"+\
                        'Sequence'+"\t"+\
                        'old_inserion'+"\t"+\
                        'real_insertion'+"\t"+\
                        'forward_mh'+"\t"+\
                        'afterward_mh'+"\n")
                        
    for read in vector_sam_file:

        if read.reference_name == "vector":
            qname = read.query_name
            qname_index = transloc_file['Qname'][transloc_file['Qname'] == qname].index[0]#get qname index in transloc_file
            Sequence = transloc_file['Sequence'][qname_index]#get sequence
            
            vector_start = read.reference_start + 1 #pysam reference start from 0
            vector_end = read.reference_end

            #~~~~ find forward microhomolog ~~~~#
            num = 0
            insertion = vector_seq[vector_start-1-num:vector_end]
            old_inserion =  insertion       
            while str(insertion) in str(Sequence):
                if vector_start-1 <= 0:
                    break
                num = num + 1
                insertion = vector_seq[vector_start-1-num:vector_end]
            #get number of microhomolog
            num = num -1
            forward_mh = vector_seq[vector_start-1-num:vector_start-1]
            
            #~~~~ find afterward microhomolog ~~~~#
            num = 0
            insertion = vector_seq[vector_start-1:vector_end+num]            
            while str(insertion) in str(Sequence):
                if vector_end+num >= len(vector_seq):
                    break
                num = num + 1
                insertion = vector_seq[vector_start-1:vector_end+num]
            #get number of microhomolog
            num = num -1
            afterward_mh = vector_seq[vector_end:vector_end+num]
            
            real_insertion = forward_mh + old_inserion + afterward_mh
            
            insertion_mh_tab.write(qname+"\t"+\
                                str(old_inserion)+"\t"+\
                                str(real_insertion)+"\t"+\
                                str(forward_mh)+"\t"+\
                                str(afterward_mh)+"\t"+\
                                str(Sequence)+"\n")
                       
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='find_MH_of_insertion 1.0')
    
    kwargs = {'vector_fasta':args['<vector_fasta>'],'vector_sam':args['<vector_sam>'],'transloc_tab':args['<transloc_tab>']}
    print('[superCasQ] vector_fasta: ' + str(kwargs['vector_fasta']))
    print('[superCasQ] vector_sam: ' + str(kwargs['vector_sam']))
    print('[superCasQ] transloc_tab: ' + str(kwargs['transloc_tab']))

    find_mh_of_insertion(**kwargs)

    print("\nfind_MH_of_insertion.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
            
            
        
            