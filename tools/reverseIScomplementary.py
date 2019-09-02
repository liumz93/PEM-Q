#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.8.9

"""reverseIScomplementary 

Usage:
    reverseIScomplementary  <tab_file>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .

tab_file: generated from superQ2

Author: Mengzhu LIU
Last Update:2019.7.5

"""

from Bio.Seq import Seq
import pandas as pd
from time import time
from docopt import docopt

def judge_if_complementary(seq):
    
    sequence = Seq(seq)
    reverse_complement_sequence = sequence.reverse_complement()
    # reverse_sequence = seq[::-1]
    if reverse_complement_sequence == sequence:
        answer = 1
    else:
        answer = 0
    
    return(answer)
    

def count_condition(tab_file=None):
    
    # tab_file = "dicentricORhairpin.tab"
    data = pd.read_csv(tab_file,sep = '\t',index_col=False,low_memory=False)
    count = data["Microhomolog"].count()
    print("Total:",count)
    flag = 0 
    for i in range(0, count):
        homo_string = data["Microhomolog"][i]
        flag = flag + judge_if_complementary(homo_string)
    print("ReverseIScomplement:" ,flag)
    
    return()
    
def main():
    
    start_time = time()
    args = docopt(__doc__,version='reverseIScomplementary 1.0')
    
    kwargs = {'tab_file':args['<tab_file>']}
    print('[superCasQ] tab_file: ' + str(kwargs['tab_file']))

    count_condition(**kwargs)

    print("\nreverseIScomplementary.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
        