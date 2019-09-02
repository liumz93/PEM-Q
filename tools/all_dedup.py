#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.6.23
"""all_dedup

Usage:
    all_dedup   <basename>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script find bait break site and translocation site from bwa align file.
And also filter reads for indel define.

Input file: transloc and indel tab file / Output file: unique transloc and indel tab file

Author: Mengzhu LIU
Last Update:2019.6.24

"""
import os
from time import time
from docopt import docopt
import pandas as pd

def dedup_all_transloc(basename):
    transloc = pd.read_csv(basename + "_transloc_all.tab",sep = '\t')
    data = transloc.drop_duplicates(['Bait_start', 'Bait_end', 'Prey_rname', 'Prey_strand', 'Junction'],keep='first')
    data.to_csv(basename + "_transloc_all_unique.tab", header = True, sep = '\t', index=False)

def main():
    start_time = time()
    
    args = docopt(__doc__,version='all_dedup 1.0')
    
    kwargs = {'basename':args['<basename>']}
    print('[superCasQ] basename: ' + str(kwargs['basename']))

    dedup_all_transloc(str(kwargs['basename']))
    print("\nall_dedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
    