#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
"""repeats_dedup

I wish i could solve this problem as quickly as possible!

Usage:
    repeats_dedup [options] <file> <mapqual_cutoff>

Options:
-f feat, --feature_list=feat        example: -f 'Rname,Junction,Strand'
-h --help                           Show this screen.
-v --version                        Show version.

"""
import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np

def repeats_dedup(file,mapqual_cutoff,feature_list):
    
    file_tab = pd.read_csv("raw/"+file,sep = '\t', header=0)
    
    if file_tab['Qname'].count() > 0 :
        # split barcode into single base
        length = file_tab['Barcode'].str.len()
        bl = max(length)
        barcode_list = list(map(str, range(1,(bl+1))))
        x = range(0,bl)
        y = range(0,bl)
        couple = zip(x,y)
        for i,j in couple:
            file_tab.loc[:,barcode_list[i]] = file_tab['Barcode'].str[j]
        # structure: 5,9,13
        dedup_list_high_mapqual = feature_list + ['1','2','3','4','6','7','8','10','11','12','14','15']
        # low_feature_list = ['B_Qstart','B_Qend','Qstart','Qend','Qlen']
        # low_feature_list = ['Rname','Strand','Bait_end','Junction']
        dedup_list_low_mapqual = ['1','2','3','4','6','7','8','10','11','12','14','15']
        print(dedup_list_high_mapqual)
        file_tab['Length'] = length
    
        print(file_tab['Qname'].count())
        file_tab_high_mapqual = file_tab[file_tab['Prey_MQ'] >= int(mapqual_cutoff)]
        file_tab_high_mapqual.sort_values(by=['Length'], ascending=False)
        file_tab_high_mapqual = file_tab_high_mapqual.drop_duplicates(dedup_list_high_mapqual,keep='first')
        file_tab_low_mapqual = file_tab[file_tab['Prey_MQ'] < int(mapqual_cutoff)]
        file_tab_low_mapqual.sort_values(by=['Length'], ascending=False)
        file_tab_low_mapqual = file_tab_low_mapqual.drop_duplicates(dedup_list_low_mapqual,keep='first')
        file_tab = file_tab_high_mapqual.append(file_tab_low_mapqual)
        file_tab.to_csv("unique/"+file, header = True, sep = '\t', index=False)
        print(file_tab['Qname'].count())
    else:
        file_tab.to_csv("unique/"+file, header = True, sep = '\t', index=False)
        
def main():
    start_time = time()
    
    args = docopt(__doc__,version='repeats_dedup 1.0')
    feature_list = args['--feature_list'].rsplit(sep=',')
    kwargs = {'file':args['<file>'],'mapqual_cutoff':args['<mapqual_cutoff>'],'feature_list':feature_list}
    print('file: ' + str(kwargs['file']))
    print('mapqual_cutoff: ' + str(kwargs['mapqual_cutoff']))
    print('feature_list: ' + str(kwargs['feature_list']))
    
    repeats_dedup(**kwargs)
    
    print("\nrepeats_dedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()