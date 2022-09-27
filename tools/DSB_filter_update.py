#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.9.26
"""all_dedup

    This is a pipeline to analyze PEM-seq data or data similar, help you analyze repair outcome of your DNA library.

    Copyright (C) 2019  Mengzhu Liu

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


Author: Mengzhu LIU
Contact: liu.mengzhu128@gmail.com/liumz@pku.edu.cn

Usage:
    DSB_filter   <basename> <chromosome> <cutsite> <strand> <de_method> <yourfile>

Options:
-h --help               Show this screen.
-v --version            Show version.
<de_method>             1:sequence; 2:bait start/end & prey start/end; 3:junction

This script filter reads not come from DSB rule.

Input file: transloc tab file / Output file: transloc clean tab file

Author: Mengzhu LIU
Last Update:2019.6.24

"""
import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
from difflib import get_close_matches

def DSB_filter(basename,chromosome,cutsite,strand,de_method,yourfile):
    ### column to dedup : Bait_rname,Bait_strand,Bait_start,Bait_end,Prey_rname,Prey_strand,Prey_start,Prey_end,Rname,Strand,Junction,Sequence
    #total dataset
    transloc = pd.read_csv("results/" + basename + yourfile + ".tab",sep = '\t') 
    Full = pd.read_csv("results/" + basename + yourfile + ".tab",sep = '\t')
    transloc['Seq_len']=transloc['Sequence'].map(lambda x: 0 if x is np.nan else len(x))
    de_list = ['Rname','Strand','Junction','Prey_start','Prey_end']
    dup = transloc[transloc.duplicated(de_list,keep = False)]
    dup.to_csv(basename + "_dsb_dup.tab", header = True, sep = '\t', index=False)
    # cutsite_tab = dup[(c1) & (c2) & (c3)]
    print(dup['Qname'].count())
    #reads to be judge 
    if de_method == '1':
        de_list = ['Sequence']
    if de_method == '2':
        de_list = ['Rname','Strand','Prey_start','Prey_end','Bait_end','Barcode']
    # if de_method == '22':
    #     de_list = ['Rname','Strand','Junction','Bait_end','Prey_end']
    if de_method == '3' or de_method == '3f' or de_method == '3fe':
        de_list = ['Rname','Strand','Junction','Prey_start','Prey_end']
    
    #mark fist read of duplicated as true
    reads = dup.loc[~dup.duplicated(de_list,keep = 'first'),:]
    
    print(reads['Qname'].count())
    #drop reads around cutsite because they are true dsb reads and make this program fast
    c1 = reads['Rname'] == chromosome
    c2 = reads['Junction'] >= int(cutsite) - 20
    c3 = reads['Junction'] <= int(cutsite) + 20
    reads.drop(reads[(c1) & (c2) & (c3) ].index, inplace=True)
    reads.reset_index(level=reads.index.names, inplace=True)
    print(reads['Qname'].count())
    # #dataset without duplicated reads
    # dataset = transloc.drop_duplicates(['Rname','Strand','Junction'],keep= False)
    # dataset.reset_index(level=dataset.index.names, inplace=True)
    # print(dataset['Qname'].count())
    reads.to_csv(basename + "_dsb_reads.tab", header = True, sep = '\t', index=False)
    #set finalset
    finalset = transloc
    # finalset = dataset
    add = pd.DataFrame()
    
    for i in range(0,reads['Qname'].count()):
        
        Qname = reads['Qname'][i]
        Bait_rname = reads['Bait_rname'][i]
        Bait_strand = reads['Bait_strand'][i]
        Bait_start = reads['Bait_start'][i]
        Bait_end = reads['Bait_end'][i]
        Prey_rname = reads['Prey_rname'][i]
        Prey_strand = reads['Prey_strand'][i]
        Prey_start = reads['Prey_start'][i]
        Prey_end = reads['Prey_end'][i]
        Rname = reads['Rname'][i]
        Strand = reads['Strand'][i]
        Junction = reads['Junction'][i]
        Sequence = reads['Sequence'][i]
        Insertion = reads['Insertion'][i]
        Seq_len = reads['Seq_len'][i]
        Barcode = reads['Barcode'][i]
        
        start = reads['Junction'][i] - 5
        end = reads['Junction'][i] + 5
            
        if de_method == '1':
            de_condition1 = transloc['Sequence']==Sequence
        if de_method == '2':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Prey_start']==Prey_start) & \
                        (transloc['Prey_end']==Prey_end) & \
                        (transloc['Bait_end']==Bait_end) & \
                        (transloc['Barcode']==Barcode)
        if de_method == '3':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Prey_start']==Prey_start) & \
                        (transloc['Prey_end']==Prey_end) & \
                        (transloc['Bait_end']==Bait_end)
                        # (transloc['Bait_start']==Bait_start) & \
        if de_method == '4':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Junction']==Junction)
        if de_method == '3f':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Junction']==Junction) & \
                        (transloc['Bait_start']>=Bait_start-5) & \
                        (transloc['Bait_start']<=Bait_start+5) & \
                        (transloc['Bait_end']>=Bait_end-10) & \
                        (transloc['Bait_end']<=Bait_end+10)
        if de_method == '3fe':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Junction']==Junction) & \
                        (transloc['Bait_end']>=Bait_end-2) & \
                        (transloc['Bait_end']<=Bait_end+2)
        if de_method == '5':
            de_condition1 = (transloc['Barcode']==Barcode)

        #dataset to check
        dataset = transloc.drop(transloc[de_condition1].index)

        #condition
        if Strand == '+':
            s = '-'
        if Strand == '-':
            s = '+'
        count = dataset[(dataset['Rname'] == Rname) & (dataset['Junction'] >= start) & (dataset['Junction'] <= end)]# & (dataset['Strand'] == s)
        condition = count['Qname'].count() <= 20
        
        # condition = True
        # condition = True
        # if not condition:
        #     if Rname == "chr1":
        #         print(count['Qname'],start,end)
        #         # print(dataset[(dataset['Rname'] == Rname) & (dataset['Junction'] >= start) & (dataset['Junction'] <= end)])
        #         ccc = transloc[de_condition1]
        #         ccc.to_csv("check.tab", header = True, sep = '\t', index=False)
        
        if condition:
            if de_method == '1':
                de_condition2 = finalset['Sequence']==Sequence
            if de_method == '2':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Prey_start']==Prey_start) & \
                            (finalset['Prey_end']==Prey_end) & \
                            (finalset['Bait_end']==Bait_end) & \
                            (finalset['Barcode']==Barcode)
            if de_method == '3':
                
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Prey_start']==Prey_start) & \
                            (finalset['Prey_end']==Prey_end) & \
                            (finalset['Bait_end']==Bait_end)
                            # (finalset['Bait_start']==Bait_start) & \
                
            if de_method == '4':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Junction']==Junction)
            if de_method == '3f':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Junction']==Junction) & \
                            (finalset['Bait_start']>=Bait_start-5) & \
                            (finalset['Bait_start']<=Bait_start+5) & \
                            (finalset['Bait_end']>=Bait_end-10) & \
                            (finalset['Bait_end']<=Bait_end+10)
            if de_method == '3fe':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Junction']==Junction) & \
                            (finalset['Bait_end']>=Bait_end-2) & \
                            (finalset['Bait_end']<=Bait_end+2)
            if de_method == '5':
                de_condition2 = (finalset['Barcode']==Barcode)
            
            # check_barcode = finalset[de_condition2]
            # close_barcode = get_close_matches(Barcode, check_barcode['Barcode'], n=100, cutoff=0.6)
            # print(close_barcode)
            # de_barcode = (finalset['Barcode'].isin(close_barcode))
            de_barcode = (finalset['Barcode'].str.contains(Barcode[0:15]))
            de_barcode1 =(finalset['Rname']==Rname) & \
                         (finalset['Strand']==Strand) & \
                         (finalset['Prey_start']==Prey_start) & \
                         (finalset['Prey_end']==Prey_end) & \
                         (finalset['Bait_end']==Bait_end)
            
            finalset.drop(finalset[(de_condition2)].index, inplace=True)
            # finalset.drop(finalset[(de_barcode1) & (de_barcode)].index, inplace=True)
            # print(finalset['Qname'].count())
            
            # finalset = finalset.drop_duplicates(['Rname','Strand','Barcode'],keep= 'first')
            ######
            
            add = add.append(reads.loc[i])
    finalset = finalset.append([add])
    finalset.reset_index(level=finalset.index.names, inplace=True)
    
    #Barcode filter
    finalset2 = finalset.copy()
    i_list = []
    f_list = []
    for i in range(0,finalset2['Qname'].count()):
        if i in i_list:
            continue
        Qname = finalset2['Qname'][i]
        Rname = finalset2['Rname'][i]
        Strand = finalset2['Strand'][i]
        Junction = finalset2['Junction'][i]
        Sequence = finalset2['Sequence'][i]
        Seq_len = finalset2['Seq_len'][i]
        Barcode = finalset2['Barcode'][i]

        #add barcode dedup method
        de_barcode_method = (finalset2['Barcode']==Barcode) & \
                            (finalset2['Rname']==Rname) & \
                            (finalset2['Strand']==Strand) & \
                            (finalset2['Junction']>=Junction-5) & \
                            (finalset2['Junction']<=Junction+5)

        add_list = finalset2.index[de_barcode_method].tolist()
        f_list = list(set(add_list) - set([i]))
        i_list = i_list + f_list
    # print(i_list)
    finalset.drop(index=i_list, inplace=True)
        
    # finalset.reset_index(level=finalset.index.names, inplace=True)
    # add.reset_index(level=add.index.names, inplace=True)
    
    # split barcode into single base
    length = finalset['Barcode'].str.len()
    finalset['Length'] = length
    bl = max(length)
    barcode_list = list(map(str, range(1,(bl+1))))
    x = range(0,bl)
    y = range(0,bl)
    couple = zip(x,y)
    for i,j in couple:
        finalset.loc[:,barcode_list[i]] = finalset['Barcode'].str[j]
    # get barcode length
    length = finalset['Barcode'].astype('str')
    length = finalset['Barcode'].str.len()
    finalset['Length'] = length
    finalset = finalset.sort_values(by=['Rname','Junction','Length'], ascending=False)
    finalset = finalset.drop_duplicates(['Rname','Strand','Prey_start','Prey_end','1','2','3','4','5',\
                                        '6','7','8','9','10','11','12','13','14','15'],keep='first')
    
    
    print(Full['Qname'].count())
    removed = Full[~Full['Qname'].isin(finalset['Qname'])]
    removed.to_csv("results/" + basename + yourfile + "_dsbremove.tab", header = True, sep = '\t', index=False)
    
    finalset = finalset[~finalset["Rname"].str.contains("_")]
    finalset = finalset.drop_duplicates(['Qname'],keep= 'first')
    finalset = finalset.sort_values(by=['Rname','Junction'], ascending=False)
    # finalset.reset_index(level=finalset.index.names, inplace=True)
    finalset.to_csv("results/" + basename + yourfile + "_dsb.tab", header = True, sep = '\t', index=False)
    # print(finalset[finalset['Rname'] == "chr15"]['Qname'].count())
    # add.to_csv(basename + "_transloc_dup_keeped.tab", header = True, sep = '\t', index=False)
    # reads.to_csv(basename + "_transloc_reads.tab", header = True, sep = '\t', index=False)
    # print(add['Qname'].count())
    # print(finalset['Qname'].count())
    
    # ### plot html file
    # check_range = (finalset['Rname'] == "chr2")
    # check = finalset[check_range]
    # check = check.sort_values(by=[u'Prey_end'], ascending=True)
    # check.to_csv(basename + yourfile + "_check3f.tab", header = True, sep = '\t', index=False, columns = [u'Qname',u'Bait_rname',u'Bait_start',u'Bait_end',
    #                                                                                                     u'Prey_rname',u'Prey_start',u'Prey_end',u'Insertion'])
    # check_range = (removed['Rname'] == "chr2")
    # checkremoved = removed[check_range]
    # checkremoved = checkremoved.sort_values(by=[u'Prey_end'], ascending=True)
    # checkremoved.to_csv(basename + yourfile + "_checkremoved3f.tab", header = True, sep = '\t', index=False, columns = [u'Qname',u'Bait_rname',u'Bait_start',u'Bait_end',
                                                                                                       # u'Prey_rname',u'Prey_start',u'Prey_end',u'Insertion'])
    
def main():
    start_time = time()
    
    args = docopt(__doc__,version='DSB_filter 1.0')
    
    kwargs = {'basename':args['<basename>'],'cutsite':args['<cutsite>'],'chromosome':args['<chromosome>'],'de_method':args['<de_method>'],'strand':args['<strand>'],'yourfile':args['<yourfile>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] chromosome: ' + str(kwargs['chromosome']))
    print('[PEM-Q] strand: ' + str(kwargs['strand']))
    print('[PEM-Q] de_method: ' + str(kwargs['de_method']))
    print('[PEM-Q] yourfile: ' + str(kwargs['yourfile']))
    
    DSB_filter(**kwargs)
    
    print("\nDSB_filter.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()