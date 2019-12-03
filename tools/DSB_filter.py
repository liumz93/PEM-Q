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

def DSB_filter(basename,chromosome,cutsite,strand,de_method,yourfile):
    ### column to dedup : Bait_rname,Bait_strand,Bait_start,Bait_end,Prey_rname,Prey_strand,Prey_start,Prey_end,Rname,Strand,Junction,Sequence
    #total dataset
    transloc = pd.read_csv("unique/"+basename+yourfile+".tab",sep = '\t') 
    Full = pd.read_csv("unique/"+basename+yourfile+".tab",sep = '\t')
    if de_method == '1':
        de_list = ['Sequence']
    if de_method == '2':
        de_list = ['Bait_rname','Bait_strand','Bait_start','Bait_end','Prey_rname','Prey_strand','Prey_start','Prey_end']
    # if de_method == '22':
    #     de_list = ['Rname','Strand','Junction','Bait_end','Prey_end']
    if de_method == '3' or de_method == '3f' or de_method == '22':
        de_list = ['Rname','Strand','Junction']
    de_list = ['Rname','Strand','Junction']
    dup = transloc[transloc.duplicated(de_list,keep = False)]
    dup.to_csv(basename + "_dsb_dup.tab", header = True, sep = '\t', index=False)
    # cutsite_tab = dup[(c1) & (c2) & (c3)]
    print(dup['Qname'].count())
    #reads to be judge 
    reads = dup.loc[~dup.duplicated(de_list,keep = 'first'),:]
    print(reads['Qname'].count())
    #drop reads around cutsite because they are true dsb reads
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
    #set finalset
    finalset = transloc
    # finalset = dataset
    add = pd.DataFrame()
    for i in range(0,reads['Qname'].count()):
        
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
        
        start = reads['Junction'][i] - 5
        end = reads['Junction'][i] + 5
        
        if de_method == '1':
            de_condition1 = transloc['Sequence']==Sequence
        if de_method == '2':
            de_condition1 = (transloc['Prey_rname']==Prey_rname) & \
                        (transloc['Prey_strand']==Prey_strand) & \
                        (transloc['Prey_start']==Prey_start) & \
                        (transloc['Prey_end']==Prey_end) & \
                        (transloc['Bait_rname']==Bait_rname) & \
                        (transloc['Bait_strand']==Bait_strand) & \
                        (transloc['Bait_start']==Bait_start) & \
                        (transloc['Bait_end']==Bait_end)
        if de_method == '22':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Junction']==Junction) & \
                        (transloc['Bait_end']==Bait_end) & \
                        (transloc['Prey_end']==Prey_end)
                        # (transloc['Bait_start']==Bait_start) & \
        if de_method == '3':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Junction']==Junction)
        if de_method == '3f':
            de_condition1 = (transloc['Rname']==Rname) & \
                        (transloc['Strand']==Strand) & \
                        (transloc['Junction']==Junction) & \
                        (transloc['Bait_start']>=Bait_start-1) & \
                        (transloc['Bait_start']<=Bait_start+1) & \
                        (transloc['Bait_end']>=Bait_end-1) & \
                        (transloc['Bait_end']<=Bait_end+1)
        
        #dataset to check
        dataset = transloc.drop(transloc[de_condition1].index)
        
        #condition
        if Strand == '+':
            s = '-'
        if Strand == '-':
            s = '+'
        count = dataset[(dataset['Rname'] == Rname) & (dataset['Junction'] >= start) & (dataset['Junction'] <= end)]# & (dataset['Strand'] == s)
        condition = count['Qname'].count() <= 2
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
                de_condition2 = (finalset['Prey_rname']==Prey_rname) & \
                            (finalset['Prey_strand']==Prey_strand) & \
                            (finalset['Prey_start']==Prey_start) & \
                            (finalset['Prey_end']==Prey_end) & \
                            (finalset['Bait_rname']==Bait_rname) & \
                            (finalset['Bait_strand']==Bait_strand) & \
                            (finalset['Bait_start']==Bait_start) & \
                            (finalset['Bait_end']==Bait_end)
            if de_method == '22':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Junction']==Junction) & \
                            (finalset['Bait_end']==Bait_end) & \
                            (finalset['Prey_end']==Prey_end)
                            # (finalset['Bait_start']==Bait_start) & \
            if de_method == '3':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Junction']==Junction)
            if de_method == '3f':
                de_condition2 = (finalset['Rname']==Rname) & \
                            (finalset['Strand']==Strand) & \
                            (finalset['Junction']==Junction) & \
                            (transloc['Bait_start']>=Bait_start-1) & \
                            (transloc['Bait_start']<=Bait_start+1) & \
                            (transloc['Bait_end']>=Bait_end-1) & \
                            (transloc['Bait_end']<=Bait_end+1)
            # print(Junction)
            finalset.drop(finalset[de_condition2].index, inplace=True)
            # print(finalset['Qname'].count())
            add = add.append(reads.loc[i])
    
    
    finalset = finalset.append([add])
    # finalset.reset_index(level=finalset.index.names, inplace=True)
    # add.reset_index(level=add.index.names, inplace=True)
    print(Full['Qname'].count())
    removed = Full[~Full['Qname'].isin(finalset['Qname'])]
    removed.to_csv("unique/" + basename + yourfile + "_dsbremoved.tab", header = True, sep = '\t', index=False)
    finalset.to_csv("unique/" + basename + yourfile + "_dsb.tab", header = True, sep = '\t', index=False)
    # print(finalset[finalset['Rname'] == "chr15"]['Qname'].count())
    # add.to_csv(basename + "_transloc_dup_keeped.tab", header = True, sep = '\t', index=False)
    # reads.to_csv(basename + "_transloc_reads.tab", header = True, sep = '\t', index=False)
    # print(add['Qname'].count())
    # print(finalset['Qname'].count())
def main():
    start_time = time()
    
    args = docopt(__doc__,version='DSB_filter 1.0')
    
    kwargs = {'basename':args['<basename>'],'cutsite':args['<cutsite>'],'chromosome':args['<chromosome>'],'de_method':args['<de_method>'],'strand':args['<strand>'],'yourfile':args['<yourfile>']}
    print('[superCasQ] basename: ' + str(kwargs['basename']))
    print('[superCasQ] cutsite: ' + str(kwargs['cutsite']))
    print('[superCasQ] chromosome: ' + str(kwargs['chromosome']))
    print('[superCasQ] strand: ' + str(kwargs['strand']))
    print('[superCasQ] de_method: ' + str(kwargs['de_method']))
    print('[superCasQ] yourfile: ' + str(kwargs['yourfile']))
    
    DSB_filter(**kwargs)
    
    print("\nDSB_filter.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()