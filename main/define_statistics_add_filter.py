#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2021.01.04
"""define_substitution

    This is part of PEM-Q pipeline to analyze PEM-seq data or data similar, help you analyze repair outcome of your DNA library.

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
    define_statistics    <basename> <genome> <cutsite> <transloc_range> <primer> <primer_chr> <primer_strand> <adapter>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .

Input file: indel alignment bam file / Output file: a informative tab file

Author: Mengzhu LIU
Last Update:2021.01.14

"""


import os
import pysam
import re
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
    
def statistics_add_filter(basename=None,genome=None,cutsite=None,transloc_range=None,primer=None,primer_chr=None,primer_strand=None,adapter=None):
    
    
    print("Processing junction filter and statistics...")
    os.system('mkdir results')
    
    # ~~~~ Germline ~~~~ #
    
    germ_file = "unique/" + basename + "_Germline_final.tab"
    germ = pd.read_csv(germ_file, sep = '\t', index_col=False, low_memory=False)
    germ.to_csv("results/" + basename + "_Germline.tab", header = True, sep = '\t', index=False)
    germ_count = germ["Qname"].count()
    
    # ~~~~ Deletion  ~~~~ #
    
    deletion_file = "unique/" + basename+"_Deletion.tab"
    deletion_all = pd.read_csv(deletion_file, sep = '\t', index_col=False, low_memory=False)
    # print(deletion_all['Qname'].count())
    # Deletion misprime filter
    condition_misprime = deletion_all['Bait_end'] - deletion_all['Bait_start'] > len(primer) + 10
    deletion_all = deletion_all[condition_misprime]
    # print(deletion_all['Qname'].count())
    
    # Deletion (within transloc_range) 
    deletion_all['Microhomolog_len']=deletion_all['Microhomolog'].map(lambda x: 0 if x is np.nan else len(x))
    
    if primer_strand == "+":
        #deletion length
        deletion_all['deletion_length'] = deletion_all['Prey_start'] - deletion_all['Bait_end'] + deletion_all['Microhomolog_len'] - 1
        condition_del = deletion_all['deletion_length'] <= int(transloc_range)
    else:
        deletion_all['deletion_length'] = deletion_all['Bait_start'] - deletion_all['Prey_end'] + deletion_all['Microhomolog_len'] - 1
        condition_del = deletion_all['deletion_length'] <= int(transloc_range)
    deletion = deletion_all[condition_del]
    deletion.to_csv("results/" + basename+"_Deletion.tab", header = True, sep = '\t', index=False)
    deletion_count = deletion["Qname"].count()
    
    small_del = deletion[deletion['deletion_length']<=100]['Qname'].count()
    large_del = deletion[deletion['deletion_length']>100]['Qname'].count()
    # print(small_del,large_del)
    del_stats_file = open("results/" + basename+"_del_len_statistics.txt","w")
    del_stats_file.write("small_del(<=100bp)"+"\t"+str(small_del)+"\n")
    del_stats_file.write("large_del(>100bp)"+"\t"+str(large_del)+"\n")
    del_stats_file.close()
    
    deletion_length_freq = deletion['deletion_length'].value_counts()
    deletion_length_freq = deletion_length_freq.sort_index()
    deletion_length_freq.to_csv("results/" + basename+"_deletion_length.txt", header = True, sep = '\t', index=True)
    
    
    # ~~~~ Insertion ~~~~ #
    
    insertion_file = "unique/" + basename+"_Insertion.tab"
    insertion = pd.read_csv(insertion_file, sep = '\t', index_col=False, low_memory=False)
    # insertion misprime filter
    condition_misprime = insertion['Bait_end'] - insertion['Bait_start'] > len(primer) + 10
    insertion = insertion[condition_misprime]
    
    insertion.to_csv("results/" + basename + "_Insertion.tab", header = True, sep = '\t', index=False)
    insertion_count = insertion["Qname"].count()
    
    insertion['insertion_length']=insertion['Insertion'].map(lambda x: 0 if x is np.nan else len(x))
    insertion_length_freq = insertion['insertion_length'].value_counts() 
    insertion_length_freq = insertion_length_freq.sort_index()
    insertion_length_freq.to_csv("results/" + basename+"_insertion_length.txt", header = True, sep = '\t', index=True)
    
    one_bp_inser = insertion[insertion['insertion_length']==1]['Qname'].count()
    small_inser = insertion[insertion['insertion_length']<20]['Qname'].count()
    large_inser = insertion[insertion['insertion_length']>=20]['Qname'].count()
    # print(small_inser,large_inser)
    inser_stats_file = open("results/" + basename+"_inser_len_statistics.txt","w")
    inser_stats_file.write("1_bp"+"\t"+str(one_bp_inser)+"\n")
    inser_stats_file.write("small_inser(<20bp)"+"\t"+str(small_inser)+"\n")
    inser_stats_file.write("large_inser(>=20bp)"+"\t"+str(large_inser)+"\n")
    inser_stats_file.close()
    
    # ~~~~ Inversion ~~~~ #
    
    inverion_file = "unique/" + basename+"_invertion.tab"
    inverion = pd.read_csv(inverion_file, sep = '\t', index_col=False, low_memory=False)
    # inverion misprime filter
    condition_misprime = inverion['Bait_end'] - inverion['Bait_start'] > len(primer) + 10
    inverion = inverion[condition_misprime]
    
    inverion['Microhomolog_len']=inverion['Microhomolog'].map(lambda x: 0 if x is np.nan else len(x))
    inverion.to_csv("results/" + basename + "_Inversion.tab", header = True, sep = '\t', index=False)
    inverion_count = inverion["Qname"].count()
    if primer_strand == "+":
        condition_close_inver = abs(inverion['Prey_start'] - inverion['Bait_end']) + inverion['Microhomolog_len'] - 1 < int(transloc_range)
    else:
        condition_close_inver = abs(inverion['Bait_start'] - inverion['Prey_end']) + inverion['Microhomolog_len'] - 1 < int(transloc_range)
    close_inver = inverion[condition_close_inver]
    close_inver.to_csv("results/" + basename+"_close.inver.tab", header = True, sep = '\t', index=False)
    close_inver_count = close_inver["Qname"].count()
    
    # ~~~~ Translocation ~~~~ #
    
    # Part1. Intra_translocations_from_deletions
    
    if primer_strand == "+":
        condition_transloc = (deletion_all['Prey_start'] - deletion_all['Bait_end'] + deletion_all['Microhomolog_len'] - 1) > int(transloc_range)
    else:
        condition_transloc = (deletion_all['Bait_start'] - deletion_all['Prey_end'] + deletion_all['Microhomolog_len'] - 1) > int(transloc_range)
    transloc_del = deletion_all[condition_transloc]
    transloc_del.to_csv("results/" + basename+"_intra_Translocation.del.tab", header = True, sep = '\t', index=False)
    transloc_del_count = transloc_del["Qname"].count()
    
    # Part2. Intra_transloctions_from_inversion
    
    if primer_strand == "+":
        condition_inver = abs(inverion['Prey_start'] - inverion['Bait_end']) + inverion['Microhomolog_len'] - 1 > int(transloc_range)
    else:
        condition_inver = abs(inverion['Bait_start'] - inverion['Prey_end']) + inverion['Microhomolog_len'] - 1 > int(transloc_range)
    transloc_inver = inverion[condition_inver]
    transloc_inver.to_csv("results/" + basename+"_intra_Translocation.inver.tab", header = True, sep = '\t', index=False)
    transloc_inver_count = transloc_inver["Qname"].count()
    
    # Part3. Inter_translocations
    
    intersv_file = "unique/" + basename+"_SV.tab"
    intersv = pd.read_csv(intersv_file, sep = '\t', index_col=False, low_memory=False)
    intersv_con = (intersv["Prey_rname"] != primer_chr)
    intersv = intersv[intersv_con]
    
    # intersv misprime filter
    condition_misprime = intersv['Bait_end'] - intersv['Bait_start'] > len(primer) + 10
    intersv = intersv[condition_misprime]
    # # insertin(gap) filter
    # intersv['insertion_len']=intersv['Insertion'].map(lambda x: 0 if x is np.nan else len(x))
    # condition_insertion= intersv['insertion_len'] < 20
    # print(intersv['Qname'].count())
    # intersv = intersv[condition_insertion]
    # print(intersv['Qname'].count())
    
    
    intersv.to_csv("results/" + basename + "_inter_Translocation.tab", header = True, sep = '\t', index=False)
    intersv_count = intersv["Qname"].count()
    
    # Total translocations
    translocations = transloc_del.append([transloc_inver])
    translocations = translocations.append([intersv])
    translocations = translocations[['Qname','Bait_rname','Bait_strand','Bait_start','Bait_end','Prey_rname',\
                                    'Prey_strand','Prey_start','Prey_end','Rname','Strand','Junction','Sequence',\
                                    'B_Qstart','B_Qend','Qstart','Qend','Qlen','Prey_MQ','Insertion','Microhomolog','Barcode']]
    translocations.to_csv("results/" + basename + "_Translocation.tab", header = True, sep = '\t', index=False)
    
    # Translocation filter
    cmd = "DSB_filter_update.py {} {} {} {} {} {}".format(basename,primer_chr,cutsite,primer_strand,"3fe","_Translocation")
    print(cmd)
    os.system(cmd)

    # translocations = pd.read_csv("results/" + basename + "_Translocation_dsb.tab", sep = '\t', index_col=False, low_memory=False)
    # translocations = translocations[['Qname','Bait_rname','Bait_strand','Bait_start','Bait_end','Prey_rname',\
    #                                 'Prey_strand','Prey_start','Prey_end','Rname','Strand','Junction','Sequence',\
    #                                 'B_Qstart','B_Qend','Qstart','Qend','Qlen','Insertion','Microhomolog','Barcode']]
    #
    # translocations.to_csv("results/" + basename + "_Translocation_filtered.tab", header = True, sep = '\t', index=False)
    translocations_count = translocations["Qname"].count()
    
    total_events = germ_count + deletion_count + insertion_count + transloc_del_count + inverion_count + intersv_count
    edit_events = deletion_count + insertion_count + transloc_del_count + transloc_inver_count + intersv_count
    edit_effi = edit_events/total_events
    
    stats_file = open("results/" + basename+"_statistics.txt","w")
    stats_file.write("NoJunction"+"\t"+str(germ_count)+"\n")
    stats_file.write("Deletion"+"\t"+str(deletion_count)+"\n")
    stats_file.write("Insertion"+"\t"+str(insertion_count)+"\n")
    stats_file.write("Intra_Translocation.del"+"\t"+str(transloc_del_count)+"\n")
    stats_file.write("Close_inversion.del"+"\t"+str(close_inver_count)+"\n")
    stats_file.write("Intra_Translocation.inver"+"\t"+str(transloc_inver_count)+"\n")
    stats_file.write("Inter_Translocation"+"\t"+str(intersv_count)+"\n")
    stats_file.write("Translocation"+"\t"+str(translocations_count)+"\n")
    stats_file.write("Editing Events"+"\t"+str(edit_events)+"\n")
    stats_file.write("Total Events"+"\t"+str(total_events)+"\n")
    stats_file.write("Editing Efficiency"+"\t"+str(edit_effi)+"\n")
    stats_file.close()
    
    # Merge editing events
    deletion.drop(columns=['Microhomolog_len','deletion_length'], inplace=True)
    insertion.drop(columns=['insertion_length'], inplace=True)
    close_inver.drop(columns=['Microhomolog_len'], inplace=True)
    #concatenate dataframes
    frames = [deletion, insertion, translocations, close_inver]
    merge_df = pd.concat(frames, sort=False)
    merge_df = merge_df.sort_values(by=['Rname','Junction','Insertion','Barcode'], ascending=True)
    merge_df.reset_index(drop=True, inplace=True)
    merge_df.to_csv("results/" + basename + "_Editing_events.tab", header = True, sep = '\t', index=False)
    
    #Plot dot plot
    cmd = "TranslocPlot.R results/{}_Editing_events.tab results/{}_Editing_events_dot_plot.pdf \
           binsize=2000000 assembly={} plotshape=octogon strand=0".format(basename,basename,genome)
    print(cmd)
    os.system(cmd)
    
    #generate html file
    cmd = "TranslocHTMLReads_PEMQ.pl results/{}_Editing_events.tab results/{}_Editing_events.html \
           --primer {} -adapter {}".format(basename,basename,primer,adapter)
    print(cmd)
    os.system(cmd)
    
def main():

    start_time = time()

    args = docopt(__doc__,version='define_substitution 1.0')

    kwargs = {'basename':args['<basename>'],'genome':args['<genome>'],'cutsite':args['<cutsite>'],'transloc_range':args['<transloc_range>'],\
              'primer':args['<primer>'],'primer_chr':args['<primer_chr>'],'primer_strand':args['<primer_strand>'],'adapter':args['<adapter>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] transloc_range: ' + str(kwargs['transloc_range']))
    print('[PEM-Q] primer: ' + str(kwargs['primer']))
    print('[PEM-Q] primer_chr: ' + str(kwargs['primer_chr']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    print('[PEM-Q] adapter: ' + str(kwargs['adapter']))
    
    
    statistics_add_filter(**kwargs)

    print("\ndefine_substitution.py Done in {}s".format(round(time()-start_time, 3)))

if __name__ == '__main__':
    main()
