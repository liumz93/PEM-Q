#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.18
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
    define_substitution    <basename>  <cutsite>  <cutoff> <primer_strand>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .

Input file: indel alignment bam file / Output file: a informative tab file

Author: Mengzhu LIU
Last Update:2019.6.10

"""


import os
import pysam
import re
from time import time
from docopt import docopt
import pandas as pd
from Bio import SearchIO
import numpy as np

def cal_soft_clipping_number(cigar):
    
    cigar_list = re.split('(\d+)',cigar)
    if cigar_list[2] == 'S':
        number = int(cigar_list[1])
    else:
        number = 0
    return(number)
    
def define_substitution(basename=None,cutsite=None,cutoff=None,primer_strand=None):
    
    substitution_file = "unique/" + basename + "_Substitution.tab"
    data = pd.read_csv(substitution_file, sep = '\t', index_col=False, low_memory=False)
    data_cutoff = pd.DataFrame()
    germline_addup_st = pd.DataFrame()
    
    mdstring = data["MDstring"]
    count = data["MDstring"].count()
    n = 0
    for i in range(0, count):
        md_string = data["MDstring"][i]
        md_list = re.split('(\d+)',md_string)
        md_string_list = md_list[1:(len(md_list)-1)]
        len_list = []
        seq_length = 0
        sub_length = 0
        for j in range(0, len(md_string_list)):
            if md_string_list[j].isdigit():
                len_list.append(int(md_string_list[j]))
                seq_length = sum(len_list)
            else:
                sub_length = sub_length + 1
                seq_length = seq_length + sub_length
                bed_start = (int(data["Position"][i]) + seq_length - 1) - 1
                bed_end = bed_start + 1#real position
            
                clip_number=cal_soft_clipping_number(data["Cigar"][i])
                mismatch_base = data["Sequence"][i][(clip_number+seq_length-1):(clip_number+seq_length)]
                # print(bed_start,bed_end)
                if bed_start >= (int(cutsite) - int(cutoff)) and bed_end <= (int(cutsite) + int(cutoff) - 1):
                    define = "Y"
                    n = n + 1
                    # print([data.loc[[i]]])
                    data_cutoff = data_cutoff.append([data.loc[[i]]])
                    break
                   
    data_cutoff.to_csv("unique/" + basename + "_Substitution_cutoff.tab", header = True, sep = '\t', index=False)
    print("Substitutions in cutsite +- cutoff:",n)
    if n == 0:
        germline_addup_st = data
    else:
        germline_addup_st = data[~data['Qname'].isin(data_cutoff['Qname'])]
    # germline_addup_st = data.append([data_cutoff])
    germline_addup_st = germline_addup_st.drop_duplicates(keep=False)
    germline_addup_st.to_csv("unique/" + basename + "_Germline_addup_st.tab", header = True, sep = '\t', index=False)
    print(germline_addup_st["Qname"].count())
    
    #statistics
    print("Processing statistics...")
    deletion_file = "unique/" + basename+"_Deletion.tab"
    data = pd.read_csv(deletion_file, sep = '\t', index_col=False, low_memory=False)
    deletion_count = data["Qname"].count()
    
    insertion_file = "unique/" + basename+"_Insertion.tab"
    data = pd.read_csv(insertion_file, sep = '\t', index_col=False, low_memory=False)
    insertion_count = data["Qname"].count()
    
    sv_file = "unique/" + basename+"_SV.tab"
    data = pd.read_csv(sv_file, sep = '\t', index_col=False, low_memory=False)
    sv_count = data["Qname"].count()
    
    #invertions
    chrom = data["Bait_rname"][0]
    invertions_con = (data["Prey_rname"] == chrom)
    invertions = data[invertions_con]
    invertions.to_csv("unique/" + basename+"_invertion.tab", header = True, sep = '\t', index=False)
    invertions_count = invertions["Qname"].count()
    #translocations(inter)
    interTransloc_con = (data["Prey_rname"] != chrom)
    interTransloc = data[interTransloc_con]
    interTransloc.to_csv("unique/" + basename+"_interSV.tab", header = True, sep = '\t', index=False)
    interTransloc_count = interTransloc["Qname"].count()
    
    germ_file = "unique/" + basename+"_Germline.tab"
    germ_addup_st = pd.read_csv("unique/" + basename+"_Germline_addup_st.tab", sep = '\t', index_col=False, low_memory=False)
    germ_addup_indel = pd.read_csv("unique/" + basename+"_Germline_addup_indel.tab", sep = '\t', index_col=False, low_memory=False)
    data = pd.read_csv(germ_file, sep = '\t', index_col=False, low_memory=False)
    data = data.append([germ_addup_st])
    data = data.append([germ_addup_indel])
    data.to_csv("unique/" + basename + "_Germline_final.tab", header = True, sep = '\t', index=False)
    germ_count = data["Qname"].count()
    
    # total_events = n + deletion_count + insertion_count + sv_count + germ_count
    # edit_events = n + deletion_count + insertion_count + sv_count
    # edit_effi = edit_events/total_events
    #
    # stats_file = open(basename+"_stats.txt","w")
    # stats_file.write("Substitution:"+str(n)+"\n")
    # stats_file.write("Deletion:"+str(deletion_count)+"\n")
    # stats_file.write("Insertion:"+str(insertion_count)+"\n")
    # stats_file.write("Structure Variation(SV):"+str(sv_count)+"\n")
    # stats_file.write("Inversion:"+str(invertions_count)+"\n")
    # stats_file.write("interSV:"+str(interTransloc_count)+"\n")
    # stats_file.write("Germline:"+str(germ_count)+"\n")
    # stats_file.write("Total Events:"+str(total_events)+"\n")
    # stats_file.write("Editing Events:"+str(edit_events)+"\n")
    # stats_file.write("Editing Efficiency:"+str(edit_effi)+"\n")
    # stats_file.close()
    
def main():

    start_time = time()

    args = docopt(__doc__,version='define_substitution 1.0')

    kwargs = {'basename':args['<basename>'],'cutsite':args['<cutsite>'],'cutoff':args['<cutoff>'],'primer_strand':args['<primer_strand>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] cutoff: ' + str(kwargs['cutoff']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    
    
    define_substitution(**kwargs)

    print("\ndefine_substitution.py Done in {}s".format(round(time()-start_time, 3)))

if __name__ == '__main__':
    main()
