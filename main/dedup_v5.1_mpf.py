#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.21
"""rmb_dedup

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
    dedup   <basename>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script help to dedup reads according its barcode attatched.

Input file:  / 

Author: Mengzhu LIU
Last Update:2019.10.21

"""
import os
import pysam
from time import time
from docopt import docopt
import pandas as pd
import numpy as np

def dedup(basename=None):
    # barcode = pd.read_csv("barcode/"+basename+"_barcode_uniq.txt",sep = '\t',low_memory=False)
    barcode = pd.read_csv("barcode/"+basename+"_barcode_list.txt",sep = '\t',names = ["Qname", "Barcode"])
    sid = pd.read_csv("raw/"+basename+"_SID_all.tab",sep = '\t')
    sv = pd.read_csv("raw/"+basename+"_SV.tab",sep = '\t')
    substitution = pd.read_csv("raw/"+basename+"_Substitution.tab",sep = '\t')
    insertion = pd.read_csv("raw/"+basename+"_Insertion.tab",sep = '\t')
    deletion = pd.read_csv("raw/"+basename+"_Deletion.tab",sep = '\t')
    germline = pd.read_csv("raw/"+basename+"_Germline.tab",sep = '\t')
    transloc = pd.read_csv("transloc/"+basename+"_mut.tab",sep = '\t')
    indel = pd.read_csv("indel/"+basename+"_indel_mut.tab",sep = '\t')
    # all_insertion = pd.read_csv(basename+"_All_Insertions.tab",sep = '\t')
    germline_addup = pd.read_table("raw/" + basename + "_Germline_addup_indel.tab",sep = '\t')
    indel_cutoff = pd.read_table("raw/"+basename + "_indel_cutoff.tab",sep = '\t')
    insertion_sv = pd.read_csv("raw/"+basename+"_Insertions_inSV.tab",sep = '\t')

    sid_merge = pd.merge(sid, barcode, on='Qname', how='inner')
    sv_merge = pd.merge(sv, barcode, on='Qname', how='inner')
    substitution_merge = pd.merge(substitution, barcode, on='Qname', how='inner')
    insertion_merge = pd.merge(insertion, barcode, on='Qname', how='inner')
    deletion_merge = pd.merge(deletion, barcode, on='Qname', how='inner')
    germline_merge = pd.merge(germline, barcode, on='Qname', how='inner')
    transloc_merge = pd.merge(transloc, barcode, on='Qname', how='inner')
    indel_merge = pd.merge(indel, barcode, on='Qname', how='inner')
    germline_addup_merge = pd.merge(germline_addup, barcode, on='Qname', how='inner')
    indel_cutoff_merge = pd.merge(indel_cutoff, barcode, on='Qname', how='inner')
    insertion_sv_merge = pd.merge(insertion_sv, barcode, on='Qname', how='inner')
    # print(sid_merge['Qname'].count(),sv_merge['Qname'].count(),substitution_merge['Qname'].count(),insertion_merge['Qname'].count(),deletion_merge['Qname'].count(),germline_merge['Qname'].count())
    # all_insertion_merge = pd.merge(all_insertion, barcode, on='Qname', how='inner')


    sid_merge.to_csv("raw/"+basename + "_SID_all.tab", header = True, sep = '\t', index=False)
    sv_merge.to_csv("raw/"+basename + "_SV.tab", header = True, sep = '\t', index=False)
    insertion_merge.to_csv("raw/"+basename + "_Insertion.tab", header = True, sep = '\t', index=False)
    insertion_sv_merge.to_csv("raw/"+basename + "_Insertion_sv.tab", header = True, sep = '\t', index=False)
    deletion_merge.to_csv("raw/"+basename + "_Deletion.tab", header = True, sep = '\t', index=False)
    germline_merge.to_csv("raw/"+basename + "_Germline.tab", header = True, sep = '\t', index=False)
    germline_addup_merge.to_csv("raw/"+basename + "_Germline_addup_indel.tab", header = True, sep = '\t', index=False)
    substitution_merge.to_csv("raw/"+basename + "_Substitution.tab", header = True, sep = '\t', index=False)
    indel_cutoff_merge.to_csv("raw/"+basename + "_smallindel_cutoff.tab", header = True, sep = '\t', index=False)
    transloc_merge.to_csv("raw/"+basename + "_transloc.tab", header = True, sep = '\t', index=False)
    indel_merge.to_csv("raw/"+basename + "_smallindel.tab", header = True, sep = '\t', index=False)

    os.system("mkdir unique")

    germline_dedup = germline_merge.drop_duplicates(['Cigar','Barcode','Position'],keep= 'first')
    germline_dedup.to_csv("unique/"+basename + "_Germline.tab", header = True, sep = '\t', index=False)

    germline_addup_dedup = germline_addup_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Barcode'],keep='first')
    germline_addup_dedup.to_csv("unique/"+basename + "_Germline_addup_indel.tab", header = True, sep = '\t', index=False)

    substition_dedup = substitution_merge.drop_duplicates(['MDstring','Barcode'],keep='first')
    substition_dedup.to_csv("unique/"+basename + "_Substitution.tab", header = True, sep = '\t', index=False)
    
     
    # 'B_Qstart,B_Qend,Qstart,Qend,Qlen'
    # 'Rname,Strand,Bait_end,Junction'
    cmd = "repeats_dedup.py {}_SID_all.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_SV.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_Insertion.tab 30 -f 'Rname,Strand,Bait_end,Insertion,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_Insertion_sv.tab 30 -f 'Rname,Strand,Bait_end,Insertion,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_Deletion.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
#     cmd = "repeats_dedup.py {}_Germline.tab 30 -f 'Cigar,Position'".format(basename)
#     print(cmd)
#     os.system(cmd)
#     cmd = "repeats_dedup.py {}_Germline_addup_indel.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
#     print(cmd)
#     os.system(cmd)
#     cmd = "repeats_dedup.py {}_Substitution.tab 30 -f 'MDstring'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_smallindel_cutoff.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_transloc.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    cmd = "repeats_dedup.py {}_smallindel.tab 30 -f 'Rname,Strand,Bait_end,Junction'".format(basename)
    print(cmd)
    os.system(cmd)
    
    
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='dedup 1.0')
    
    kwargs = {'basename':args['<basename>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    
    dedup(**kwargs)

    print("\ndedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()

    