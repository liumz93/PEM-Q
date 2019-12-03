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

def dedup(basename=None):
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
    
    sid_merge = pd.merge(sid, barcode, on='Qname', how='inner')
    sv_merge = pd.merge(sv, barcode, on='Qname', how='inner')
    substitution_merge = pd.merge(substitution, barcode, on='Qname', how='inner')
    insertion_merge = pd.merge(insertion, barcode, on='Qname', how='inner')
    deletion_merge = pd.merge(deletion, barcode, on='Qname', how='inner')
    germline_merge = pd.merge(germline, barcode, on='Qname', how='inner')
    transloc_merge = pd.merge(transloc, barcode, on='Qname', how='inner')
    indel_merge = pd.merge(indel, barcode, on='Qname', how='inner')
    # print(sid_merge['Qname'].count(),sv_merge['Qname'].count(),substitution_merge['Qname'].count(),insertion_merge['Qname'].count(),deletion_merge['Qname'].count(),germline_merge['Qname'].count())
    # all_insertion_merge = pd.merge(all_insertion, barcode, on='Qname', how='inner')
    
    sid_dedup = sid_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Barcode'],keep='first')
    sv_dedup = sv_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Barcode'],keep='first')
    insertion_dedup = insertion_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Insertion','Barcode'],keep='first')
    deletion_dedup = deletion_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Microhomolog','Barcode'],keep='first')
    
    germline_dedup = germline_merge.drop_duplicates(['Cigar','Barcode','Position'],keep= 'first')
    substition_dedup = substitution_merge.drop_duplicates(['MDstring','Barcode'],keep='first')
    
    transloc_dedup = transloc_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Barcode'],keep= 'first')
    indel_dedup = indel_merge.drop_duplicates(['Rname','Strand','Bait_end','Junction','Barcode'],keep='first')
    
    os.system("mkdir unique")
    sid_dedup.to_csv("unique/"+basename + "_SID_all.tab", header = True, sep = '\t', index=False)
    sv_dedup.to_csv("unique/"+basename + "_SV.tab", header = True, sep = '\t', index=False)
    insertion_dedup.to_csv("unique/"+basename + "_Insertion.tab", header = True, sep = '\t', index=False)
    deletion_dedup.to_csv("unique/"+basename + "_Deletion.tab", header = True, sep = '\t', index=False)
    germline_dedup.to_csv("unique/"+basename + "_Germline.tab", header = True, sep = '\t', index=False)
    substition_dedup.to_csv("unique/"+basename + "_Substitution.tab", header = True, sep = '\t', index=False)
    
    transloc_dedup.to_csv("unique/"+basename + "_transloc.tab", header = True, sep = '\t', index=False)
    indel_dedup.to_csv("unique/"+basename + "_smallindel.tab", header = True, sep = '\t', index=False)

def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='dedup 1.0')
    
    kwargs = {'basename':args['<basename>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    
    dedup(**kwargs)

    print("\ndedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()

    