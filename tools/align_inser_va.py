#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.8.16

"""align_inser

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
    align_inser <basename> [-i <vector_fa>]

Options:
-i <vector_fa>
-h --help               Show this screen.
-v --version            Show version.

This script help you align insetions sequences to vector.

Input file: transloc/basename_transloc.tab
            vector/vector.fa
Output file: alignment sam/bam file

Author: Mengzhu LIU
Last Update:2019.10.15

"""

import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO

def align_inser_to_vector(basename=None,vector_fa=None):
    fa_index = vector_fa.split(".")[0]
    # fa_index = "px330_sp_myc1"
    insertion = basename + "_insertion.fa"
    insertion_sam = basename + '_inser_vector.sam'
    insertion_bam = basename + '_inser_vector.bam'
    insertion_bam_sort = basename + '_inser_vector.sort.bam'
    
    #~~~~~~~~generate fasta file of insertions~~~~~~~~~~#

    insertion_file = "unique/" + basename + "_SID_all.tab"
    data = pd.read_csv(insertion_file, sep = '\t', index_col=False, low_memory=False)
    datb = data[data["Insertion"].str.len() > 0]
    # datb.to_csv(basename+"_transloc_insertion_all.tab",header=True,sep='\t',index=False)

    insertion_f = open(insertion, "w")
    print(len(data['Qname']))
    n = 0
    for i in range(0,len(data['Qname'])):
        if data['Insertion'][i] is np.nan:            
            continue
        else:
            n = n+1
            # print(data['Insertion'][i])
            insertion_f.write(">"+data['Qname'][i]+"\n")
            insertion_f.write(data['Insertion'][i]+"\n")
    insertion_f.close()
    print(n)
    
    #~~~~~~~~alignment insertions/unalign to vector~~~~~~~~~~#

    # #build bwa index of vector
    # cmd = "samtools faidx vector/vector.fa"
    # os.system(cmd)
    # cmd = "bwa index -a bwtsw -p vector/vector vector/vector.fa 1>vector/build_index.o 2>vector/build_index.e"
    # os.system(cmd)

    #alignment insertions
    print("[PEM-Q]  align insertions to vector...")
                  
    cmd = "bwa mem -Y -t 4 vector/{} {} > vector/{} 2>vector/bwa_align_vector.log".format(fa_index,insertion, 
                                                      insertion_sam)
    print(cmd )
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(insertion_sam, insertion_bam, insertion_bam, insertion_bam_sort , insertion_bam_sort)
    print("[PEM-Q]  sort and index bam...")
    os.system(cmd) 

def find_MH_end(Vector):
    
    Seq = Vector
    MH_end = ''
    n = 1
    while Seq[0:n] == Seq[-n:]:
        MH_end = Seq[0:n]
        n = n+1
        if n >len(Seq):
            break
           
    return(MH_end)
        
def generate_vector_align_tab(basename=None,vector_fa=None):
    print("[PEM-Q]  generating vector tab file...")
    fa_index = vector_fa.split(".")[0]
    # fa_index = "px330_sp_myc1"
    vector_seq = SeqIO.read('vector/'+fa_index+'.fa', "fasta").seq
    # vector_seq = SeqIO.read('PE2.fa', "fasta").seq
    
    insertion_bam = pysam.AlignmentFile("vector/" + basename + '_inser_vector.sort.bam', "rb")
    # insertion_bam = pysam.AlignmentFile(basename + '.sort.bam', "rb")
    vector_tab = open(basename + "_vector_confident_inser.tab", "w")
    vector_tab.write('Qname'+"\t"+\
                        'Vector_start'+"\t"+\
                        'Vector_end'+"\t"+\
                        'Vector_strand'+"\t"+\
                        'Align_sequence'+"\t"+\
                        'Align_sequence_R2'+"\t"+\
                        'Vector_inser_size'+"\t"+\
                        'Type'+"\n")

    for read in insertion_bam:
        if read.reference_name == "vector":
        # if read.reference_name == "PE2":
            pair_flag = "Confident"
            Align_sequence_R2 = ""
            Start = read.reference_start + 1
            End = read.reference_end
            size = End - Start + 1
            Insertion_Seq = read.query_sequence
            Vector = vector_seq[Start-1:End]
            MH_end = find_MH_end(Vector)
            if read.is_reverse:
                Strand = '-'
                Junction = End
            else:
                Strand = '+'
                Junction = Start
                
            vector_tab.write(read.query_name+"\t"+\
                                str(Start)+"\t"+\
                                str(End)+"\t"+\
                                str(Strand)+"\t"+\
                                Insertion_Seq+"\t"+\
                                Align_sequence_R2+"\t"+\
                                str(size)+"\t"+
                                pair_flag+"\n")
    vector_tab.close()
    
    vector_tab_file = pd.read_csv(basename+"_vector_confident_inser.tab",sep = '\t', index_col=False)
    insertion_file = pd.read_csv("unique/" + basename + "_Insertion.tab",sep = '\t', index_col=False)
    vector_merge = pd.merge(vector_tab_file, insertion_file, on='Qname', how='inner')
    vector_merge.to_csv(basename + "_vector_confident_inser.tab", header = True, sep = '\t', index=False, columns = [u'Qname',
    u'Vector_start', u'Vector_end', u'Vector_strand',u'Vector_inser_size', u'Bait_start', u'Bait_end', u'Prey_rname',
    u'Prey_start', u'Prey_end', u"Align_sequence", u"Align_sequence_R2",u"Type",u'Barcode'])

def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='align_inser 1.0')
    
    kwargs = {'basename':args['<basename>'],'vector_fa':args['-i']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] vector_fa: ' + str(kwargs['vector_fa']))
    

    align_inser_to_vector(**kwargs)
    generate_vector_align_tab(**kwargs)
    # align_vector_to_genome(**kwargs)

    print("\nalign_inser.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
