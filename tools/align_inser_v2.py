#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.15

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
    align_inser <basename>

Options:
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

def load(inputfile):
    
    if not os.path.exists(inputfile):
        raise ValueError('[PEM-Q] The {} file does not exist.'.format(inputfile))

def align_merged_fq_to_vector(basename):

    vector_fa = "vector/vector.fa"
    merged_fq = "flash_out/" + basename + ".merge.fastq.gz"
    load(vector_fa)
    load(merged_fq)
    
    merge_sam = basename + '_merge_vector.sam'
    merge_bam = basename + '_merge_vector.bam'
    merge_bam_sort = basename + '_merge_vector.sort.bam'
    
    os.system("mkdir vector")
    #~~~~~~~~align merged_fq to vector~~~~~~~~~~#

    #build bwa index of vector
    cmd = "samtools faidx vector/vector.fa"
    os.system(cmd)
    cmd = "bwa index -a bwtsw -p vector/vector vector/vector.fa 1>vector/build_index.o 2>vector/build_index.e"
    os.system(cmd)

    #alignment
    print("[PEM-Q]  align merged_fq to vector...")
                  
    cmd = "bwa mem -t 8 vector/vector -K 20  -L 0 -T 15 {} > vector/{} 2>vector/bwa_align_vector.log".format(merged_fq, 
                                                      merge_sam)
    print(cmd )
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(merge_sam, merge_bam, merge_bam, merge_bam_sort , merge_bam_sort)
    print("[PEM-Q]  sort and index bam...")
    os.system(cmd) 

def rmb_dedup(basename):
    
    merge_bam_sort = "vector/" + basename + '_merge_vector.sort.bam'
    merge_bam_sort_dedup = "vector/" + basename + '_merge_vector_dedup.bam'
    merge_bam_sort_dedup_sort = "vector/" + basename + '_merge_vector_dedup.sort.bam'
    
    unique_barcode = "barcode/"+basename+"_barcode_uniq.txt"
    load(unique_barcode)
    # get uniq qname list
    data = pd.read_csv(unique_barcode,sep = '\t',names = [u'Qname',u'Barcode',u'Freq',u'Length'],low_memory=False)
    uniq_qname_list = data['Qname'].tolist()
    
    #generate dedup bam
    merge_bam = pysam.AlignmentFile(merge_bam_sort, 'rb')
    dedup_bam = pysam.AlignmentFile(merge_bam_sort_dedup, "wb", template=merge_bam)
    
    #index bam name by pysam to generate dedup bam
    name_indexed = pysam.IndexedReads(merge_bam)
    name_indexed.build()
    for name in uniq_qname_list:
            try:
                name_indexed.find(name)
            except KeyError:
                pass
            else:
                iterator = name_indexed.find(name)
                for x in iterator:
                    dedup_bam.write(x)
    merge_bam.close()
    dedup_bam.close()
    pysam.sort("-o", merge_bam_sort_dedup_sort, merge_bam_sort_dedup)

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
    
def generate_vector_align_tab(basename):
    print("[PEM-Q]  generating vector tab file...")

    vector_seq = SeqIO.read('vector/vector.fa', "fasta").seq

    vector_bam = pysam.AlignmentFile("vector/" + basename + '_merge_vector_dedup.sort.bam', "rb")
    vector_tab = open("vector/" + basename + "_vector.tab", "w")
    vector_tab.write('Qname'+"\t"+\
                        'Vector'+"\t"+\
                        'Strand'+"\t"+\
                        'Start'+"\t"+\
                        'End'+"\t"+\
                        'Junction'+"\t"+\
                        'MH_end'+"\n")

    for read in vector_bam:
        
        if read.reference_name == "vector":
            Start = read.reference_start + 1
            End = read.reference_end
            Vector = vector_seq[Start-1:End]
            MH_end = find_MH_end(Vector)
            if read.is_reverse:
                Strand = '-'
                Junction = End
            else:
                Strand = '+'
                Junction = Start
            vector_tab.write(read.query_name+"\t"+\
                                str(Vector)+"\t"+\
                                str(Strand)+"\t"+\
                                str(Start)+"\t"+\
                                str(End)+"\t"+\
                                str(Junction)+"\t"+\
                                str(MH_end)+"\n")
    
    # #where come from
    # vector_tab_file = pd.read_csv("vector/" + basename + "_vector.tab",sep = '\t')
    # vector_list = vector_tab_file['Qname'].tolist()
    # discard_file = "indel/"+basename+"_discard.tab"
    # sid_file = pd.read_csv("YJ049a_SID_all.tab",sep = '\t')
    # sid_list = sid_file['Qname'].tolist()
    # insertion_vector_file = pd.read_csv(basename+"_vector.tab",sep = '\t')
    # insertion_vector_list = insertion_vector_file['Qname'].tolist()
    # with open(discard_file, "r") as f:
    #    discard = f.read().splitlines()
    # n = 0
    # m = 0
    # j = 0
    # for i in vector_list[0:len(vector_list)]:
    #     if str(i) in sid_list:
    #         j = j + 1
    # print("sid:   "+ str(j))
    # for i in vector_list[0:len(vector_list)]:
    #     if str(i) in discard:
    #         n = n + 1
    #     else:
    #         if str(i) in insertion_vector_list:
    #             m = m + 1
    #         else:
    #             print(i)
    # print("discard:   "+ str(n))
    # print("insertion_vector:   "+ str(m))
    
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='align_inser 2.0')
    
    kwargs = {'basename':args['<basename>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))

    align_merged_fq_to_vector(**kwargs)
    rmb_dedup(**kwargs)
    generate_vector_align_tab(**kwargs)
    
    print("\nalign_inser_v2.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
