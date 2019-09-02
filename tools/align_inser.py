#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.8.16

"""align_inser

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
Last Update:2019.5.8

"""

import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO

def align_inser_to_vector(basename):

    insertion = basename + "_insertion.fa"
    insertion_sam = basename + '_inser_vector.sam'
    insertion_bam = basename + '_inser_vector.bam'
    insertion_bam_sort = basename + '_inser_vector.sort.bam'
    inser_genom_sam = basename + '_inser_genom.sam'
    inser_genom_bam = basename + '_inser_genom.bam'
    inser_genom_bam_sort = basename + '_inser_genom.sort.bam'
    
    #~~~~~~~~generate fasta file of insertions~~~~~~~~~~#

    transloc_file = basename + "_transloc_all.tab"
    data = pd.read_csv(transloc_file, sep = '\t', index_col=False, low_memory=False)
    datb = data[data["Insertion"].str.len() > 0]
    datb.to_csv(basename+"_transloc_insertion_all.tab",header=True,sep='\t',index=False) 

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
    
    #~~~~~~~~alignment insertions to vector~~~~~~~~~~#

    #build bwa index of vector
    cmd = "samtools faidx vector/vector.fa"
    os.system(cmd)
    cmd = "bwa index -a bwtsw -p vector/vector vector/vector.fa 1>vector/build_index.o 2>vector/build_index.e"
    os.system(cmd)

    #alignment
    print("[superCasQ]  align insertions to vector...")
                  
    cmd = "bwa mem -t 4 vector/vector -k 10 -L 0 -T 10 {} > bwa_align/{} 2>bwa_align/bwa_align_vector.log".format(insertion, 
                                                      insertion_sam)
    print(cmd )
    os.system(cmd)

    cmd = "samtools view -S -b -h bwa_align/{} > bwa_align/{} \
           && samtools sort bwa_align/{} > bwa_align/{} \
           && samtools index bwa_align/{}".format(insertion_sam, insertion_bam, insertion_bam, insertion_bam_sort , insertion_bam_sort)
    print("[superCasQ]  sort and index bam...")
    os.system(cmd) 
    
    #~~~~~~~~alignment insertions to genome~~~~~~~~~~#
    
    bwa_index_path = "~/database/bwa_indexes/{}/{}".format("hg38", "hg38")
    print("[superCasQ] index file used "+bwa_index_path)
    print("[superCasQ] align insertions to genome...")
                      
    cmd = "bwa mem -Y -t 4 {} {} > bwa_align/{} 2>bwa_align/bwa_align_inser_genom.log".format(bwa_index_path, 
                                         insertion,
                                         inser_genom_sam)
    os.system(cmd)
    
    cmd = "samtools view -S -b -h bwa_align/{} > bwa_align/{} \
           && samtools sort bwa_align/{} > bwa_align/{} \
           && samtools index bwa_align/{}".format(inser_genom_sam, 
                                        inser_genom_bam, 
                                        inser_genom_bam, 
                                        inser_genom_bam_sort, 
                                        inser_genom_bam_sort)

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
    print("[superCasQ]  generating vector tab file...")
    
    vector_seq = SeqIO.read('vector/vector.fa', "fasta").seq
    
    insertion_bam = pysam.AlignmentFile("bwa_align/" + basename + '_inser_vector.sort.bam', "rb")
    vector_tab = open(basename + "_vector.tab", "w")
    vector_tab.write('Qname'+"\t"+\
                        'Vector'+"\t"+\
                        'Strand'+"\t"+\
                        'Start'+"\t"+\
                        'End'+"\t"+\
                        'Junction'+"\t"+\
                        'MH_end'+"\n")

    for read in insertion_bam:
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
    vector_tab.close()
    
def align_vector_to_genome(basename):
    
    vector_fa = basename + "_vector.fa"
    vector_genom_sam = basename + '_vector_genom.sam'
    vector_genom_bam = basename + '_vector_genom.bam'
    vector_genom_bam_sort = basename + '_vector_genom.sort.bam'
    
    #~~~~~~~~generate fasta file of vectors~~~~~~~~~~#

    vector_file = basename + "_vector.tab"
    data = pd.read_csv(vector_file, sep = '\t', index_col=False, low_memory=False)

    vector_f = open(vector_fa, "w")
    print(len(data['Qname']))
    n = 0
    for i in range(0,len(data['Qname'])):
        if data['Vector'][i] is np.nan:            
            continue
        else:
            n = n+1
            # print(data['Insertion'][i])
            vector_f.write(">"+data['Qname'][i]+"\n")
            vector_f.write(data['Vector'][i]+"\n")
    vector_f.close()
    print(n)
    
    #~~~~~~~~alignment vector to genome~~~~~~~~~~#
    
    bwa_index_path = "~/database/bwa_indexes/{}/{}".format("hg38", "hg38")
    print("[superCasQ] index file used "+bwa_index_path)
    print("[superCasQ] align vector to genome...")
                      
    cmd = "bwa mem -Y -t 4 {} {} > bwa_align/{} 2>bwa_align/bwa_align_vector_genom.log".format(bwa_index_path, 
                                         vector_fa,
                                         vector_genom_sam)
    os.system(cmd)
    
    cmd = "samtools view -S -b -h bwa_align/{} > bwa_align/{} \
           && samtools sort bwa_align/{} > bwa_align/{} \
           && samtools index bwa_align/{}".format(vector_genom_sam, 
                                        vector_genom_bam, 
                                        vector_genom_bam, 
                                        vector_genom_bam_sort, 
                                        vector_genom_bam_sort)
    
    os.system(cmd)
    
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='align_inser 1.0')
    
    kwargs = {'basename':args['<basename>']}
    print('[superCasQ] basename: ' + str(kwargs['basename']))

    align_inser_to_vector(**kwargs)
    generate_vector_align_tab(**kwargs)
    align_vector_to_genome(**kwargs)

    print("\nalign_inser.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
