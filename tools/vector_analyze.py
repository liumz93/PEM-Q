#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.23
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
    vector_analyze   <basename> <vector_type>

Options:
-h --help               Show this screen.
-v --version            Show version.
<vector_type>           pX330_asCpf1; pX330_spCas9; lentiVirus

This script help to dedup reads according its barcode attatched, using levenshstein 
distance method. A bam file of adapter alignment is needed for this script. Both adapter 
and barcode are allowed 2 mismatches.

Input file: adapter alignment bam file / Output file: unique reads's query name list

Author: Mengzhu LIU
Last Update:2019.10.23

"""

import os
import pysam
from time import time
from docopt import docopt
import pandas as pd

def load(inputfile):
    
    if not os.path.exists(inputfile):
        raise ValueError('[PEM-Q] The {} file does not exist.'.format(inputfile))

def align_to_vector(basename,vector_type):
    
    os.system("mkdir vector")
    #~~~~~~~~align ped_fq to vector~~~~~~~~~~#
    
    vector_fa = "/home/mengzhu/database/{}/vector.fa".format(vector_type)
    pe_fq_r1 = basename + "_R1.fq.gz"
    pe_fq_r2 = basename + "_R2.fq.gz"
    
    load(vector_fa)
    load(pe_fq_r1)
    load(pe_fq_r1)
    
    pe_sam = basename + '_pe_vector.sam'
    pe_bam = basename + '_pe_vector.bam'
    pe_bam_sort = basename + '_pe_vector.sort.bam'
    
    #build bwa index of vector
    cmd = "samtools faidx /home/mengzhu/database/{}/vector.fa".format(vector_type)
    os.system(cmd)
    cmd = "bwa index -a bwtsw -p /home/mengzhu/database/{}/vector /home/mengzhu/database/{}/vector.fa 1>vector/build_index.o 2>vector/build_index.e".format(vector_type,vector_type)
    os.system(cmd)

    #alignment
    print("[PEM-Q]  align pe_fq to vector...")
                  
    cmd = "bwa mem -t 8 /home/mengzhu/database/{}/vector {} {} > vector/{} 2>vector/bwa_align_pe_vector.log".format(vector_type, pe_fq_r1, pe_fq_r2, pe_sam)
    print(cmd )
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(pe_sam, pe_bam, pe_bam, pe_bam_sort , pe_bam_sort)
    print("[PEM-Q]  sort and index bam...")
    print(cmd)
    os.system(cmd) 

def primer_filter(basename,vector_type):
    
    print("[PEM-Q]  processing primer filter...")
    pe_bam_sort = basename + '_pe_vector.sort.bam'
    pe_primer_bam = basename + '_primer_vector.bam'
    pe_primer_bam_sort = basename + '_primer_vector.sort.bam'
    
    primer_list_file = pd.read_csv("primer/bamlist_stitch.txt",sep = ' ',names = ["Qname", "Bait_start", "Bait_end"])
    primer_list = primer_list_file["Qname"]
    vector_pe_bam = pysam.AlignmentFile("vector/"+pe_bam_sort,'rb')
    vector_primer_bam = pysam.AlignmentFile("vector/"+pe_primer_bam, "wb", template=vector_pe_bam)
    
    vector_pe_bam_indexed = pysam.IndexedReads(vector_pe_bam)
    vector_pe_bam_indexed.build()
    n = 0
    for name in primer_list:
        try:
            vector_pe_bam_indexed.find(name)
        except KeyError:
            pass
        else:
            
            iterator = vector_pe_bam_indexed.find(name)
            for x in iterator:
                n = n + 1
                vector_primer_bam.write(x)
    print("primer filter left:",n)
    vector_pe_bam.close()
    vector_primer_bam.close()
    pysam.sort("-o", "vector/"+pe_primer_bam_sort, "vector/"+pe_primer_bam)
    cmd = "samtools index vector/{}".format(pe_primer_bam_sort)
    os.system(cmd)
    
def proper_pair_tab(basename,vector_type):
    
    print("[PEM-Q]  generating proper pair tab...")
    
    total_vector_len = 8484
    
    pe_primer_bam_sort = "vector/"+ basename + '_primer_vector.sort.bam'
    vector_pe_bam = pysam.AlignmentFile(pe_primer_bam_sort,'rb')
    pairedreads = pysam.AlignmentFile("vector/"+basename+"_pe_vector.paired.bam", "wb", template=vector_pe_bam)
    r1 = pysam.AlignmentFile("vector/"+basename+"_r1.paired.bam", "wb", template=vector_pe_bam)
    r2 = pysam.AlignmentFile("vector/"+basename+"_r2.paired.bam", "wb", template=vector_pe_bam)
    
    # get paired and both mapped reads
    n = 0
    m = 0
    k = 0
    for read in vector_pe_bam:
        if read.is_paired and (not read.is_unmapped) and (not read.is_supplementary):
            pairedreads.write(read)
            n = n + 1
            if read.is_read1:
                m = m + 1
                read_type = "read1"
                r1.write(read)
            else:
                k = k + 1
                read_type = "read2"
                r2.write(read)
    pairedreads.close()
    print("paired:",n,"r1:",n,"r2:",n)
    r1.close()
    r2.close()
    pysam.sort("-o", "vector/"+basename+"_pe_vector.paired.sort.bam", "vector/"+basename+"_pe_vector.paired.bam")
    pysam.sort("-o", "vector/"+basename+"_r1.paired.sort.bam", "vector/"+basename+"_r1.paired.bam")
    pysam.sort("-o", "vector/"+basename+"_r2.paired.sort.bam", "vector/"+basename+"_r2.paired.bam")
    
    #extract r1,r2 end of vector
    r1_bam = pysam.AlignmentFile("vector/"+basename+"_r1.paired.sort.bam",'rb')
    r2_bam = pysam.AlignmentFile("vector/"+basename+"_r2.paired.sort.bam",'rb')
    r1_name_indexed = pysam.IndexedReads(r1_bam)
    r2_name_indexed = pysam.IndexedReads(r2_bam)
    r1_name_indexed.build()
    r2_name_indexed.build()
    
    vector_tab = open("vector/"+basename+"_vector.tab", "w")
    vector_tab.write("Qname"+"\t"+
                        "Vector_start"+"\t"+
                        "Vector_end"+"\t"+
                        "Vector_strand"+"\t"+
                        "Align_sequence"+"\t"+
                        "Vector_inser_size"+"\n")
    for read in r1_bam:
        name = read.query_name
        
        #find end from read2
        try:
            r2_name_indexed.find(name)
        except KeyError:
            um = "unmapped"
            pass
        else:
            iterator = r2_name_indexed.find(name)
            for x in iterator:
                read2 = x
                um = "mapped"
                
        if um == "unmapped":
            vector_start = read.reference_start + 1
            vector_end = read.reference_end
            inser_len = vector_end-vector_start + 1
            if read.is_reverse:
                strand = "-"
            else:
                strand = "+"
        else:
            if read.is_reverse and (not read2.is_reverse):
                strand = "-"
                vector_start = read2.reference_start + 1
                vector_end = read.reference_end
                check2 = read.reference_start - read2.reference_start
            elif read2.is_reverse and (not read.is_reverse):
                strand = "+"
                vector_start = read.reference_start + 1
                vector_end = read2.reference_end
                check2 = read2.reference_start - read.reference_start
            else:
                continue
            
            check = vector_end-vector_start + 1
            if check < 0:
                inser_len = total_vector_len - abs(check)
            else:
                if check2 < 0:
                    inser_len = 0
                else:
                    inser_len = check
        qname = read.query_name
        vector_start = vector_start
        vector_end = vector_end
        align_seq = read.query_alignment_sequence
        # aslign_seq2 = read2.query_alignment_sequence
        inser_len = inser_len
        vector_tab.write(qname+"\t"+
                            str(vector_start)+"\t"+
                            str(vector_end)+"\t"+
                            str(strand)+"\t"+
                            align_seq+"\t"+
                            str(inser_len)+"\n")
    vector_tab.close()    
    # merge vector tab with primer info and rmb info
    vector_tab_file = pd.read_csv("vector/" + basename + "_vector.tab", sep = '\t')
    primer_list_file = pd.read_csv("primer/bamlist_stitch.txt",sep = ' ', names = ["Qname", "Bait_start", "Bait_end"])
    primer_list_file.to_csv("vector/bamlist_stitch.txt",sep = '\t', index=False, header = True)
    primer_list_file = pd.read_csv("vector/bamlist_stitch.txt",sep = '\t', names = ["Qname", "Bait_start", "Bait_end"],low_memory=False)
    rmb = pd.read_csv("barcode/" + basename + "_barcode_list.txt",sep = '\t', names = ["Qname", "Barcode"])
    vector_merge1_tab = pd.merge(vector_tab_file, primer_list_file, on='Qname', how='inner')
    vector_merge2_tab = pd.merge(vector_merge1_tab, rmb, on='Qname', how='inner')
    vector_merge2_tab.to_csv("vector/" + basename + "vector.tab", header = True, sep = '\t', index=False)
    
    print("[PEM-Q]  rmb dedup...")
    os.system("mkdir unique")
    # rmb dedup
    vector_dedup_tab = vector_merge2_tab.drop_duplicates(['Bait_start','Bait_end','Vector_start','Vector_end','Barcode'],keep='first')
    vector_dedup_tab.to_csv("unique/" + basename + "_vector_uni.tab", header = True, sep = '\t', index=False)
        
def remove_those_align_to_genome(basename,vector_type):
    
    #generate fa file to align to genome
    vector_file = pd.read_csv("unique/" + basename + "_vector_uni.tab",sep = '\t')
    vector_fa = open("vector/" + basename + "_vector.fa","w")
    n = 0
    for i in range(0,len(vector_file['Align_sequence'])):
        vector_fa.write(">"+vector_file['Qname'][i]+"\n")
        vector_fa.write(vector_file['Align_sequence'][i]+"\n")
        n += 1
    vector_fa.close()
    print("Vector fa:",n)
    
    #alignment
    print("[PEM-Q]  align vector fa to genome...")
    bwa_index_path = "/home/mengzhu/database/bwa_indexes/{}/{}".format("hg38", "hg38")
    cmd = "bwa mem -t 2 {} {} > {} 2>vector/bwa_align_vector_fa.log".format(bwa_index_path,"vector/" + basename + "_vector.fa", "vector/" + basename + "_vector_fa.sam")
    print(cmd)
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(basename + "_vector_fa.sam", 
                                               basename + "_vector_fa.bam", 
                                               basename + "_vector_fa.bam", 
                                               basename + "_vector_fa.sort.bam", 
                                               basename + "_vector_fa.sort.bam")
    print("[PEM-Q]  sort and index bam...")
    print(cmd)
    os.system(cmd)
    
    #remove reads align to genome
    print("[PEM-Q]  remove vector sequence align to genome...")
    vector_fa_bam = pysam.AlignmentFile("vector/"+basename + "_vector_fa.sort.bam",'rb')
    vector_genome_file = open("unique/" + basename + "_vector_genome_uni.tab","w")
    vector_genome_file.write("Qname" + "\t" +
                            "In_Genome" + "\t" +
                            "Gname" + "\t" +
                            "Gstart" + "\t" +
                            "Gend" + "\n")
    n = 0
    m = 0
    for read in vector_fa_bam:
        m += 1
        if read.reference_name is None:
            qname = read.query_name
            gname = ""
            genome_flag = "N"
            gstart = ""
            gend = ""
        else:
            qname = read.query_name
            genome_flag = "Y"
            gname = read.reference_name
            gstart = read.reference_start + 1
            gend = read.reference_end
            n += 1
        vector_genome_file.write(qname + "\t" +
                                genome_flag + "\t" +
                                gname + "\t" +
                                str(gstart) + "\t" +
                                str(gend) + "\n")
    print("Total vector:", m, "Align to genome:", n)
    vector_genome_file.close()
    vector_genome_file = pd.read_csv("unique/" + basename + "_vector_genome_uni.tab",sep = '\t')
    # merge vector_genome_file and vector_file
    vector_merge_file = pd.merge(vector_file, vector_genome_file, on='Qname', how='inner')
    vector_merge_dedup = vector_merge_file.drop_duplicates(['Qname'],keep='first')
    vector_merge_dedup.to_csv("unique/" + basename + "_vector.tab", header = True, sep = '\t', index=False)
    #generate no genome vector file
    vector_no_genome_file = vector_merge_dedup[vector_merge_file["In_Genome"] == "N"]
    vector_no_genome_file.to_csv("unique/" + basename + "_vector_NoGenome.tab", header = True, sep = '\t', index=False)
    #generate genome vector file
    vector_genome_file = vector_merge_dedup[vector_merge_file["In_Genome"] == "Y"]
    vector_genome_file.to_csv("unique/" + basename + "_vector_Genome.tab", header = True, sep = '\t', index=False)
    #clean file
    os.system("rm " + "unique/" + basename + "_vector_genome_uni.tab " + "unique/" + basename + "_vector_uni.tab")
def main():
    start_time = time()
    
    args = docopt(__doc__,version='vector_analyze 1.0')
    
    kwargs = {'basename':args['<basename>'],'vector_type':args['<vector_type>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] vector_type: ' + str(kwargs['vector_type']))
    
    ## function ##
    align_to_vector(**kwargs)
    primer_filter(**kwargs)
    proper_pair_tab(**kwargs)
    remove_those_align_to_genome(**kwargs)

    print("\nvector_analyze.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()