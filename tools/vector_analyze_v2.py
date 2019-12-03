#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.23
"""vector_analyze_v2

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
    vector_analyze   <basename> <vector_fa>

Options:
-h --help               Show this screen.
-v --version            Show version.
<vector_type>           pX330_asCpf1; pX330_spCas9; lentiVirus

This script .

Input file:  / Output file: 
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
        raise ValueError('[PEM-Q Vector Analysis] The {} file does not exist.'.format(inputfile))
    
    
def align_discard_to_vector(basename,vector_fa):
    
    #extract_reads_from_discard
    discard_list = "indel/" + basename + "_discard.tab"
    # discard_list = "unique/u6_list"
    outfq_r1 = basename+"_discard_R1.fq"
    outfq_r2 = basename+"_discard_R2.fq"
    # outfq_r1 = basename+"_u6_R1.fq"
    # outfq_r2 = basename+"_u6_R2.fq"
    cmd = "seqtk subseq {} {} > {}".format(basename+"_R1.fq.gz", discard_list, outfq_r1)
    print("[PEM-Q Vector Analysis]" + cmd)
    os.system(cmd)
    cmd = "seqtk subseq {} {} > {}".format(basename+"_R2.fq.gz", discard_list, outfq_r2)
    print("[PEM-Q Vector Analysis]" + cmd)
    os.system(cmd)
    
    os.system("mkdir vector")
    #~~~~~~~~align ped_fq to vector~~~~~~~~~~#
    print("[PEM-Q Vector Analysis]  check file...")
    # vector_fa = "/home/mengzhu/database/{}/vector.fa".format(vector_type)
    pe_fq_r1 = outfq_r1
    pe_fq_r2 = outfq_r2
    
    load(vector_fa)
    load(pe_fq_r1)
    load(pe_fq_r1)
    
    os.system("cp {} vector/".format(vector_fa))
    pe_sam = basename + '_pe_vector.sam'
    pe_bam = basename + '_pe_vector.bam'
    pe_bam_sort = basename + '_pe_vector.sort.bam'
    
    print("[PEM-Q Vector Analysis]  building vector index...")
    #build bwa index of vector
    cmd = "samtools faidx vector/{}".format(vector_fa)
    os.system(cmd)
    vector_fa_index = vector_fa.split(".")[0]
    cmd = "bwa index -a bwtsw -p vector/{} vector/{} 1>vector/build_index.o 2>vector/build_index.e".format(vector_fa_index,vector_fa)
    os.system(cmd)

    #alignment
    print("[PEM-Q Vector Analysis]  align pe_fq to vector...")
                  
    cmd = "bwa mem -t 8 vector/{} {} {} > vector/{} 2>vector/bwa_align_pe_vector.log".format(vector_fa_index, pe_fq_r1, pe_fq_r2, pe_sam)
    print(cmd )
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(pe_sam, pe_bam, pe_bam, pe_bam_sort , pe_bam_sort)
    print("[PEM-Q Vector Analysis]  sort and index bam...")
    print(cmd)
    os.system(cmd) 
    
    #align r2 to genome
    r2_sam = basename + '_r2_genome.sam'
    r2_bam = basename + '_r2_genome.bam'
    r2_bam_sort = basename + '_r2_genome.sort.bam'
    
    print("[PEM-Q Vector Analysis]  align r2 to genome...")
                  
    cmd = "bwa mem -t 8 -k 10 /home/mengzhu/database/hg38/hg38 {} > vector/{} 2>vector/bwa_align_pe_vector.log".format(pe_fq_r2, r2_sam)
    print(cmd )
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(r2_sam, r2_bam, r2_bam, r2_bam_sort , r2_bam_sort)
    print("[PEM-Q Vector Analysis]  sort and index bam...")
    print(cmd)
    os.system(cmd)
    
def primer_filter(basename,vector_fa):
    
    print("[PEM-Q Vector Analysis]  processing primer filter...")
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
    
def proper_pair_tab(basename,vector_fa):
    
    print("[PEM-Q Vector Analysis]  generating proper pair tab...")
    vector_fa_index = vector_fa.split(".")[0]
    
    vector_file = pd.read_csv(vector_fa,names=["v1"])
    total_vector = vector_file["v1"][1]
    total_vector_len = len(total_vector)
    vector_genome = open("vector/"+vector_fa_index+".genome","w")
    vector_genome.write("vector"+"\t"+str(total_vector_len))
    vector_genome.close()
    
    
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
    print("paired:",n,"r1:",m,"r2:",k)
    r1.close()
    r2.close()
    pysam.sort("-o", "vector/"+basename+"_pe_vector.paired.sort.bam", "vector/"+basename+"_pe_vector.paired.bam")
    pysam.sort("-o", "vector/"+basename+"_r1.paired.sort.bam", "vector/"+basename+"_r1.paired.bam")
    pysam.sort("-o", "vector/"+basename+"_r2.paired.sort.bam", "vector/"+basename+"_r2.paired.bam")
    
    #extract r1,r2 end of vector
    r1_bam = pysam.AlignmentFile("vector/"+basename+"_r1.paired.sort.bam",'rb')
    r2_bam = pysam.AlignmentFile("vector/"+basename+"_r2.paired.sort.bam",'rb')
    r2_genome_bam = pysam.AlignmentFile("vector/"+basename+"_r2_genome.sort.bam",'rb')
    r1_name_indexed = pysam.IndexedReads(r1_bam)
    r2_name_indexed = pysam.IndexedReads(r2_bam)
    r2_genome_name_indexed = pysam.IndexedReads(r2_genome_bam)
    r1_name_indexed.build()
    r2_name_indexed.build()
    r2_genome_name_indexed.build()
    
    vector_tab = open("vector/"+basename+"_vector.tab", "w")
    vector_tab.write("Qname"+"\t"+
                        "Vector_start"+"\t"+
                        "Vector_end"+"\t"+
                        "Vector_strand"+"\t"+
                        "Align_sequence"+"\t"+
                        "Align_sequence_R2"+"\t"+
                        "Vector_inser_size"+"\t"+
                        "Prey_rname"+"\t"+
                        "Prey_start"+"\t"+
                        "Prey_end"+"\t"+
                        "Type"+"\n")
    for read in r1_bam:
        name = read.query_name
        
        #find vector end from read2
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
        #find genome end from read2
        try:
            r2_genome_name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator_g = r2_genome_name_indexed.find(name)
            for x in iterator_g:
                read2_genome = x
                if read2_genome.reference_name is not None:
                     r2_genome_flag = "mapped"
                else:
                    r2_genome_flag = "unmapped"
                break
                
        if r2_genome_flag == "mapped":
            Prey_rname = read2_genome.reference_name
            if read.is_reverse:
                strand = "-"
                if read2_genome.is_reverse:
                    Prey_start = read.reference_start
                    Prey_end = read.reference_end
                else:
                    Prey_start = read.reference_end
                    Prey_end = read.reference_start
            else:
                strand = "+"
                if read2_genome.is_reverse:
                    Prey_start = read.reference_end
                    Prey_end = read.reference_start
                else:
                    Prey_start = read.reference_start
                    Prey_end = read.reference_end
        else:
            Prey_rname = ""
            Prey_start = ""
            Prey_end = ""
            
        if um == "unmapped":
            # print(read2_genome)
            align_seq_R2 = "unmapped"
            pair_flag = "Medium"
            pair_flag = "Suspected"
            vector_start = read.reference_start + 1
            vector_end = read.reference_end
            inser_len = vector_end-vector_start + 1
            if read.is_reverse:
                strand = "-"
            else:
                strand = "+"
            if r2_genome_flag == "mapped":
                pair_flag = "Half"
            else:
                pair_flag = "Discard"
        else:
            align_seq_R2 = read2.query_alignment_sequence
            pair_flag = "Large"
            if r2_genome_flag == "mapped":
                pair_flag = "Confident"
            else:
                pair_flag = "Suspected"
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
                # inser_len = total_vector_len - abs(check)
                inser_len = 0
            else:
                if check2 < 0:
                    inser_len = 0
                else:
                    inser_len = check
        qname = read.query_name
        vector_start = vector_start
        vector_end = vector_end
        align_seq_R1 = read.query_alignment_sequence
        # aslign_seq2 = read2.query_alignment_sequence
        inser_len = inser_len
        vector_tab.write(qname+"\t"+
                            str(vector_start)+"\t"+
                            str(vector_end)+"\t"+
                            str(strand)+"\t"+
                            align_seq_R1+"\t"+
                            align_seq_R2+"\t"+
                            str(inser_len)+"\t"+
                            Prey_rname+"\t"+
                            str(Prey_start)+"\t"+
                            str(Prey_end)+"\t"+
                            pair_flag+"\n")
    vector_tab.close()
    r1_bam.close()
    r2_bam.close()
    r2_genome_bam.close()
    # merge vector tab with primer info and rmb info
    vector_tab_file = pd.read_csv("vector/" + basename + "_vector.tab", sep = '\t')
    primer_list_file = pd.read_csv("primer/bamlist_stitch.txt",sep = ' ', names = ["Qname", "Bait_start", "Bait_end"])
    primer_list_file.to_csv("vector/bamlist_stitch.txt",sep = '\t', index=False, header = True)
    primer_list_file = pd.read_csv("vector/bamlist_stitch.txt",sep = '\t', names = ["Qname", "Bait_start", "Bait_end"],low_memory=False)
    rmb = pd.read_csv("barcode/" + basename + "_barcode_list.txt",sep = '\t', names = ["Qname", "Barcode"])
    vector_merge1_tab = pd.merge(vector_tab_file, primer_list_file, on='Qname', how='inner')
    vector_merge2_tab = pd.merge(vector_merge1_tab, rmb, on='Qname', how='inner')
    vector_merge2_tab.to_csv("vector/" + basename + "_vector.tab", header = True, sep = '\t', index=False)
    
    print("[PEM-Q Vector Analysis]  rmb dedup...")
    os.system("mkdir unique")
    # rmb dedup
    vector_dedup_tab = vector_merge2_tab.drop_duplicates(['Bait_start','Bait_end','Vector_start','Vector_end','Barcode'],keep='first')
    vector_dedup_tab.to_csv(basename + "_vector_baitonly_inser.tab", header = True, sep = '\t', index=False, columns = [u'Qname',
    u'Vector_start', u'Vector_end', u'Vector_strand',u'Vector_inser_size', u'Bait_start', u'Bait_end', u'Prey_rname',
    u'Prey_start', u'Prey_end', u"Align_sequence", u"Align_sequence_R2",u"Type",u'Barcode'])
    # vector_dedup_tab.to_csv(basename + "_vector_u6_inser.tab", header = True, sep = '\t', index=False, columns = [u'Qname',
    # u'Vector_start', u'Vector_end', u'Vector_strand', u'Vector_inser_size', u'Bait_start', u'Bait_end', u'Prey_rname',u'Prey_start', u'Prey_end', u'Barcode'])
    
    os.system("align_inser_va.py {} -i {}".format(basename,vector_fa))
    vector_baitonly = pd.read_csv(basename + "_vector_baitonly_inser.tab", sep = '\t')
    vector_insertion = pd.read_csv(basename + "_vector_confident_inser.tab", sep = '\t')
    vector_final = vector_baitonly.append([vector_insertion])
    vector_final.to_csv(basename + "_all_vector.tab", header = True, sep = '\t', index=False, columns = [u'Qname',
    u'Vector_start', u'Vector_end', u'Vector_strand',u'Vector_inser_size', u'Bait_start', u'Bait_end', u'Prey_rname',
    u'Prey_start', u'Prey_end', u"Align_sequence", u"Align_sequence_R2",u"Type",u'Barcode'])
    
def main():
    start_time = time()
    
    args = docopt(__doc__,version='vector_analyze 1.0')
    
    kwargs = {'basename':args['<basename>'],'vector_fa':args['<vector_fa>']}
    print('[PEM-Q Vector Analysis] basename: ' + str(kwargs['basename']))
    print('[PEM-Q Vector Analysis] vector_fa: ' + str(kwargs['vector_fa']))
    
    ## function ##
    align_discard_to_vector(**kwargs)
    primer_filter(**kwargs)
    proper_pair_tab(**kwargs)

    print("\nvector_analyze.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()