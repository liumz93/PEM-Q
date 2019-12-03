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
    find_mh_of_vector   <basename> <primer_strand>

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
import numpy as np

def alignment_extraction(basename, primer_strand):
    
    vector_all_file_name = basename + "_vector_suspected_inser.tab"
    vector_all_file = pd.read_csv(vector_all_file_name, sep = '\t')
    print(vector_all_file['Qname'].count())
    condition = (vector_all_file['Vector_start'] >= 7925) & \
                (vector_all_file['Vector_end'] <= 7965)
    vector_all_file = vector_all_file[~condition]
    print(vector_all_file['Qname'].count())
    type_of_vector = "Confident"
    condition = vector_all_file['Type'] == type_of_vector
    # vector_all_file = vector_all_file[condition]
    qname_list = vector_all_file['Qname'][1:]
    qname_list.to_csv(type_of_vector + "_qname.list", header = False,  index=False)
    
    #### bait ####
    
    #extract bait bam
    bait_bam_file_name = "primer/" + basename + "_sti.p.sort.bam"
    qname_list_name = type_of_vector + "_qname.list"
    name = type_of_vector + "_bait"
    cmd = "extract_reads_of_qname.py {} {} {}".format(bait_bam_file_name, qname_list_name, name)
    print(cmd)
    os.system(cmd)
    extract_bait_bam_file_name = name + "_extrac.sort.bam"
    
    #catch sequence pos info
    bait_sequence_info = open("vector/"+basename+"_bait_info.tab","w")
    bait_sequence_info.write("Qname" + "\t" +
                                "Bait_Cigar" + "\t" +
                                "Bait_start_sp" + "\t" +
                                "Bait_end_sp" + "\n")
    extract_bait_bam_file = pysam.AlignmentFile(extract_bait_bam_file_name,'rb')
    for read in extract_bait_bam_file:
        qname = read.query_name
        # print(qname)
        bait_cigar = read.cigarstring
        qstart = read.query_alignment_start + 1
        qend = read.query_alignment_end
        bait_sequence_info.write(qname + "\t" +
                                    bait_cigar + "\t" +
                                    str(qstart) + "\t" +
                                    str(qend) + "\n")
    bait_sequence_info.close()
    
    #### prey ####
    
    #extract prey bam
    prey_bam_file_name = "vector/" + basename + "_r2_genome.sort.bam"
    qname_list_name = type_of_vector + "_qname.list"

    #Realign R2 to genome
    outfq_r2 = basename+"_vector_R2.fq"
    cmd = "seqtk subseq {} {} > {}".format(basename+"_R2.fq.gz", qname_list_name, outfq_r2)
    os.system(cmd)
   
    r2_sam = basename + '_vector_r2_genome.sam'
    r2_bam = basename + '_vector_r2_genome.bam'
    r2_bam_sort = basename + '_vector_r2_genome.sort.bam'
    
    print("[PEM-Q Vector Analysis]  align r2 to genome...")
                  
    cmd = "bwa mem -t 8 -k 10 /home/mengzhu/database/hg38/hg38 {} > vector/{} 2>vector/bwa_align_pe_vector.log".format(outfq_r2, r2_sam)
    print(cmd )
    os.system(cmd)
    
    cmd = "samtools view -S -b -h vector/{} > vector/{} \
           && samtools sort vector/{} > vector/{} \
           && samtools index vector/{}".format(r2_sam, r2_bam, r2_bam, r2_bam_sort , r2_bam_sort)
    print("[PEM-Q Vector Analysis]  sort and index bam...")
    print(cmd)
    os.system(cmd)
    
    
    name = type_of_vector + "_prey"
    cmd = "extract_reads_of_qname.py vector/{} {} {}".format(r2_bam_sort, qname_list_name, name)
    print(cmd)
    os.system(cmd)
    extract_prey_bam_file_name = name + "_extrac.sort.bam"
    
    #catch sequence pos info
    prey_sequence_info = open("vector/"+basename+"_prey_info.tab","w")
    prey_sequence_info.write("Qname" + "\t" +
                                "Prey_Cigar" + "\t" +
                                "Prey_start_sp" + "\t" +
                                "Prey_end_sp" + "\n")
    extract_prey_bam_file = pysam.AlignmentFile(extract_prey_bam_file_name,'rb')
    for read in extract_prey_bam_file:
        qname = read.query_name
        # print(qname)
        if read.reference_name is None:
            prey_cigar = ""
            qstart = ""
            qend = ""
        else:
            prey_cigar = read.cigarstring
            qstart = read.query_alignment_start + 1
            qend = read.query_alignment_end
        prey_sequence_info.write(qname + "\t" +
                                    prey_cigar + "\t" +
                                    str(qstart) + "\t" +
                                    str(qend) + "\n")
    prey_sequence_info.close()
    
    #### Vector R1 ####
    
    #extract vector R1
    vector_r1_bam_file_name = "vector/" + basename + "_r1.paired.sort.bam"
    qname_list_name = type_of_vector + "_qname.list"
    name = type_of_vector + "_vector_r1"
    cmd = "extract_reads_of_qname.py {} {} {}".format(vector_r1_bam_file_name, qname_list_name, name)
    print(cmd)
    os.system(cmd)
    extract_vector_r1_bam_file_name = name + "_extrac.sort.bam"
    
    #catch sequence pos info
    vector_r1_sequence_info = open("vector/"+basename+"_vector_r1_info.tab","w")
    vector_r1_sequence_info.write("Qname" + "\t" +
                                "Vector_R1_Cigar" + "\t" +
                                "Vector_R1_Strand" + "\t" +
                                "Vector_R1_start_sp" + "\t" +
                                "Vector_R1_end_sp" + "\n")
    extract_vector_r1_bam_file = pysam.AlignmentFile(extract_vector_r1_bam_file_name,'rb')
    for read in extract_vector_r1_bam_file:
        qname = read.query_name
        # print(read.is_reverse)
        if read.is_reverse:
            vector_r1_strand = "-"
        else:
            vector_r1_strand = "+"
        # print(qname)
        vector_r1_cigar = read.cigarstring
        # print(vector_r1_cigar)
        qstart = read.query_alignment_start + 1
        qend = read.query_alignment_end
        vector_r1_sequence_info.write(qname + "\t" +
                                    vector_r1_cigar + "\t" +
                                    vector_r1_strand + "\t" +
                                    str(qstart) + "\t" +
                                    str(qend) + "\n")
    vector_r1_sequence_info.close()
    
    #### Vector R2 ####
    
    #extract vector R2
    vector_r2_bam_file_name = "vector/" + basename + "_r2.paired.sort.bam"
    qname_list_name = type_of_vector + "_qname.list"
    name = type_of_vector + "_vector_r2"
    cmd = "extract_reads_of_qname.py {} {} {}".format(vector_r2_bam_file_name, qname_list_name, name)
    print(cmd)
    os.system(cmd)
    extract_vector_r2_bam_file_name = name + "_extrac.sort.bam"
    
    #catch sequence pos info
    vector_r2_sequence_info = open("vector/"+basename+"_vector_r2_info.tab","w")
    vector_r2_sequence_info.write("Qname" + "\t" +
                                "Vector_R2_Cigar" + "\t" +
                                "Vector_R2_Strand" + "\t" +
                                "Vector_R2_start_sp" + "\t" +
                                "Vector_R2_end_sp" + "\n")
    extract_vector_r2_bam_file = pysam.AlignmentFile(extract_vector_r2_bam_file_name,'rb')
    for read in extract_vector_r2_bam_file:
        qname = read.query_name
        # vector_r2_sequence= read.query_sequence
        if read.is_reverse:
            vector_r2_strand = "-"
        else:
            vector_r2_strand = "+"
        # print(qname)
        vector_r2_cigar = read.cigarstring
        qstart = read.query_alignment_start + 1
        qend = read.query_alignment_end
        vector_r2_sequence_info.write(qname + "\t" +
                                    vector_r2_cigar + "\t" +
                                    vector_r2_strand + "\t" +
                                    str(qstart) + "\t" +
                                    str(qend) + "\n")
    vector_r2_sequence_info.close()
    
def find_mh(basename, primer_strand):
    
    bait_sequence_info = "vector/"+basename+"_bait_info.tab"
    prey_sequence_info = "vector/"+basename+"_prey_info.tab"
    vector_r1_sequence_info = "vector/"+basename+"_vector_r1_info.tab"
    vector_r2_sequence_info = "vector/"+basename+"_vector_r2_info.tab"
    
    bait_sequence_info_file = pd.read_csv(bait_sequence_info, sep = '\t')
    prey_sequence_info_file = pd.read_csv(prey_sequence_info, sep = '\t')
    vector_r1_sequence_info_file = pd.read_csv(vector_r1_sequence_info, sep = '\t')
    vector_r2_sequence_info_file = pd.read_csv(vector_r2_sequence_info, sep = '\t')
    # print(vector_r1_sequence_info_file['Vector_R1_Strand'].head())
    vector_mh_file = pd.merge(bait_sequence_info_file, vector_r1_sequence_info_file, on='Qname', how='inner')
    # print(vector_mh_file['Vector_R1_Strand'].head())
    vector_mh_file2 = pd.merge(vector_mh_file, prey_sequence_info_file, on='Qname', how='inner')
    vector_mh_file3 = pd.merge(vector_mh_file2, vector_r2_sequence_info_file, on='Qname', how='inner')
    
    vector_mh_file3.loc[:,'VS']= np.where(vector_mh_file3['Vector_R1_Strand'] == primer_strand, vector_mh_file3.Vector_R1_start_sp, 150 - vector_mh_file3.Vector_R1_start_sp)
    vector_mh_file3.loc[:,'VE']= np.where(vector_mh_file3['Vector_R1_Strand'] == primer_strand, vector_mh_file3.Vector_R1_end_sp, 150 - vector_mh_file3.Vector_R1_end_sp)
    vector_mh_file3['Bait_MH'] = vector_mh_file3['Bait_end_sp'] - vector_mh_file3['VS'] + 1
    
    
    vector_mh_file3.to_csv(basename + "_vector_mh.tab", header = True, sep = '\t', index=False, columns = [u'Qname',u'Bait_MH',
    u'Bait_start_sp', u'Bait_end_sp', u'Bait_Cigar',u'Vector_R1_start_sp',u'Vector_R1_end_sp', u'VS', u'VE', u'Vector_R1_Cigar', u'Vector_R1_Strand', u'Prey_start_sp',u'Prey_end_sp',
    u'Prey_Cigar', u'Vector_R2_start_sp',u'Vector_R2_end_sp', u'Vector_R2_Cigar'])
    
def main():
    start_time = time()
    
    args = docopt(__doc__,version='find_mh_of_vector 1.0')
    
    kwargs = {'basename':args['<basename>'],'primer_strand':args['<primer_strand>']}
    print('[PEM-Q Vector Analysis] basename: ' + str(kwargs['basename']))
    print('[PEM-Q Vector Analysis] primer_strand: ' + str(kwargs['primer_strand']))
    
    ## function ##
    alignment_extraction(**kwargs)
    find_mh(**kwargs)

    print("\nfind_mh_of_vector.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()