#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.4.23
"""define_transloc

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
    define_transloc   <basename> <cutsite>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script find bait break site and translocation site from bwa align file.
And also filter reads for indel define.

Input file: reads alignment bam file / Output file: a informative tab file

Author: Mengzhu LIU
Last Update:2019.5.8

"""

    ######################################
    ## i)software dependencies:
    ##    None
    ## ii)data needed:
    ##    1.adapter alignment bam file
    ## iii)python package:
    ##    1.os
    ##    2.pysam
    ##    3.time
    ##    4.docopt
    ##    5.pandas
    ##    6.re
    #######################################
    
import os
import pysam
import re
from time import time
from docopt import docopt
import pandas as pd
from Bio import SearchIO
from Bio.Seq import Seq

class Define_transloc(object):
    
    def __init__(self, basename=None,cutsite=None):
        
        self.basename = basename
        bam_sort = "barcode/" + basename + "_sti.filter.sort.bam"
        self.bam_sort = bam_sort
        self.cutsite = int(cutsite)
    
    def load(self, f_file):
        
        if not os.path.exists(f_file):
            raise ValueError('[PEM-Q] The {} file does not exist.'.format(f_file))   
        
    def reverse_cigar_value(self,cigar):

        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)

        reverse_letter = letter[::-1]
        reverse_number = number[::-1]


        merge_cigar = [None]*(len(reverse_number)+len(reverse_letter))
        merge_cigar[::2] = reverse_number
        merge_cigar[1::2] = reverse_letter

        reverse_cigar = ''.join(merge_cigar)

        return(reverse_cigar)

    def cigar_map_len(self,map_type,cigar):

        #map_type contain: seq/referece
        #reference map len = total - insertion
        #sequence map len = total - deletion
        
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int
        
        #here only want to know mapping length
        #so don't consider strand info and reverse cigar 
        #map_len = (start M + .. + end M) - D
        
        match_list = []
        dele_list = []
        inser_list = []
        for i in range(0,len(letter)):
            if letter[i] == 'M':
                match_list.append(i)
            if letter[i] == 'D':
                dele_list.append(number[i])
            if letter[i] == 'I':
                inser_list.append(number[i])
                
        index_min = min(match_list)
        index_max = max(match_list)

        map_sum = sum(number[index_min:index_max+1])
        del_len = sum(dele_list)
        inser_len = sum(inser_list)
        
        if map_type == 'reference':
            map_len = map_sum - inser_len
        if map_type == 'sequence':
            map_len = map_sum - del_len
        
        return(map_len)

    def cigar_map_seq_start(self,cigar):

        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int

        for i in range(0,len(letter)):
            if letter[i] == 'M':
                map_seq_start = sum(number[0:i]) + 1
                break
        return(map_seq_start)

    def cigar_map_seq_end(self,cigar):

        map_seq_start = self.cigar_map_seq_start(cigar)
        map_seq_len = self.cigar_map_len('sequence',cigar)
        map_seq_end = map_seq_start + map_seq_len - 1

        return(map_seq_end)

    def transloc_microhomo(self, rep_sequence, rep_strand, rep_cigarstring, consistant_sup_cigarstring):

        # be careful when use this function because
        # rep_cigarstring have consistant sequential order with sup_cigarstring when have the same trand
        # and vice-versa

        rep_map_seq_start = self.cigar_map_seq_start(rep_cigarstring)
        rep_map_seq_end = self.cigar_map_seq_end(rep_cigarstring)

        sup_map_seq_start = self.cigar_map_seq_start(consistant_sup_cigarstring)
        sup_map_seq_end = self.cigar_map_seq_end(consistant_sup_cigarstring)

        microhomo = ''
        if rep_strand == '-':
            start = rep_map_seq_start
            end = sup_map_seq_end
        else:
            start = sup_map_seq_start
            end = rep_map_seq_end

        if end >= start:
            microhomo = rep_sequence[start-1 : end]
        
        return(microhomo)


    def transloc_find_insertion(self, rep_sequence, rep_strand, rep_cigarstring, consistant_sup_cigarstring):

        rep_map_seq_start = self.cigar_map_seq_start(rep_cigarstring)
        rep_map_seq_end = self.cigar_map_seq_end(rep_cigarstring)

        sup_map_seq_start = self.cigar_map_seq_start(consistant_sup_cigarstring)
        sup_map_seq_end = self.cigar_map_seq_end(consistant_sup_cigarstring)

        insertion = ''
        if rep_strand == '-':
            start = rep_map_seq_start
            end = sup_map_seq_end
        else:
            start = sup_map_seq_start
            end = rep_map_seq_end

        if (start - end) > 1:
            insertion = rep_sequence[end : start-1]
        
        return(insertion)

    def generate_transloc_tab(self):

        os.system("mkdir transloc")
        os.system("mkdir indel")
        
        bam_file = pysam.AlignmentFile(self.bam_sort, 'rb')
        transloc_tab = open("transloc/" + self.basename + "_mut.tab", "w")
        discard_tab = open("transloc/" + self.basename + "_discard.tab", "w")
        indel_bam = pysam.AlignmentFile("indel/" + self.basename + "_indel.bam", "wb", template=bam_file)

        transloc_tab.write('Qname'+"\t"+\
                            'Bait_rname'+"\t"+\
                            'Bait_strand'+"\t"+\
                            'Bait_start'+"\t"+\
                            'Bait_end'+"\t"+\
                            'Prey_rname'+"\t"+\
                            'Prey_strand'+"\t"+\
                            'Prey_start'+"\t"+\
                            'Prey_end'+"\t"+\
                            'Rname'+"\t"+\
                            'Strand'+"\t"+\
                            'Junction'+"\t"+\
                            'Sequence'+"\t"+\
                            'B_Qstart'+"\t"+\
                            'B_Qend'+"\t"+\
                            'Qstart'+"\t"+\
                            'Qend'+"\t"+\
                            'Qlen'+"\t"+\
                            'Insertion'+"\t"+\
                            'Microhomolog'+"\t"+\
                            'Prey_MQ'+"\n")
        c = 0
        for read in bam_file:

            # because in align_make, reads are already filtered by primer
            # here we only need to classify translocation reads/indel reads/germline reads
            
            #define strand
            
            reference_start = read.reference_start + 1
            reference_end = read.reference_end
            if read.is_reverse:
                read_strand = '-'
            else:
                read_strand = '+'
            
            primer_length = 20   
            #reads contain 'SA' tag have chimeric alignment are potential translocation reads
            
            transloc_check = abs(reference_start - read.reference_end)

            condition1 = any('SA' == tg[0] for tg in read.get_tags())
            condition2 = transloc_check > primer_length + 10
            # condition3 = transloc_check > self.cutsite + 10
            # if condition3:
            #     indel_bam.write(read)
            # else:
            if condition2:
                if condition1:
                    sa_tag = read.get_tag('SA')
                    # multiple transloclation filter (2020.12.29)
                    sa_number = len(sa_tag.split(';'))

                    if sa_number > 2:
                        # print(sa_tag)
                        c = c + 1
                        continue
                    
                    sa_tag_list = sa_tag.split(',')
                    read_query_name = read.query_name
                    rep_rname = read.reference_name
                    rep_reference_start = reference_start
                    rep_reference_end = read.reference_end
                    rep_strand = read_strand
                    rep_cigarstring = read.cigarstring

                    sup_cigarstring = sa_tag_list[3]
                    sup_map_ref_len = self.cigar_map_len('reference',sup_cigarstring)
            
                    sup_rname = sa_tag_list[0]
                    sup_reference_start = sa_tag_list[1]
                    sup_strand = sa_tag_list[2]
                    sup_reference_end = int(sup_reference_start) + sup_map_ref_len - 1
                    
                    #define transloc junction
                    if sup_strand == '+':
                        transloc_junction = sup_reference_start
                    else:
                        transloc_junction = sup_reference_end

                    #define microhomo sequence

                    # rep_cigarstring have consistant sequential order with sup_cigarstring when have the same trand
                    # and vice-versa
            
                    consistant_sup_cigarstring = sup_cigarstring
                    if sup_strand != rep_strand:
                        consistant_sup_cigarstring = self.reverse_cigar_value(sup_cigarstring)#reverse cigarstring

                    rep_sequence = read.query_sequence
                    microhomo = self.transloc_microhomo(rep_sequence,rep_strand,rep_cigarstring,consistant_sup_cigarstring)

                    # define insertion sequence
                    insertion = self.transloc_find_insertion(rep_sequence,rep_strand,rep_cigarstring,consistant_sup_cigarstring)
                    
                    #sequence info
                    if read.is_reverse:
                        read_strand = '-'
                        sequence = Seq(read.query_sequence)
                        sequence = str(sequence.complement())
                    else:
                        read_strand = '+'
                        sequence = read.query_sequence
                        
                    B_Qstart = self.cigar_map_seq_start(rep_cigarstring)
                    B_Qend = self.cigar_map_seq_end(rep_cigarstring)
                    Qstart = self.cigar_map_seq_start(consistant_sup_cigarstring)
                    Qend = self.cigar_map_seq_end(consistant_sup_cigarstring)
                    Qlen = len(rep_sequence)
                    sup_mapqual = sa_tag_list[4]
                        
                    transloc_tab.write(read_query_name+"\t"+\
                                       rep_rname+"\t"+\
                                       rep_strand+"\t"+\
                                       str(rep_reference_start)+"\t"+\
                                       str(rep_reference_end)+"\t"+\
                                       sup_rname+"\t"+\
                                       sup_strand+"\t"+\
                                       str(sup_reference_start)+"\t"+\
                                       str(sup_reference_end)+"\t"+\
                                       sup_rname+"\t"+\
                                       sup_strand+"\t"+\
                                       str(transloc_junction)+"\t"+\
                                       sequence+"\t"+\
                                       str(B_Qstart)+"\t"+\
                                       str(B_Qend)+"\t"+\
                                       str(Qstart)+"\t"+\
                                       str(Qend)+"\t"+\
                                       str(Qlen)+"\t"+\
                                       insertion+"\t"+\
                                       microhomo+"\t"+\
                                       str(sup_mapqual)+"\n")
                    #reads do not contain 'SA' tag are potential indel/germline reads
                else:
                    indel_bam.write(read)
            else:
                discard_tab.write(read.query_name+"\n")
        print('multiple junctions remove',c)    
        bam_file.close()
        transloc_tab.close()
        indel_bam.close()
        discard_tab.close()
        pysam.sort("-o", "indel/" + self.basename + "_indel.sort.bam", "indel/" + self.basename + "_indel.bam")
        
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='denfine_transloc 1.0')
    
    kwargs = {'basename':args['<basename>'], 'cutsite':args['<cutsite>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    
    define_transloc = Define_transloc(**kwargs)
    define_transloc.generate_transloc_tab()
    

    print("\ndefine_transloc.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
