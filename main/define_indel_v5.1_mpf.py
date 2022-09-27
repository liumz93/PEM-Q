#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.5.8
"""define_indel

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
    define_indel    <basename>  <cutsite> <primer_strand>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .

Input file: indel alignment bam file / Output file: a informative tab file

Author: Mengzhu LIU
Last Update:2019.6.10

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
import numpy as np

class Define_indel(object):
    
    def __init__(self, basename=None,cutsite=None,primer_strand=None):
        
        self.basename = basename
        bam_sort = "indel/" + basename + "_indel.sort.bam"
        self.bam_sort = bam_sort
        self.cutsite = int(cutsite)
        self.primer_strand = primer_strand
    
    def load(self, f_file):
        
        if not os.path.exists(f_file):
            raise ValueError('[PEM-Q] The {} file does not exist.'.format(f_file))   


    def cigar_find_delelen(self,cigar):
        
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int
        
        #here only want to know mapping length
        #so don't consider strand info and reverse cigar 
        
        dele_len = 0
        for i in range(0,len(letter)):
            if letter[i] == 'D':
                dele_len = dele_len + number[i]
        
        return(dele_len)
        
    def cigar_find_inserlen(self,cigar):
        
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int
        
        #here only want to know mapping length
        #so don't consider strand info and reverse cigar 
        
        inser_len = 0
        for i in range(0,len(letter)):
            if letter[i] == 'I':
                inser_len = inser_len + number[i]
        
        return(inser_len)
        
    def find_microhomo_from_indel(self,n1,n2,md,cigar,sequence):
        microhomo_seq = ''

        #decode cigar and find max dele and extract seq after
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        number = list(map(int, number))#convert string to int
        #which D is the D cross cutsite
        if n2 != 0:
            #find D after seq position
            dele_len = 0
            dele_sum = 0 
            for i in range(0,n1):
                if letter[i] == 'D':
                    dele_len = number[i]
                    sum1 = sum(number[0:i]) - dele_sum
                    dele_seq_after_start = sum1 + 1 
                    dele_sum = dele_sum + sum(number[i:i+1])
            if dele_len > 0:
                dele_seq_after_end = dele_seq_after_start + dele_len - 1
                dele_seq_after = sequence[dele_seq_after_start-1:dele_seq_after_end]

                #decode md and find dele sequence
                md_letter = re.findall('\D', md)
                md0 =''
                md_sub = md0.join(md_letter)
                md_sub_list = md_sub.split('^')
                md_sub_list = md_sub_list[1:]
                #the n string is the D cross cutsite
                dele_seq = md_sub_list[n2-1]
                md_dele_len = md_sub_list[n2-1]

                #find microhomo seq
                micro_mark = 0
                mark = 0 
                for i in range(0,dele_len):
                    if dele_seq[i] == dele_seq_after[i]:
                        mark = mark + 1
                    if dele_seq[i] != dele_seq_after[i]:
                        micro_mark = i
                        break
                if micro_mark == 0:
                    micro_mark = mark
                if micro_mark > 0:
                    microhomo_seq = dele_seq[0:micro_mark]
        return(microhomo_seq)
        
    def define_corss_deletion_len(self,n1,cigar):
        
        #decode cigar and find max dele and extract seq after
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        number = list(map(int, number))#convert string to int
                
        #find max dele
        dele_len = number[n1-1]
                    
        return(dele_len)            
        
    def decide_indel_in_cutsite(self, cigar, cutsite, reference_start):
        # print(cigar)
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        number = list(map(int, number))#convert string to int
        
        match_list = []
        dele_list = []
        inser_list = []
        
        del_start_list = []
        del_end_list = []
        inser_start_list = []
        inser_end_list = []
        
        del_distance_list = []
        inser_distance_list = []
        
        deletion_index_list = []
        insertion_index_list = []
        
        for i in range(0, len(letter)):
            if letter[i] == 'M':
                match_list.append(i)
            if letter[i] == 'D':
                deletion_index_list.append(i)
                if match_list == []:
                    map_sum = 0
                else:
                    index_min = min(match_list)
                    index_max = max(match_list)
                    map_sum = sum(number[index_min:index_max+1])
                del_len = sum(dele_list)
                inser_len = sum(inser_list)
                map_len = map_sum - inser_len
                
                # this 'D' range
                dele_start = reference_start + map_len - 1
                dele_end = dele_start + number[i] + 1
                del_start_list.append(dele_start)
                del_end_list.append(dele_end)
                del_mid = (dele_start + dele_end)/2
                del_distance = abs(del_mid - cutsite)
                del_distance_list.append(del_distance)
                
                dele_list.append(number[i])
                
            if letter[i] == 'I':
                insertion_index_list.append(i)
                if match_list == []:
                    map_sum = 0
                else:
                    index_min = min(match_list)
                    index_max = max(match_list)
                    map_sum = sum(number[index_min:index_max+1])
                del_len = sum(dele_list)
                inser_len = sum(inser_list)
                map_len = map_sum - inser_len
                
                # this 'I' pos
                inser_start = reference_start + map_len - 1
                inser_end = inser_start + 1
                inser_start_list.append(inser_start)
                inser_end_list.append(inser_end)
                inser_mid = (inser_start + inser_end)/2
                inser_distance = abs(inser_mid - cutsite)
                inser_distance_list.append(inser_distance)
                
                inser_list.append(number[i])
        #deletion
        if len(del_distance_list) != 0:
            del_min = min(del_distance_list)
            count = del_distance_list.count(del_min)
            if count > 1 :
                del_distance_list_m = del_distance_list
                del_min_index1 = del_distance_list.index(del_min)
                del_distance_list_m[del_min_index1] = 'N'
                del_min_index2 = del_distance_list_m.index(del_min)
                # print(inser_start_list,inser_end_list,del_min_index1,del_min_index2)
                del_arm1 = (abs(del_start_list[del_min_index1] - del_end_list[del_min_index1])-1)/2
                del_arm2 = (abs(del_start_list[del_min_index2] - del_end_list[del_min_index2])-1)/2
                if del_arm1 > del_arm2:
                    del_min_index = del_min_index2
                else:
                    del_min_index = del_min_index1
            else:
                del_min_index = del_distance_list.index(del_min)
        else:
            del_min = float("inf")# infinity
        #insertion    
        if len(inser_distance_list) != 0:
            inser_min = min(inser_distance_list)
            inser_min_index = inser_distance_list.index(inser_min)
        else:
            inser_min = float("inf")# infinity
        if inser_min == del_min == float("inf"):
            indel_type = "Other"
            start = ''
            end = ''
            index = ''
            n = ''
        else:
            if del_min == inser_min:
                start = del_start_list[del_min_index]
                end = del_end_list[del_min_index]
                del_arm = (abs(start - end)-1)/2
                if del_arm > 0.5:
                    indel_type = "D"
                    start = del_start_list[del_min_index]
                    end = del_end_list[del_min_index]
                    index = del_min_index
                    n = deletion_index_list[del_min_index]
                else:
                    indel_type = "I"
                    start = inser_start_list[inser_min_index]
                    end = inser_end_list[inser_min_index]
                    index = inser_min_index
                    n = insertion_index_list[inser_min_index]
            if del_min < inser_min:
                indel_type = "D"
                start = del_start_list[del_min_index]
                end = del_end_list[del_min_index]
                index = del_min_index
                n = deletion_index_list[del_min_index]
            if del_min > inser_min:
                indel_type = "I"
                start = inser_start_list[inser_min_index]
                end = inser_end_list[inser_min_index]
                index = inser_min_index
                n = insertion_index_list[inser_min_index]
        
        return([indel_type,start,end,n,index])
    
    def decide_letter_in_cutsite(self,letter_type,cigar,cutsite,reference_start,cutoff,extract_seq_mark=''):

        #map_type contain: seq/referece
        #reference map len = total - insertion
        #sequence map len = total - deletion
        
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int
        
        #here only want to know mapping length
        #so don't consider strand info and reverse cigar 
        #map_len = (start M + .. + end M) - D
        
        answer = 'N'
        match_list = []
        dele_list = []
        inser_list = []
        for i in range(0,len(letter)):
            if letter[i] == 'M':
                match_list.append(i)
            
            if letter[i] == 'D':
                
                # if find one 'D' in cutsite then jump out
                if letter_type == 'D':
                    # map reference len
                    index_min = min(match_list)
                    index_max = max(match_list)
                    map_sum = sum(number[index_min:index_max+1])
                    del_len = sum(dele_list)
                    inser_len = sum(inser_list)
                    map_len = map_sum - inser_len
                    
                    # this 'D' range
                    dele_start = reference_start + map_len
                    dele_end = reference_start + map_len + number[i]
                    
                    # condition that 'D' in range of cutsite
                    condition1 = dele_start <= (cutsite + cutoff)
                    condition2 = dele_end >= (cutsite - 0)
                    if condition1 and condition2:
                        answer = 'Y'
                        if extract_seq_mark == 'Deletion':
                            answer = i
                        break
                    else:
                        answer = 'N'

                dele_list.append(number[i])

            if letter[i] == 'I':
                
                # if find one 'I' in cutsite then jump out
                if letter_type == 'I':
                    # map reference len
                    index_min = min(match_list)
                    index_max = max(match_list)
                    map_sum = sum(number[index_min:index_max+1])
                    del_len = sum(dele_list)
                    inser_len = sum(inser_list)
                    map_len = map_sum - inser_len
                    
                    # this 'I' pos
                    inser_pos = reference_start + map_len
                    
                    # condition that 'I' in range of cutsite
                    condition1 = inser_pos <= (cutsite + cutoff)
                    condition2 = inser_pos >= (cutsite - 0)
                    if condition1 and condition2:
                        answer = 'Y'
                        if extract_seq_mark == 'Insertion':
                            answer = i
                        break
                    else:
                        answer = 'N'

                inser_list.append(number[i])
                
        return(answer)
        
    def cigar_seq_pos(self,pos_type,cigar):#calculating map len according to cigar

        #map_type contain: seq/referece
        #reference map len = total - insertion
        #sequence map len = total - deletion
        
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int
        
        #here only want to know mapping length
        #so don't consider strand info and reverse cigar 
        #map_len = (start M + .. + end M) - D
        
        soft_list = []
        match_list = []
        inser_list = []
        dele_list = []
        for i in range(0,len(letter)):
            if letter[i] == 'S':
                soft_list.append(number[i])
            if letter[i] == 'M':
                match_list.append(number[i])
            if letter[i] == 'I':
                inser_list.append(number[i])
            if letter[i] == 'D':
                dele_list.append(number[i])
                
        soft_len = sum(soft_list)
        match_len = sum(match_list)
        inser_len = sum(inser_list)
        dele_len = sum(dele_list)
        if pos_type == 'seq':
            seq_pos = soft_len + match_len + inser_len
        if pos_type == 'ref':
            seq_pos = match_len + dele_len
        
        return(seq_pos)

    def extract_insertion_from_indel(self, sequence, cigar, cutsite, reference_start, index):

        number = re.findall('\d+', cigar)
        number = list(map(int, number))

        if index == 'N':
            inser_seq = ''
        else:
            i = index
            inser_len = number[int(i)]
            cigat_list = re.split('(\d+)',cigar)
            cigar_list_inser = cigat_list[1:i*2+1]
            cigar_for_inser = ''.join(cigar_list_inser)#convert list to string
            inser_seq_start = self.cigar_seq_pos('seq', cigar_for_inser) + 1
            inser_seq_end = inser_seq_start + inser_len - 1
            inser_seq = sequence[inser_seq_start-1:inser_seq_end]

        return(inser_seq)


    def mismatch_pos_of_md(self, md_string, reference_start, cutsite, cutoff):
        
        # find　if there is mismatch in cutsite
        answer = 'N'
        # #remove deletion seuqence
        # partition_list = md_string.partition('^')
        # part1 = partition_list[0]
        # part2 = partition_list[2]
        # split_list = re.split('(\d+)',part2)
        # splist_list_remain = split_list[1:len(split_list)]
        # #list of new md string
        # new_md_string = part1 + part2 #md string remove '^'
        # new_md_string = md_string.replace('^', '')#remove '^'in the string
        new_md_list = re.split('(\d+)',md_string)#convert string to list
        
        su = 0
        digit_list = []
        for i in new_md_list:
            if i.isdigit():#int then put into list
                i = int(i)
                digit_list.append(i)
            elif len(i) == 1:#only one letter then sum up list and check then add 1 to the list
                su = sum(digit_list[0:len(digit_list)])
                condition1 = reference_start + su >= cutsite
                condition2 = reference_start + su <= (cutsite + cutoff)
                if condition1 and condition2:
                    answer = 'Y'
                    break
                else:
                    answer = 'N' 
                digit_list.append(1)
            else:
                digit_list.append(len(i)-1)#multiple letter means deletion then add len(i)-1('^')
        return(answer)
    
    def cigar_map_seq_start(self,cigar):

        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        
        number = list(map(int, number))#convert string to int

        for i in range(0,len(letter)):
            if letter[i] == 'M':
                map_seq_start = sum(number[0:i]) + 1
                break
        return(map_seq_start)

    def generate_indel_tab_file(self):
        
        os.system("mkdir raw")
        bam_file = pysam.AlignmentFile(self.bam_sort, 'rb')
        indel_tab = open("indel/" + self.basename + "_indel.tab", "w")
        transloc_tab = open("indel/" + self.basename + "_indel_mut.tab", "w")
        deletion_tab = open("indel/" + self.basename + "_deletion_mut.tab", "w")
        insertion_tab = open("indel/" + self.basename + "_insertion_mut.tab", "w")
        substitution_tab = open("raw/"+self.basename + "_Substitution.tab", "w")
        germline_tab = open("raw/"+self.basename + "_Germline.tab", "w")
        discard_tab = open("indel/" + self.basename + "_discard.tab", "w")
        
        #all reads
        indel_tab.write('Qname'+"\t"+\
                        'Start'+"\t"+\
                        'End'+"\t"+\
                        'Rname'+"\t"+\
                        'Qseq'+"\t"+\
                        'Strand'+"\t"+\
                        'Cigar'+"\t"+\
                        'MDstring'+"\t"+\
                        'Deletion'+"\t"+\
                        'Decide_deletion'+"\t"+\
                        'Insertion'+"\t"+\
                        'Decide_insertion'+"\t"+\
                        'Substitution_num'+"\t"+\
                        'Decide_Substitution'+"\t"+\
                        'Deletion_start'+"\t"+\
                        'Deletion_end'+"\t"+\
                        'Insertion_seq'+"\t"+\
                        'Microhomolog'+"\t"+\
                        'Prey_MQ'+"\n")
        #deletion or insertion cross cutsite
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
        #deletion
        deletion_tab.write('Qname'+"\t"+\
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
        #only insertion around cutsite
        insertion_tab.write('Qname'+"\t"+\
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
        #substitution around cutstie
        substitution_tab.write('Qname'+"\t"+\
                                'Cigar'+"\t"+\
                                'MDstring'+"\t"+\
                                'Position'+"\t"+\
                                'Sequence'+"\t"+\
                                'Prey_MQ'+"\n")
        #germline reads
        germline_tab.write('Qname'+"\t"+\
                                'Cigar'+"\t"+\
                                'MDstring'+"\t"+\
                                'Position'+"\t"+\
                                'Sequence'+"\t"+\
                                'Prey_MQ'+"\n")
        
        for read in bam_file:
            # print(read.query_name)
            cutsite = self.cutsite
            reference_start = read.reference_start + 1 #pysam reference start from 0
            reference_end = read.reference_end #pysam reference end is exactly the end base
            prey_mapqual = read.mapping_quality
            if read.is_reverse: 
                strand = '-'
                read_check = cutsite - reference_start
            else: 
                strand = '+'
                read_check = reference_end - cutsite
            if read_check <= 10:
                discard_tab.write(read.query_name+"\n")#sonication
            # reads have longer sonication
            else:
                #count insertion and deletion number
                if any('MD' == tg[0] for tg in read.get_tags()):
                    md_string = read.get_tag('MD')
                    delelen = str(self.cigar_find_delelen(read.cigarstring))
                    inserlen = str(self.cigar_find_inserlen(read.cigarstring))
                #count mismatch number
                if any('NM' == tg[0] for tg in read.get_tags()):
                    nm_num = str(read.get_tag('NM'))
                    mismatch_num = int(nm_num) - int(delelen) - int(inserlen)
                
                ###~~~~~~~ classify indel ~~~~~~~~~####
                # cutoff = 5
                # decide_insertion = self.decide_letter_in_cutsite('I',read.cigarstring,self.cutsite,reference_start,cutoff)
                # decide_mismatch = self.mismatch_pos_of_md(md_string, reference_start, self.cutsite, cutoff)
                # decide_deletion =  self.decide_letter_in_cutsite('D',read.cigarstring,self.cutsite,reference_start,cutoff)
                #indel info list: "I"/"D",start,end,n,index
                indel_info_list = self.decide_indel_in_cutsite(read.cigarstring, self.cutsite, reference_start)
                dele_junc_start = indel_info_list[1]
                dele_junc_end = indel_info_list[2]
                indel_type = indel_info_list[0]
                n = indel_info_list[3]
                index = indel_info_list[4]
                #~~~~~~~~ General Info ~~~~~~~~~~#
                qname = read.query_name
                rname = read.reference_name
                bait_rname = read.reference_name
                prey_rname = read.reference_name
                bait_strand = strand
                prey_strand = strand
                sequence = read.query_sequence
                microhomolog = ''
                insertion = ''
                if strand == '-':
                    bait_start = dele_junc_end
                    bait_end = read.reference_end
                    prey_start = reference_start
                    prey_end = dele_junc_start
                    junction = prey_end
                else:
                    bait_start = reference_start
                    bait_end = dele_junc_start
                    prey_start = dele_junc_end
                    prey_end = read.reference_end
                    junction = prey_start
                #~~~~~~~~ class:Germline/substitution or Deletion/Insertion ~~~~~~~~~~#
                #~~~~~~~~ Germline/substitution ~~~~~~~~~~#
                if indel_type == "Other":
                    # print(md_string)
                    if strand == '-':
                        bait_start = reference_end
                    else:
                        bait_start = reference_start
                    # no any mutation
                    if md_string.isdigit():
                        germline_tab.write(qname+"\t"+\
                                            read.cigarstring+"\t"+\
                                            md_string+"\t"+\
                                            str(reference_start)+"\t"+\
                                            sequence+"\t"+\
                                            str(prey_mapqual)+"\n")
                    # substitution
                    else:
                        substitution_tab.write(qname+"\t"+\
                                                read.cigarstring+"\t"+\
                                                md_string+"\t"+\
                                                str(reference_start)+"\t"+\
                                                sequence+"\t"+\
                                            str(prey_mapqual)+"\n")
                #~~~~~~~~ Deletion/Insertion ~~~~~~~~~~#
                else:
                    if indel_type == "D":
                        ### print microhomo sequence ###
                        microhomo = ''
                        i = n
                        #which letter is the D cross cutsite
                        n1 = i + 1
                        #which D is the D cross cutsite
                        letter = re.findall('\D', read.cigarstring)
                        number = re.findall('\d+', read.cigarstring)
                        number = list(map(int, number))
                        n2 = 0 
                        for j in letter[0:i+1]:
                            if j == 'D':
                                n2 = n2 + 1
                        microhomo = self.find_microhomo_from_indel(n1,n2,md_string,read.cigarstring,read.query_sequence)
                        microhomolog = microhomo
                        if strand == '+':
                            bait_end_cor = bait_end + len(microhomolog)
                            prey_end_cor = prey_end
                        else:
                            prey_end_cor = prey_end + len(microhomolog)
                            bait_end_cor = bait_end
                        # if bait_end_cor - bait_end > cutsite - bait_end:
                        #     print(read)
                        
                        # sequence info
                        if read.is_reverse:
                            read_strand = '-'
                            sequence = Seq(read.query_sequence)
                            sequence = str(sequence.complement())
                            Qstart = self.cigar_map_seq_start(read.cigarstring)
                            Qend = Qstart + abs(bait_end_cor - bait_start)
                            B_Qstart = Qend - len(microhomolog) - 1
                            B_Qend = self.cigar_seq_pos('seq', read.cigarstring)
                            Qlen = len(sequence)
                        else:
                            read_strand = '+'
                            sequence = read.query_sequence
                            B_Qstart = self.cigar_map_seq_start(read.cigarstring)
                            B_Qend = B_Qstart + abs(bait_end_cor - bait_start)
                            Qstart = B_Qend - len(microhomolog) - 1
                            Qend = self.cigar_seq_pos('seq', read.cigarstring)
                            Qlen = len(sequence)
                        
                        # print(qname,B_Qstart,B_Qend,Qstart,Qend)
                        
                        deletion_tab.write(qname+"\t"+\
                                            str(bait_rname)+"\t"+\
                                            str(bait_strand)+"\t"+\
                                            str(bait_start)+"\t"+\
                                            str(bait_end_cor)+"\t"+\
                                            prey_rname+"\t"+\
                                            prey_strand+"\t"+\
                                            str(prey_start)+"\t"+\
                                            str(prey_end_cor)+"\t"+\
                                            rname+"\t"+\
                                            prey_strand+"\t"+\
                                            str(junction)+"\t"+\
                                            sequence+"\t"+\
                                            str(B_Qstart)+"\t"+\
                                            str(B_Qend)+"\t"+\
                                            str(Qstart)+"\t"+\
                                            str(Qend)+"\t"+\
                                            str(Qlen)+"\t"+\
                                            str(insertion)+"\t"+\
                                            str(microhomolog)+"\t"+\
                                            str(prey_mapqual)+"\n")
                    #~~~~~~~~ Insertion specific ~~~~~~~~~~#
                    if indel_type == "I":
                        #### extract insertion seq ###
                        # print(read.query_sequence, read.cigarstring, self.cutsite, reference_start,index)
                        inser_seq = self.extract_insertion_from_indel(read.query_sequence, read.cigarstring, self.cutsite, reference_start, n)
                        insertion = inser_seq
                        
                        # sequence info
                        if read.is_reverse:
                            read_strand = '-'
                            sequence = Seq(read.query_sequence)
                            sequence = str(sequence.complement())
                        else:
                            read_strand = '+'
                            sequence = read.query_sequence
                            
                        B_Qstart = self.cigar_map_seq_start(read.cigarstring)
                        B_Qend = B_Qstart + abs(bait_end - bait_start)
                        Qstart = B_Qend + len(insertion) + 1
                        Qend = self.cigar_seq_pos('seq', read.cigarstring)
                        Qlen = len(sequence)
                        
                        insertion_tab.write(qname+"\t"+\
                                            str(bait_rname)+"\t"+\
                                            str(bait_strand)+"\t"+\
                                            str(bait_start)+"\t"+\
                                            str(bait_end)+"\t"+\
                                            prey_rname+"\t"+\
                                            prey_strand+"\t"+\
                                            str(prey_start)+"\t"+\
                                            str(prey_end)+"\t"+\
                                            rname+"\t"+\
                                            prey_strand+"\t"+\
                                            str(junction)+"\t"+\
                                            sequence+"\t"+\
                                            str(B_Qstart)+"\t"+\
                                            str(B_Qend)+"\t"+\
                                            str(Qstart)+"\t"+\
                                            str(Qend)+"\t"+\
                                            str(Qlen)+"\t"+\
                                            str(insertion)+"\t"+\
                                            str(microhomolog)+"\t"+\
                                            str(prey_mapqual)+"\n")
                        bait_end_cor = bait_end
                        prey_end_cor = prey_end
                    ###. all class 1 reads output. ###
                    # indel_tab.write(read.query_name+"\t"+\
                    #                 str(reference_start)+"\t"+\
                    #                 str(read.reference_end)+"\t"+\
                    #                 str(read.reference_name)+"\t"+\
                    #                 read.query_sequence+"\t"+\
                    #                 strand+"\t"+\
                    #                 read.cigarstring+"\t"+\
                    #                 md_string+"\t"+\
                    #                 delelen+"\t"+\
                    #                 decide_deletion+"\t"+\
                    #                 inserlen+"\t"+\
                    #                 decide_insertion+"\t"+\
                    #                 str(mismatch_num)+"\t"+\
                    #                 decide_mismatch+"\t"+\
                    #                 str(dele_junc_start)+"\t"+\
                    #                 str(dele_junc_end)+"\t"+\
                    #                 str(insertion)+"\t"+\
                    #                 str(microhomolog)+"\n")
                    transloc_tab.write(qname+"\t"+\
                                        str(bait_rname)+"\t"+\
                                        str(bait_strand)+"\t"+\
                                        str(bait_start)+"\t"+\
                                        str(bait_end_cor)+"\t"+\
                                        prey_rname+"\t"+\
                                        prey_strand+"\t"+\
                                        str(prey_start)+"\t"+\
                                        str(prey_end_cor)+"\t"+\
                                        rname+"\t"+\
                                        prey_strand+"\t"+\
                                        str(junction)+"\t"+\
                                        sequence+"\t"+\
                                        str(B_Qstart)+"\t"+\
                                        str(B_Qend)+"\t"+\
                                        str(Qstart)+"\t"+\
                                        str(Qend)+"\t"+\
                                        str(Qlen)+"\t"+\
                                        str(insertion)+"\t"+\
                                        str(microhomolog)+"\t"+\
                                        str(prey_mapqual)+"\n")
        bam_file.close()
        indel_tab.close()
        deletion_tab.close()
        insertion_tab.close()
        transloc_tab.close()
        substitution_tab.close()
        germline_tab.close()
        discard_tab.close()

        #cutoff indel
        indel_transloc = pd.read_table("indel/" + self.basename + "_indel_mut.tab",sep = '\t')
        print("indel/" + self.basename + "_indel_mut.tab",indel_transloc['Qname'].count())
        if self.primer_strand == '+':
            indel_s = 'Bait_end'
            indel_e = 'Prey_start'
        else:
            indel_s = 'Prey_end'
            indel_e = 'Bait_start'
        print(indel_s,indel_e,self.cutsite)
        cutoff_condition = (indel_transloc[indel_s] <= self.cutsite + 10) & \
                            (indel_transloc[indel_e] >= self.cutsite - 10)
        indel_transloc_cutoff = indel_transloc[cutoff_condition]
        indel_transloc_cutoff.to_csv("raw/"+self.basename + "_indel_cutoff.tab", header = True, sep = '\t', index=False)
        germline_addup = indel_transloc[~(cutoff_condition)]
        
        germline_addup.to_csv("raw/"+self.basename + "_Germline_addup_indel.tab", header = True, sep = '\t', index=False)
        
        # merge two transloc file !!!!
        transloc = pd.read_table("transloc/" + self.basename + "_mut.tab",sep = '\t')
        print("transloc/" + self.basename + "_mut.tab",transloc['Qname'].count())
        transloc = transloc.append([indel_transloc_cutoff])
        transloc.to_csv("raw/"+self.basename + "_SID_all.tab", header = True, sep = '\t', index=False)
        print("raw/"+self.basename + "_SID_all.tab",transloc['Qname'].count())
        
        # classify into deletions,insertions,SV: ####***!!!!****###
        if self.primer_strand == '+':
            sd_con = transloc['Bait_end'] < transloc['Prey_start']
        else:
            sd_con = transloc['Bait_start'] > transloc['Prey_end']
        #Deletion
        
        
        Deletions_condition = (transloc['Bait_rname'] == transloc['Prey_rname']) & \
                             (transloc['Bait_strand'] == transloc['Prey_strand']) & \
                             (sd_con) & \
                             (transloc['Insertion'].isna() ) #!!!!
        Deletions = transloc[Deletions_condition]
        Deletions.to_csv("raw/"+self.basename + "_Deletion.tab", header = True, sep = '\t', index=False)
        #Insertion(indel)
        Insertion_condition = (transloc['Bait_rname'] == transloc['Prey_rname']) & \
                             (transloc['Bait_strand'] == transloc['Prey_strand']) & \
                             (sd_con) & \
                             (transloc['Insertion'].str.len() > 0)
        Insertion = transloc[Insertion_condition]
        Insertion.to_csv("raw/"+self.basename + "_Insertion.tab", header = True, sep = '\t', index=False)
        #SV
        SVs = transloc[~transloc['Qname'].isin(Deletions['Qname'])]
        SVf = SVs[~SVs['Qname'].isin(Insertion['Qname'])]
        SVf.to_csv("raw/"+self.basename + "_SV.tab", header = True, sep = '\t', index=False)
        
        #Insertions(SV)
        Insertion_condition_sv = (SVf['Insertion'].str.len() > 0)
        Insertion_sv = SVf[Insertion_condition_sv]
        Insertion_sv.to_csv("raw/"+self.basename + "_Insertions_inSV.tab", header = True, sep = '\t', index=False)
        #ALL Insertions
        all_insertion = Insertion.append(Insertion_sv)
        all_insertion.to_csv("raw/"+self.basename + "_All_Insertions.tab", header = True, sep = '\t', index=False)
    # def generate_other_tab(self):

    #     data = pd.read_table("indel/" + self.basename + "_indel.tab",sep = '\t',index_col=False,low_memory=False,names = ['Qname',
    #                                                                                                                     'Start',
    #                                                                                                                     'End',
    #                                                                                                                     'Rname',
    #                                                                                                                     'Qseq',
    #                                                                                                                     'Strand',
    #                                                                                                                     'Cigar',
    #                                                                                                                     'MDstring',
    #                                                                                                                     'Deletion',
    #                                                                                                                     'Decide_deletion',
    #                                                                                                                     'Insertion',
    #                                                                                                                     'Decide_insertion',
    #                                                                                                                     'Mismatch_num',
    #                                                                                                                     'Decide_mismatch',
    #                                                                                                                     'Deletion_start',
    #                                                                                                                     'Deletion_end',
    #                                                                                                                     'Insertion_seq',
    #                                                                                                                     'Microhomolog'])

def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='denfine_indel 1.0')
    
    kwargs = {'basename':args['<basename>'], 'cutsite':args['<cutsite>'], 'primer_strand':args['<primer_strand>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    
    define_indel = Define_indel(**kwargs)
    define_indel.generate_indel_tab_file()
    # define_indel.generate_other_tab()

    print("\ndefine_indel.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
