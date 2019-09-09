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
    define_indel    <basename>  <cutsite>

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

class Define_indel(object):
    
    def __init__(self, basename=None,cutsite=None):
        
        self.basename = basename
        bam_sort = "indel/" + basename + "_indel.sort.bam"
        self.bam_sort = bam_sort
        self.cutsite = int(cutsite)
    
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
        
    def find_microhomo_from_indel(self,md,cigar,sequence):
        
        microhomo_seq = ''

        #decode cigar and find max dele and extract seq after
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        number = list(map(int, number))#convert string to int
        
        #find max dele from cigar and decide after seq position
        dele_len = 0
        dele_sum = 0 
        for i in range(0,len(letter)):
            if letter[i] == 'D':
                if number[i] > dele_len:
                    dele_len = number[i]
                    sum1 = sum(number[0:i]) - dele_sum
                    dele_seq_after_start = sum1 + 1 
                    dele_sum = sum(number[i:i+1])
        
        if dele_len > 0:
            dele_seq_after_end = dele_seq_after_start + dele_len - 1
            dele_seq_after = sequence[dele_seq_after_start-1:dele_seq_after_end]

            #decode md and find max dele sequence
            md_letter = re.findall('\D', md)
            md0 =''
            md_sub = md0.join(md_letter)
            md_sub_list = md_sub.split('^')
            md_dele_len = 0
            for i in range(0,len(md_sub_list)):
                if len(md_sub_list[i]) > md_dele_len:
                    md_dele_len = len(md_sub_list[i])
                    dele_seq = md_sub_list[i]
                    
            #fine microhomo seq
            micro_mark = 0
            for i in range(0,dele_len):
                if dele_seq[i] != dele_seq_after[i]:
                    micro_mark = i
                    break
            if micro_mark > 0:
                microhomo_seq = dele_seq[0:micro_mark]
        
        return(microhomo_seq)
            
            
        
    def define_deletion_junction(self,reference_start,cigar):
        
        #decode cigar and find max dele and extract seq after
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        number = list(map(int, number))#convert string to int
        
        clip_len = 0 
        for i in range(0,len(letter)):
            if letter[i] == 'M':
                if i >=1:
                    clip_len = number[i-1]
                break
                
        define_deletion_junction = ''
        #find max dele from cigar and decide after reference pos
        max_dele_len = 0
        dele_sum = 0 
        inser_len = 0
        for i in range(0,len(letter)):
            if letter[i] == 'I':
                inser_len = inser_len + number[i]
            if letter[i] == 'D':
                if number[i] > max_dele_len:
                    max_dele_len = number[i]
                    dele_ref_after_start = sum(number[0:i+1]) - inser_len - clip_len + 1
                    
        if max_dele_len > 0:
            define_deletion_junction = reference_start + dele_ref_after_start - 1
            
        return(define_deletion_junction)
        
    def define_max_deletion_len(self,cigar):
        
        #decode cigar and find max dele and extract seq after
        letter = re.findall('\D', cigar)
        number = re.findall('\d+', cigar)
        number = list(map(int, number))#convert string to int
                
        #find max dele
        max_dele_len = 0 
        for i in range(0,len(letter)):
            if letter[i] == 'D':
                if number[i] > max_dele_len:
                    max_dele_len = number[i]
                    
        return(max_dele_len)            
        
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
                    condition2 = dele_end >= (cutsite - cutoff)
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
                    condition2 = inser_pos >= (cutsite - cutoff)
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

    def extract_insertion_from_indel(self, sequence, cigar, cutsite, reference_start, cutoff):

        number = re.findall('\d+', cigar)
        number = list(map(int, number))

        i = self.decide_letter_in_cutsite('I',cigar,cutsite,reference_start,cutoff,'Insertion')
        if i == 'N':
            inser_seq = ''
        else:
            inser_len = number[int(i)]
            cigat_list = re.split('(\d+)',cigar)
            cigar_list_inser = cigat_list[1:i*2+1]
            cigar_for_inser = ''.join(cigar_list_inser)#convert list to string
            inser_seq_start = self.cigar_seq_pos('seq', cigar_for_inser) + 1
            inser_seq_end = inser_seq_start + inser_len - 1
            inser_seq = sequence[inser_seq_start-1:inser_seq_end]

        return(inser_seq)

    def extract_insertion_pos(self, sequence, cigar, cutsite, reference_start, cutoff):

        number = re.findall('\d+', cigar)
        number = list(map(int, number))

        i = self.decide_letter_in_cutsite('I',cigar,cutsite,reference_start,cutoff,'Insertion')
        if i == 'N':
            inser_seq = ''
        else:
            inser_len = number[int(i)]
            cigat_list = re.split('(\d+)',cigar)
            cigar_list_inser = cigat_list[1:(i+1)*2+1]
            cigar_for_inser = ''.join(cigar_list_inser)#convert list to string
            inser_seq_start = reference_start + self.cigar_seq_pos('ref', cigar_for_inser)

        return(inser_seq_start)


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
                condition1 = reference_start + su >= (cutsite - cutoff)
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


    def generate_indel_tab_file(self):
        
        bam_file = pysam.AlignmentFile(self.bam_sort, 'rb')
        indel_tab = open("indel/" + self.basename + "_indel.tab", "w")
        transloc_tab = open("indel/" + self.basename + "_indel_transloc.tab", "w")
        insertion_tab = open("indel/" + self.basename + "_insertion_transloc.tab", "w")
        substitution_tab = open("indel/" + self.basename + "_indel_substitution.tab", "w")
        germline_tab = open("indel/" + self.basename + "_indel_germline.tab", "w")
        
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
                        'Microhomolog'+"\n")
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
                            'Insertion'+"\t"+\
                            'Microhomolog'+"\n")
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
                            'Insertion'+"\t"+\
                            'Microhomolog'+"\n")
        #substitution around cutstie
        substitution_tab.write('Qname'+"\t"+\
                                'Cigar'+"\t"+\
                                'MDstring'+"\t"+\
                                'Position'+"\t"+\
                                'Sequence'+"\n")
        #germline reads
        germline_tab.write('Qname'+"\t"+\
                            'Cigar'+"\t"+\
                            'MDstring'+"\t"+\
                            'Sequence'+"\n")
        
        for read in bam_file:
            reference_start = read.reference_start + 1 #pysam reference start from 0
            reference_end = read.reference_end #pysam reference end is exactly the end base
            #condition to keep reads cross cutsite as indel reads
            condition1 = reference_start <= self.cutsite
            condition2 = read.reference_end >= self.cutsite
            if  condition1 and condition2:
                if read.is_reverse: strand = '-'
                else: strand = '+'
                #count insertion and deletion number
                if any('MD' == tg[0] for tg in read.get_tags()):
                    md_string = read.get_tag('MD')
                    delelen = str(self.cigar_find_delelen(read.cigarstring))
                    inserlen = str(self.cigar_find_inserlen(read.cigarstring))
                #count mismatch number
                if any('NM' == tg[0] for tg in read.get_tags()):
                    nm_num = str(read.get_tag('NM'))
                    mismatch_num = int(nm_num) - int(delelen) - int(inserlen)
                #print microhomo sequence
                microhomo = self.find_microhomo_from_indel(md_string,read.cigarstring,read.query_sequence)
                #deletion_junction
                dele_junc_end = self.define_deletion_junction(reference_start,read.cigarstring)
                if dele_junc_end == '':
                    dele_junc_start = ''
                else:
                    max_delelen = self.define_max_deletion_len(read.cigarstring)
                    dele_junc_start = int(dele_junc_end) - int(max_delelen) - 1
                
                #classify indel
                cutoff = 5
                decide_insertion = self.decide_letter_in_cutsite('I',read.cigarstring,self.cutsite,reference_start,cutoff)
                decide_deletion =  self.decide_letter_in_cutsite('D',read.cigarstring,self.cutsite,reference_start,cutoff)
                decide_mismatch = self.mismatch_pos_of_md(md_string, reference_start, self.cutsite, cutoff)
                
                #extract insetion seq
                inser_seq = self.extract_insertion_from_indel(read.query_sequence, read.cigarstring, self.cutsite, reference_start, cutoff)

                indel_tab.write(read.query_name+"\t"+\
                                str(reference_start)+"\t"+\
                                str(read.reference_end)+"\t"+\
                                str(read.reference_name)+"\t"+\
                                read.query_sequence+"\t"+\
                                strand+"\t"+\
                                read.cigarstring+"\t"+\
                                md_string+"\t"+\
                                delelen+"\t"+\
                                decide_deletion+"\t"+\
                                inserlen+"\t"+\
                                decide_insertion+"\t"+\
                                str(mismatch_num)+"\t"+\
                                decide_mismatch+"\t"+\
                                str(dele_junc_start)+"\t"+\
                                str(dele_junc_end)+"\t"+\
                                str(inser_seq)+"\t"+\
                                str(microhomo)+"\n")
                

                # generate transloc tab from indel info
                qname = read.query_name
                rname = read.reference_name
                bait_rname = read.reference_name
                prey_rname = read.reference_name
                bait_strand = strand
                prey_strand = strand
                sequence = read.query_sequence
                insertion = inser_seq
                microhomolog = microhomo

                if decide_deletion == 'Y': #deletion cross cutsite
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
                    transloc_tab.write(qname+"\t"+\
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
                                        str(insertion)+"\t"+\
                                        str(microhomolog)+"\n")
                if decide_deletion == 'N' and decide_insertion == 'Y': # no deletion but has insetion
                    if strand == '-':
                        bait_start = self.extract_insertion_pos(read.query_sequence, read.cigarstring, self.cutsite, reference_start, cutoff)
                        bait_end = read.reference_end
                        prey_start = reference_start
                        prey_end = bait_start - 1
                        junction = prey_end
                    else:
                        bait_start = reference_start
                        bait_end = self.extract_insertion_pos(read.query_sequence, read.cigarstring, self.cutsite, reference_start, cutoff) - 1 
                        prey_start =  bait_end + 1
                        prey_end = read.reference_end
                        junction = prey_start
                    transloc_tab.write(qname+"\t"+\
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
                                        str(insertion)+"\t"+\
                                        str(microhomolog)+"\n")
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
                                        str(insertion)+"\t"+\
                                        str(microhomolog)+"\n")
                if decide_deletion == 'N' and decide_insertion == 'N' and decide_mismatch == 'Y':
                    if strand == '-':
                        bait_start = reference_end
                    else:
                        bait_start = reference_start
                    substitution_tab.write(qname+"\t"+\
                                            read.cigarstring+"\t"+\
                                            md_string+"\t"+\
                                            str(bait_start)+"\t"+\
                                            sequence+"\n")
                if decide_deletion == 'N' and decide_insertion == 'N' and decide_mismatch == 'N':
                    if strand == '-':
                        bait_start = reference_end
                    else:
                        bait_start = reference_start
                    germline_tab.write(qname+"\t"+\
                                            md_string+"\t"+\
                                            read.cigarstring+"\t"+\
                                            str(bait_start)+"\t"+\
                                            sequence+"\n")

                    

        bam_file.close()
        indel_tab.close()
        transloc_tab.close()
        substitution_tab.close()
        germline_tab.close()
        
        # merge two transloc file
        transloc = pd.read_table("transloc/" + self.basename + "_transloc.tab",sep = '\t')
        indel_transloc = pd.read_table("indel/" + self.basename + "_indel_transloc.tab",sep = '\t')
        transloc = transloc.append([indel_transloc])
        transloc.to_csv(self.basename + "_transloc_all.tab", header = True, sep = '\t', index=False)

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
    
    kwargs = {'basename':args['<basename>'], 'cutsite':args['<cutsite>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    
    define_indel = Define_indel(**kwargs)
    define_indel.generate_indel_tab_file()
    # define_indel.generate_other_tab()

    print("\ndefine_indel.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
