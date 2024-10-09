#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.4.23
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
    rmb_dedup   <basename>  <barcode_length> <adapter>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script help to dedup reads according its barcode attatched, using levenshstein 
distance method. A bam file of adapter alignment is needed for this script. Both adapter 
and barcode are allowed 2 mismatches.

Input file: adapter alignment bam file / Output file: unique reads's query name list

Author: Mengzhu LIU
Last Update:2019.5.9

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
    #######################################
    
import os
import pysam
import gzip
from time import time
from docopt import docopt
import pandas as pd

class Dedup(object):
    
    def __init__(self, basename=None, barcode_length=None, adapter=None):
        
        self.basename = basename
        self.barcode_length = int(barcode_length)
        
        # init adapter 
        adapter_f = "adapter/adapter.fa" 
        self.load(adapter_f)  
        adapt_f = pysam.FastaFile(adapter_f)
        adapter = adapt_f.fetch(reference = 'adapter')
        adapt_f.close()
        print("[PEM-Q] adapter_seqeunce: " + adapter)
        self.adapter = adapter
        
        # init bam file
        fastq_check = "flash_out/" + basename + ".merge.fastq.gz"
        adapter_sam_check = "barcode/" + basename + "_check.adpt.sam"
        adapter_bam_check = "barcode/" + basename + "_check.adpt.bam"
        adapter_bam_check_sort = "barcode/" + basename + "_check.adpt.sort.bam"
        
        adapter_bam = "bwa_align/" + basename + "_sti.adpt.sort.bam"#get adapter bam file
        primer_bam = "primer/" + basename + "_sti.p.sort.bam"#get primer filter stitch bam file
        dedup_bam = "barcode/" + basename + "_sti.dedup.bam"
        dedup_bam_sort = "barcode/" + basename + "_sti.dedup.sort.bam"
        dup_bam = "barcode/" + basename + "_sti.dup.bam"
        dup_bam_sort = "barcode/" + basename + "_sti.dup.sort.bam"
        filter_bam = "barcode/" + basename + "_sti.filter.bam"
        filter_bam_sort = "barcode/" + basename + "_sti.filter.sort.bam"
        
        
        self.fastq_check = fastq_check
        self.adapter_sam_check = adapter_sam_check
        self.adapter_bam_check = adapter_bam_check
        self.adapter_bam_check_sort = adapter_bam_check_sort
        
        self.adapter_bam = adapter_bam
        self.primer_bam = primer_bam
        self.dedup_bam = dedup_bam
        self.dedup_bam_sort = dedup_bam_sort
        self.dup_bam = dup_bam
        self.dup_bam_sort = dup_bam_sort
        
        self.filter_bam = filter_bam
        self.filter_bam_sort = filter_bam_sort
        
        
        
        
        # self.outputdir = self.basename
    
    def load(self, f_file):
        
        if not os.path.exists(f_file):
            raise ValueError('[PEM-Q] The {} file does not exist.'.format(f_file))
    

    def extract_barcode(self):                
            
        # only keep query with barcode length >= barcode_length-2
        os.system("mkdir barcode")
        barcode_list = open("barcode/"+self.basename+"_barcode_list.txt", "w")
        barcode_input_file = self.basename+'_I2.fq.gz'
        if os.path.exists(barcode_input_file):
            if barcode_input_file.endswith('.gz'):
                f1_in = gzip.open(barcode_input_file,'rt')
            else:
                f1_in = open(barcode_input_file,'rt')
            while (1):
                f1_id_line   = f1_in.readline().strip()
                f1_seq_line  = f1_in.readline().strip()
                f1_plus_line = f1_in.readline()
                f1_qual_line = f1_in.readline().strip()
                if not f1_qual_line : break

                query_name = f1_id_line.split(" ")[0][1:]
                barcode_seq = f1_seq_line
                bar_length = len(f1_seq_line)
                if bar_length <= self.barcode_length and bar_length >=(self.barcode_length - 2):#barcode allow 2 deletion
                    barcode_list.write(query_name + "\t" + str(barcode_seq)+"\n")




        else:
            bam_file = pysam.AlignmentFile(self.adapter_bam, "rb")
            
            for read in bam_file.fetch():
                if read.query_alignment_length >= (len(self.adapter) - 2):#aligned adapter length allow 2 deletion
                    bar_start = read.query_alignment_end + 1
                    bar_end = read.query_length
                    bar_length = bar_end - bar_start + 1
                    if bar_length <= self.barcode_length and bar_length >=(self.barcode_length - 2):#barcode allow 2 deletion
                        barcode_seq = read.query_sequence[(bar_start-1) : bar_end]
                        barcode_list.write(str(read.query_name) + "\t" + str(barcode_seq)+"\n")            
        barcode_list.close()
    
    def barcode_dedup(self):  
        
        data = pd.read_csv("barcode/"+self.basename+"_barcode_list.txt",sep = '\t',names = ["Qname", "Barcode"])
        # get barcode length
        length = data['Barcode'].astype('str')
        length = data['Barcode'].str.len()
        data['Length'] = length
        # get barcode frequency
        data['Freq'] = data.groupby('Barcode')['Barcode'].transform('count')
        data = data.sort_values(by=['Freq','Barcode','Length'], ascending=False)
        data.to_csv("barcode/"+self.basename+"_barcode_sort.txt", header = True, sep = '\t', index=False,
           columns = [u'Qname',u'Barcode',u'Freq',u'Length'])
           
        #also generate duplicated barcode
        datb = data[data.duplicated(['Barcode'])]
        datb = datb.reset_index(drop=True)
        #unique barcode
        data = data.drop_duplicates(['Barcode'], keep = 'first' )#rough dedup
        data = data.reset_index(drop=True)#reset index
        
        
        # freq max and min value
        num_max = data['Freq'].max()
        num_min = data['Freq'].min()
        # print(num_max,num_min)
        data.to_csv("barcode/"+self.basename+"_barcode_uniq.txt", header = True, sep = '\t', index=False,
           columns = [u'Qname',u'Barcode',u'Freq',u'Length'])
        datb.to_csv("barcode/"+self.basename+"_barcode_dup.txt", header = True, sep = '\t', index=False,
           columns = [u'Qname',u'Barcode',u'Freq',u'Length'])
        
        # # drop freq < freq_max*0.05 query
        # data = data[data['Freq'] >= num_max*0.05]
        # data.to_csv("barcode/"+self.basename+"_barcode_uniq.txt", header = True, sep = '\t', index=False,
        #    columns = [u'Qname',u'Barcode',u'Freq',u'Length'])
           
    # def generate_dedup_bam(self):
    #
    #     # get uniq qname list
    #     data = pd.read_csv("barcode/"+self.basename+"_barcode_uniq.txt",sep = '\t',names = [u'Qname',u'Barcode',u'Freq',u'Length'],low_memory=False)
    #     uniq_qname_list = data['Qname'].tolist()
    #
    #     #generate dedup bam
    #     primer_bam = pysam.AlignmentFile(self.primer_bam, 'rb')
    #     dedup_bam = pysam.AlignmentFile(self.dedup_bam, "wb", template=primer_bam)
    #
    #     # #too slow
    #     # for read in primer_bam:
    #     #     # print(read.query_name)
    #     #     if read.query_name in uniq_qname_list:
    #     #         dedup_bam.write(read)
    #     # primer_bam.close()
    #     # dedup_bam.close()
    #     # pysam.sort("-o", self.dedup_bam_sort, self.dedup_bam)
    #
    #     #index bam name by pysam to generate dedup bam
    #     name_indexed = pysam.IndexedReads(primer_bam)
    #     name_indexed.build()
    #     for name in uniq_qname_list:
    #             try:
    #                 name_indexed.find(name)
    #             except KeyError:
    #                 pass
    #             else:
    #                 iterator = name_indexed.find(name)
    #                 for x in iterator:
    #                     dedup_bam.write(x)
    #     primer_bam.close()
    #     dedup_bam.close()
    #     pysam.sort("-o", self.dedup_bam_sort, self.dedup_bam)
        
    # def generate_dup_bam(self):
    #
    #     # get uniq qname list
    #     data = pd.read_csv("barcode/"+self.basename+"_barcode_dup.txt",sep = '\t',names = [u'Qname',u'Barcode',u'Freq',u'Length'],low_memory=False)
    #     uniq_qname_list = data['Qname'].tolist()
    #
    #     #generate dup bam
    #     primer_bam = pysam.AlignmentFile(self.primer_bam, 'rb')
    #     dup_bam = pysam.AlignmentFile(self.dup_bam, "wb", template=primer_bam)
    #
    #     # #too slow
    #     # for read in primer_bam:
    #     #     # print(read.query_name)
    #     #     if read.query_name in uniq_qname_list:
    #     #         dup_bam.write(read)
    #     # primer_bam.close()
    #     # dup_bam.close()
    #     # pysam.sort("-o", self.dup_bam_sort, self.dup_bam)
    #
    #     #index bam name by pysam to generate dup bam
    #     name_indexed = pysam.IndexedReads(primer_bam)
    #     name_indexed.build()
    #     for name in uniq_qname_list:
    #             try:
    #                 name_indexed.find(name)
    #             except KeyError:
    #                 pass
    #             else:
    #                 iterator = name_indexed.find(name)
    #                 for x in iterator:
    #                     dup_bam.write(x)
    #     primer_bam.close()
    #     dup_bam.close()
    #     pysam.sort("-o", self.dup_bam_sort, self.dup_bam)
    
    def filter_multiple_adapter(self):
        
        chek_file = self.adapter_bam_check_sort
        
        if not os.path.exists(chek_file):
            
            if not self.fastq_check:
                raise ValueError('merged fastq is needed for adapter check alignment.')
            if not os.path.exists("adapter/adapter.fa"):
                raise ValueError('adapter/adapter.fa is needed for adapter check alignment.')
                    
            #alignment
            print("[PEM-Q]  align to check adapter...")
                              
            cmd = "bwa mem -t 8 -k 5 -L 0 -T 14 adapter/adapter {} > {} 2>barcode/bwa_align_adapter.log".format(self.fastq_check, 
                                                              self.adapter_sam_check)
            os.system(cmd)
            print("[PEM-Q] "+cmd)
            cmd = "samtools view -S -b -h {} > {} \
                   && samtools sort {} > {} \
                   && samtools index {}".format(self.adapter_sam_check, self.adapter_bam_check, self.adapter_bam_check, self.adapter_bam_check_sort, self.adapter_bam_check_sort)
            print("[PEM-Q]  sort and index bam...")
            os.system(cmd)
        else:
            print("[PEM-Q]  adapter check alignment file exist, jump...")  
        
        #keep record of multiple adapters
        multiple_adapt = open("barcode/"+self.basename+"_multiple_adapt.txt", "w")#生成储存含有multiple adapter的文件
        clean_adapt = open("barcode/"+self.basename+"_clean_adapt.txt", "w")#生成储存没有multiple adapter的文件
        bam_file = pysam.AlignmentFile(self.adapter_bam_check_sort, "rb")#读入sort后的bam文件
        multiple_adapt_list = []
        clean_adapt_list = []
        for read in bam_file:
            condition1 = any('SA' == tg[0] for tg in read.get_tags())#判断该read是否有多条supplementary的比对结果
            if condition1:
                multiple_adapt_list.append(read.query_name)
                multiple_adapt.write(read.query_name+"\n")
            else:
                clean_adapt_list.append(read.query_name)
                clean_adapt.write(read.query_name+"\n")
        multiple_adapt.close()
        clean_adapt.close()
        bam_file.close()
        #remove reads with multiple adapters
        primer_bam = pysam.AlignmentFile(self.primer_bam, 'rb')
        dedup_bam_sort = primer_bam
        filter_bam = pysam.AlignmentFile(self.filter_bam, "wb", template=dedup_bam_sort)
        name_indexed = pysam.IndexedReads(dedup_bam_sort)
        name_indexed.build()
        for name in clean_adapt_list:
                try:
                    name_indexed.find(name)
                except KeyError:
                    pass
                else:
                    iterator = name_indexed.find(name)
                    for x in iterator:
                        filter_bam.write(x)
        dedup_bam_sort.close()
        filter_bam.close()
        pysam.sort("-o", self.filter_bam_sort, self.filter_bam)
        pysam.index(self.filter_bam_sort)
        primer_bam.close()
        return()

def main():
    start_time = time()
    
    args = docopt(__doc__,version='rmb_dedup 1.0')
    
    kwargs = {'basename':args['<basename>'], 'barcode_length':args['<barcode_length>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] barcode_length: ' + str(kwargs['barcode_length']))
    
    dedup = Dedup(**kwargs)
    dedup.extract_barcode()
    dedup.barcode_dedup()
    dedup.filter_multiple_adapter()

    print("\nrmb_dedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
