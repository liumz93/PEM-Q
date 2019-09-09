#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.7.5

"""substitution_cov

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
    substitution_cov <basename> <chrom> <chrom_size>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .

basename: file basename
chrom: your target chromosome
chrom_size: a genome chromosome size annotation file

Author: Mengzhu LIU
Last Update:2019.7.5

"""

import os
import re
from time import time
from docopt import docopt
import pandas as pd
import numpy as np


def cal_soft_clipping_number(cigar):
    
    cigar_list = re.split('(\d+)',cigar)
    if cigar_list[2] == 'S':
        number = int(cigar_list[1])
    else:
        number = 0
    return(number)
    
def generate_bed_file(basename=None,chrom=None,chrom_size=None):
    
    substitution_file = "indel/" + basename + "_indel_substitution.tab"
    data = pd.read_csv(substitution_file, sep = '\t', index_col=False, low_memory=False)
    sub_bed_file = "indel/" + basename + "_substitution.bed"
    sub_bed = open(sub_bed_file, "w")
    
    mdstring = data["MDstring"]
    count = data["MDstring"].count()
    for i in range(0, count):
        md_string = data["MDstring"][i]
        md_list = re.split('(\d+)',md_string)
        md_string_list = md_list[1:(len(md_list)-1)]
        
        len_list = []
        seq_length = 0
        sub_length = 0
        for j in range(0, len(md_string_list)):
            if md_string_list[j].isdigit():
                len_list.append(int(md_string_list[j]))
                seq_length = sum(len_list)
            else:
                sub_length = sub_length + 1
                seq_length = seq_length + sub_length
                bed_start = (int(data["Position"][i]) + seq_length - 1) - 1
                bed_end = bed_start + 1
                
                clip_number=cal_soft_clipping_number(data["Cigar"][i])
                mismatch_base = data["Sequence"][i][(clip_number+seq_length-1):(clip_number+seq_length)]
                sub_bed.write(chrom+"\t"+\
                              str(bed_start)+"\t"+\
                              str(bed_end)+"\t"+\
                              md_string_list[j]+"\t"+\
                              mismatch_base+"\t"+\
                              data["Sequence"][i]+"\n")
    sub_bed.close()
    
    chrom_size_file = chrom_size
    sub_histo_file = "indel/" + basename + "_substitution_histo.txt"
    cmd = "bedtools genomecov -i {} -g {} -d > {}".format(sub_bed_file,chrom_size_file,sub_histo_file)
    # os.system(cmd)
    
def main():
    
    start_time = time()
    
    args = docopt(__doc__,version='substitution_cov 1.0')
    
    kwargs = {'basename':args['<basename>'],'chrom':args['<chrom>'],'chrom_size':args['<chrom_size>']}
    print('[superCasQ] basename: ' + str(kwargs['basename']))
    print('[superCasQ] chrom: ' + str(kwargs['chrom']))
    print('[superCasQ] chrom_size: ' + str(kwargs['chrom_size']))

    generate_bed_file(**kwargs)

    print("\nsubstitution_cov.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
