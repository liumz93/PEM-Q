#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.9.1
#Last Update:2021.6.25

"""PEM-Q
    
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
Contact: liumz@pku.edu.cn/liu.mengzhu128@gmail.com

Usage:
    PEM-Q.py <genome> <sample> <cutsite> <primer_chr> <primer_start> <primer_end> <primer_strand> <primer> <umi_length>

Options:
-h --help               Show this screen.
-v --version            Show version.
<genome>                reference genome(hg19/hg38/mm9/mm10).
<sample>                sample name of input file(Note: your input file must look like <sample>_R1/2.fq.gz).
<cutsite>               position of cutsite(3' break end of positive strand).
<primer_chr>            chromosome of red primer(eg:chr1).
<primer_start>          start of red primer.
<primer_end>            end of red primer.
<primer_strand>         strand of red primer(+/-).
<primer>                sequence of red primer.
<umi_length>            length of umi

In this program is for PEM-seq data analysis. Pair end target sequencing data is also compatible with this program.

Input file: fastq file 
Output directory: results
Last Update:2021.1.12

"""

import os
import sys
import threading
from time import time
from docopt import docopt

def run_script(sample=None, cutsite=None, genome=None, primer=None, primer_chr=None, primer_start=None, primer_end=None, primer_strand=None, umi_length=17):
    
    start_time = time()
    
    # basename = "LY015B_100"
    # cutsite = "36573328"
    # genome = "hg38"
    # primer = "AGGATCTCACCCGGAACAGC"
    basename = sample

    os.system("mkdir log")

    print("######## 01 Reads alignment... ########")
    cmd = "align_make_v5.1.py {} {}_R1.fq.gz {}_R2.fq.gz -a CCACGCGTGCTCTACA -p {} -r {} -s {} -e {} -d {}".format(genome,basename,basename,primer,primer_chr,primer_start,primer_end,primer_strand)
    print(cmd)
    os.system(cmd)

    print("######## 02 Barcode Extract... ########")
    cmd = "rmb_dedup_v4.py {} {} CCACGCGTGCTCTACA".format(basename, umi_length)
    print(cmd)
    os.system(cmd)

    print("######## 03 Define transloc... ########")
    cmd = "define_transloc_v5.1_mpf.py {} {}".format(basename, cutsite)
    print(cmd)
    os.system(cmd)

    print("######## 04 Define indels... ########")
    cmd = "define_indel_v5.1_mpf.py {} {} {}".format(basename, cutsite, primer_strand)
    print(cmd)
    os.system(cmd)

    print("######## 05 DEDUP... ########")
    cmd = "dedup_v5.1_mpf.py {}".format(basename)
    print(cmd)
    os.system(cmd)

    print("######## 06 substitutions... ########")
    cmd = "define_substitution.py {} {} {} {}".format(basename,cutsite,10,primer_strand)
    print(cmd)
    os.system(cmd)

    print("######## 07 filter and Statistics... ########")
    cmd = "define_statistics_add_filter.py {} {} {} {} {} {} {} CCACGCGTGCTCTACA".format(basename,genome,cutsite,500000,primer,primer_chr,primer_strand)
    print(cmd)
    os.system(cmd)

    # print("######## revise microhomology... ########")
    # cmd = "revise_microhomolog.py {} {} {}".format(basename,cutsite,primer_strand)
    # print(cmd)
    # os.system(cmd)
    
    # print("######## 08 processing translocation filter... ########")
    # cmd = "DSB_filter_update.py {} {} {} {} {} {}".format(basename,primer_chr,cutsite,primer_strand,"3fe","_Translocations")
    # print(cmd)
    # os.system(cmd)
    
    # print("######## 09 multiple mapping filter... ########")
    # cmd = "multiple_mapping_filter.py {} ".format(basename)
    # print(cmd)
    # os.system(cmd)
    
    print("All done, please check log files!")
    
    print("PEM-Q Done in {}s".format(round(time()-start_time, 3)))
    
def main():
    args = docopt(__doc__,version='PEM-Q v5.1s')
    
    kwargs = {'sample':args['<sample>'], 'cutsite':args['<cutsite>'],'genome':args['<genome>'],\
    'primer':args['<primer>'],'primer_chr':args['<primer_chr>'],'primer_start':args['<primer_start>'],\
    'primer_end':args['<primer_end>'],'primer_strand':args['<primer_strand>'],\
    'umi_length':args['<umi_length>']
    }
    
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] sample: ' + str(kwargs['sample']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] primer: ' + str(kwargs['primer']))
    print('[PEM-Q] primer_chr: ' + str(kwargs['primer_chr']))
    print('[PEM-Q] primer_start: ' + str(kwargs['primer_start']))
    print('[PEM-Q] primer_end: ' + str(kwargs['primer_end']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    print('[PEM-Q] umi_length: ' + str(kwargs['umi_length']))
    
    try:
        run_script(**kwargs)
    except KeyboardInterrupt:
        sys.stderr.write("See you~ :)\n")
        sys.exit(0)
    
if __name__ == "__main__":
    main()
    
    
