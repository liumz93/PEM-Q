#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.9.1

"""PEM-Q

Usage:
    PEM-Q <genome> <sample> <cutsite> <primer_chr> <primer_start> <primer_end> <primer_strand> <primer>

Options:
-h --help               Show this screen.
-v --version            Show version.
-p <sequence>           primer sequence.
-a <adapter>            adapter sequence.

In this script, reads will be mapped by bwa-mem, BOTH single end or 
pair end reads are compatible. When analyzing PEM-seq data, you should 
also provide adapter and primer sequences, so that adapter alignment 
and no primer filter can be done.

Input file: fastq file / Output file: informative tab files

Author: Mengzhu LIU
Last Update:2019.9.2

"""

import os
import sys
import threading
from time import time
from docopt import docopt

def run_script(sample=None, cutsite=None, genome=None, primer=None, primer_chr=None, primer_start=None, primer_end=None, primer_strand=None):
    
    start_time = time()
    
    # basename = "LY015B_100"
    # cutsite = "36573328"
    # genome = "hg38"
    # primer = "AGGATCTCACCCGGAACAGC"
    basename = sample
    
    os.system("mkdir log")
    
    print("######## Reads alignment... ########")
    cmd = "~/Scripts/script_mz/superCas/script/align_make.py {} {}_R1.fq.gz {}_R2.fq.gz -a CCACGCGTGCTCTACA -p {} -r {} -s {} -e {} -d {}".format(genome,basename,basename,primer,primer_chr,primer_start,primer_end,primer_strand)
    print(cmd)
    os.system(cmd)
    
    print("######## Barcode dedup... ########")
    cmd = "~/Scripts/script_mz/superCas/script/rmb_dedup.py {} 17".format(basename)
    print(cmd)
    os.system(cmd)
    
    print("######## Define transloc... ########")
    cmd = "~/Scripts/script_mz/superCas/script/define_transloc.py {} {}".format(basename, cutsite)
    print(cmd)
    os.system(cmd)
    
    print("######## Define indels... ########")
    cmd = "~/Scripts/script_mz/superCas/script/define_indel.py {} {}".format(basename, cutsite)
    print(cmd)
    os.system(cmd)
    
    print("All done, please check log files!")
    
    print("PEM-Q Done in {}s".format(round(time()-start_time, 3)))
    
def main():
    args = docopt(__doc__,version='PEM-Q 2.0')
    
    kwargs = {'sample':args['<sample>'], 'cutsite':args['<cutsite>'],'genome':args['<genome>'],\
    'primer':args['<primer>'],'primer_chr':args['<primer_chr>'],'primer_start':args['<primer_start>'],\
    'primer_end':args['<primer_end>'],'primer_strand':args['<primer_strand>']}
    
    print('[PEM-Q] genome: ' + str(kwargs['genome']))
    print('[PEM-Q] sample: ' + str(kwargs['sample']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] primer: ' + str(kwargs['primer']))
    print('[PEM-Q] primer_chr: ' + str(kwargs['primer_chr']))
    print('[PEM-Q] primer_start: ' + str(kwargs['primer_start']))
    print('[PEM-Q] primer_end: ' + str(kwargs['primer_end']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    
    try:
        run_script(**kwargs)
    except KeyboardInterrupt:
        sys.stderr.write("See you again :)\n")
        sys.exit(0)
    
if __name__ == "__main__":
    main()
    
    