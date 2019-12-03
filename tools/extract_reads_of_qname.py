#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.23



import sys
if len(sys.argv) != 4:
    sys.exit('<bam_file> <qname_list> <name>')
import os
import pysam
from time import time
from docopt import docopt
import pandas as pd

bam_file_name = sys.argv[1]
qname_list_name = sys.argv[2]
basename = sys.argv[3]

qname_list_file = pd.read_csv(qname_list_name,names = ["Qname"])
qname_list = qname_list_file["Qname"]
bam = pysam.AlignmentFile(bam_file_name,'rb')
extract_bam = pysam.AlignmentFile(basename+"_extrac.bam", "wb", template=bam)

bam_indexed = pysam.IndexedReads(bam)
bam_indexed.build()
n = 0
for name in qname_list:
    try:
        bam_indexed.find(name)
    except KeyError:
        pass
    else:

        iterator = bam_indexed.find(name)
        for x in iterator:
            n = n + 1
            extract_bam.write(x)
bam.close()
extract_bam.close()
pysam.sort("-o", basename+"_extrac.sort.bam", basename+"_extrac.bam")
cmd = "samtools index {}".format(basename+"_extrac.sort.bam")
os.system(cmd)