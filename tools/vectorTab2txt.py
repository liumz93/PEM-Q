#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.10.30

print("HELLO, Mengzhu, Remember to change vector genome file ...")
import os
import sys
from time import time
from docopt import docopt
import pandas as pd
import numpy as np

if len(sys.argv) != 3:
    sys.exit('<vector tab file> <binsize>')

output_prefix = sys.argv[1].replace(".tab", "")
binsize = int(sys.argv[2])

data = pd.read_csv(sys.argv[1], sep = "\t")
data.loc[:,'Junction'] = np.where(data['Vector_strand'] == "+", data.Vector_start, data.Vector_end)
nrow = data['Qname'].count()
# data['Rname'] = ['vector']*nrow
# data['Strand'] = data['Vector_strand']
genome = open("/home/mengzhu/database/px330/vector.genome")
for line in genome:
    l = line.strip().split()
    vector_len = int(l[1])
genome.close()
n = vector_len//binsize
m = vector_len%binsize

txtfile = open(output_prefix+".txt", "w")
txtPOSfile = open(output_prefix+"_pos.txt", "w")
txtNEGfile = open(output_prefix+"_neg.txt", "w")
i = 0
Rstart = 1
Rend = m
while i <= n:
    condition_pos = (data['Junction']>=Rstart) & (data['Junction']<Rend) & (data['Vector_strand'] == '+')
    condition_neg = (data['Junction']>=Rstart) & (data['Junction']<Rend) & (data['Vector_strand'] == '-')
    condition_total = (data['Junction']>=Rstart) & (data['Junction']<Rend)
    data_total = data[condition_total]
    Hits = data_total['Qname'].count()
    txtfile.write("vector"+"\t"+
                    str(Rstart)+"\t"+
                    str(Rend)+"\t"+
                    str(Hits)+"\n")
    data_pos = data[condition_pos]
    Hits = data_pos['Qname'].count()
    txtPOSfile.write("vector"+"\t"+
                    str(Rstart)+"\t"+
                    str(Rend)+"\t"+
                    str(Hits)+"\n")
    data_neg = data[condition_neg]
    Hits = data_neg['Qname'].count()
    txtNEGfile.write("vector"+"\t"+
                    str(Rstart)+"\t"+
                    str(Rend)+"\t"+
                    str(Hits)+"\n")
    Rstart = Rend
    Rend = Rstart + binsize
    i+=1
txtfile.close()
txtPOSfile.close()
txtNEGfile.close()
    
# data.to_csv("unique/" + output_prefix + "_plot.tab", header = True, sep = '\t', index=False)
