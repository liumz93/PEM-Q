#!/usr/bin/python

import os
import re
from time import time
from docopt import docopt
import pandas as pd
import numpy as np

print "HELLO, Mengzhu, here we are go to make a tab file to bdg file ..."
import sys, os

if len(sys.argv) != 3:
    sys.exit('<vector tab file> <vector_fa>')

vector_fa = sys.argv[2]
vector_fa_index = vector_fa.split(".")[0]
output_prefix = sys.argv[1].replace(".tab", "")
bedname = output_prefix + 'unsort.bed'
bedsort = output_prefix + '.bed'
bedgraphname = output_prefix + '.bedgraph'
bedgraphname_pos = output_prefix + '.pos.bedgraph'
bedgraphname_neg_tmp = output_prefix + '.neg.plus.bedgraph'
bedgraphname_neg = output_prefix + '.neg.bedgraph'

# convert the windows line ending to unix line ending
os.system("perl -pe 's/\r/\n/g' %s > tmpfile" % sys.argv[1])
os.system("mv tmpfile %s" % sys.argv[1])

#generate bed file
bedfile = open(bedname, 'w')

for line in open(sys.argv[1]):
    if line.startswith('Qname'): continue
    l = line.strip().split()
    if l[3] == '+': junc = l[1]
    elif l[3] == '-': junc = l[2]
    bedfile.write("%s\t%d\t%s\t.\t.\t%s\n" % ("vector", int(junc)-1, int(junc),l[3]))

bedfile.close()

# sort bed file
os.system('bedtools sort -i %s > %s' % (bedname, bedsort))

# for combo
fout = open("%s" % bedgraphname, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -bga -i %s -g %s.genome >> %s' % (bedsort, vector_fa_index, bedgraphname))
os.system('gzip %s' % (bedgraphname))

# for positive strand
fout = open("%s" % bedgraphname_pos, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=200,30,0\n' % bedgraphname_pos)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -strand + -bga -i %s -g %s.genome >> %s' % (bedsort, vector_fa_index, bedgraphname_pos))
os.system('gzip %s' % (bedgraphname_pos))

# for the minus strand
fout = open("%s" % bedgraphname_neg_tmp, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname_neg)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -strand - -bga -i %s -g %s.genome >> %s' % (bedsort, vector_fa_index, bedgraphname_neg_tmp))
commdline = "awk '{OFS=\"\t\"}{print $1, $2, $3, -$4}' %s > %s" % (bedgraphname_neg_tmp,bedgraphname_neg)
os.system(commdline)
os.system('gzip %s' % (bedgraphname_neg))

# delete files
os.system('rm -rf %s %s' % (bedgraphname_neg_tmp, bedname))
