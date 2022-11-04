#!/usr/bin/python

print "HELLO, Stranger, this scritp is for generating bdg file of PEM-Q results file ..."
import sys, os
import pandas as pd
import numpy as np
import subprocess

if len(sys.argv) != 3:
    sys.exit('<PEMQ tab file> <genome>')

output_prefix = sys.argv[1].replace(".tab", "")
genome = sys.argv[2]
bedname = output_prefix + 'unsort.bed'
bedsort = output_prefix + '.bed'
bedgraphname = output_prefix + '.bedgraph'
bedgraphname_pos = output_prefix + '.pos.bedgraph'
bedgraphname_neg_tmp = output_prefix + '.neg.plus.bedgraph'
bedgraphname_neg = output_prefix + '.neg.bedgraph'

# convert the windows line ending to unix line ending
os.system("perl -pe 's/\r/\n/g' %s > tmpfile" % sys.argv[1])
os.system("mv tmpfile %s" % sys.argv[1])

# for PEM-Q results
data = pd.read_csv(sys.argv[1], sep = '\t', index_col=False, low_memory=False)
data['1'] = 1
data['.']="."
data['Junction-1'] = data['Junction'] - data['1']
print data['Junction'].head()
print data['Junction-1'].head()
data.to_csv(bedname, header = True, sep = '\t', index=False, columns = [u'Rname',u'Junction-1',u'Junction',u'.',u'.',u'Strand'])
os.system("sed -i '1d' {}".format(bedname))

# sort bed file
os.system('bedtools sort -i %s > %s' % (bedname, bedsort))

# for combo
fout = open("%s" % bedgraphname, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -bga -i %s -g /home/ubuntu/Scripts/%s.genome >> %s' % (bedsort, genome, bedgraphname))
os.system('gzip %s' % (bedgraphname))

# for positive strand
fout = open("%s" % bedgraphname_pos, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=200,30,0\n' % bedgraphname_pos)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -strand + -bga -i %s -g /home/ubuntu/Scripts/%s.genome >> %s' % (bedsort, genome, bedgraphname_pos))
os.system('gzip %s' % (bedgraphname_pos))

# for the minus strand
fout = open("%s" % bedgraphname_neg_tmp, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname_neg)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -strand - -bga -i %s -g /home/ubuntu/Scripts/%s.genome >> %s' % (bedsort, genome, bedgraphname_neg_tmp))
commdline = "awk '{OFS=\"\t\"}{print $1, $2, $3, -$4}' %s > %s" % (bedgraphname_neg_tmp,bedgraphname_neg)
os.system(commdline)
os.system('gzip %s' % (bedgraphname_neg))

# delete files
os.system('rm -rf %s %s' % (bedgraphname_neg_tmp, bedname))
