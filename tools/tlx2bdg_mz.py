#!/usr/bin/python

print "HELLO, Mengzhu, Remember to change genome ..."
import sys, os

if len(sys.argv) != 2:
    sys.exit('<tlx file>')

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

# for pipeline version 3.0
strandconvt = {"1": "+", "-1": "-"}
bedfile = open(bedname, 'w')
for line in open(sys.argv[1]):
    if line.startswith('Qname'): continue
    l = line.strip().split()
    # bedfile.write("%s\t%d\t%s\t.\t.\t%s\n" % (l[2], int(l[3])-1, l[3], strandconvt[l[4]]))
    bedfile.write("%s\t%d\t%s\t.\t.\t%s\n" % (l[9], round(float(l[11]))-1, int(float(l[11])), l[10]))
bedfile.close()

# sort bed file
os.system('bedtools sort -i %s > %s' % (bedname, bedsort))

# for combo
fout = open("%s" % bedgraphname, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -bga -i %s -g /home/ubuntu/Scripts/hg38.genome >> %s' % (bedsort, bedgraphname))
os.system('gzip %s' % (bedgraphname))

# for positive strand
fout = open("%s" % bedgraphname_pos, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=200,30,0\n' % bedgraphname_pos)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -strand + -bga -i %s -g /home/ubuntu/Scripts/hg38.genome >> %s' % (bedsort, bedgraphname_pos))
os.system('gzip %s' % (bedgraphname_pos))

# for the minus strand
fout = open("%s" % bedgraphname_neg_tmp, 'w')
fout.write('track type=bedGraph name="%s" visibility=full color=0,30,200\n' % bedgraphname_neg)
fout.close()
####!!!!here!!!!!####
os.system('bedtools genomecov -strand - -bga -i %s -g /home/ubuntu/Scripts/hg38.genome >> %s' % (bedsort, bedgraphname_neg_tmp))
commdline = "awk '{OFS=\"\t\"}{print $1, $2, $3, -$4}' %s > %s" % (bedgraphname_neg_tmp,bedgraphname_neg)
os.system(commdline)
os.system('gzip %s' % (bedgraphname_neg))

# delete files
os.system('rm -rf %s %s' % (bedgraphname_neg_tmp, bedname))
