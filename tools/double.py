#!/bin/python3
import argparse
import os
#def argas
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",type=str,required=True,help="input insert_transloc.tab")
    parser.add_argument("-o", type=str, required=True, help="output result file")
    parser.add_argument("-chr", type=str, required=True, help="cutting chromosome")
    parser.add_argument("-lstart", type=int, required=True, help="count locus start minus 1")
    parser.add_argument("-lend", type=int, required=True, help="count locus end")
    parser.add_argument("-genome", type=str, default="/home/ubuntu/genomes/hg38/hg38.fa" ,help="cutting genome, default hg38")
    parser.add_argument("-strand",type=str,required=True,help="cutting strand:+ strand = 1 , - strand = 0")
    argas=parser.parse_args()
    return argas
argas=parse_args()

#input file;output file
path=argas.i
opath=argas.o
of=open('tmp.csv',"w")
of.write("locus,nt,num"+"\n")
#start locus
zero=argas.lstart
length=argas.lend - argas.lstart + 1

#get genome seq
os.system("samtools faidx %s %s:%d-%d > seq.fa"%(argas.genome,argas.chr,argas.lstart,argas.lend))
for line in open("seq.fa","r"):
    seq=""
    if not line.startswith(">"):
        seq=line.replace('\n','').strip()
os.system("rm -rf seq.fa")

#strand
strand=argas.strand

def count_num(nt,max,min):
    nt_num = 0
    if strand == "0":
        for line in open(path,"r"):
            l=line.split('\t')
            if l[3]==str(max) and l[11]==str(min) and l[13]== str(nt):
                nt_num += 1
    if strand == "1":
        for line in open(path,"r"):
            l=line.split('\t')
            if l[4]==str(min) and l[11]==str(max) and l[13]== str(nt):
                nt_num += 1
    return str(nt_num)


num_range=[i for i in range(0,length)]
for num in num_range:
    before_num = int(num) - 1
    a = str(zero + int(num))
    b = str(zero + int(num) - 1)
    c = str(zero + int(num) - 2)
    ntstr=list("ATGC")
    if before_num < 0: continue
    for nt in ntstr:
        if seq[int(before_num)]==nt:
            of.write(a+","+nt+","+count_num(nt,b,c)+"\n")
        else:
            of.write(a+","+nt+","+count_num(nt,a,b)+"\n")
of.close()

import pandas as pd
import numpy as np
data=pd.read_csv('tmp.csv')
creat_data=pd.pivot_table(data,values=[u'num'],index=[u'nt'],columns=[u'locus'],aggfunc=[np.sum])
creat_data.to_csv(opath)
os.system("rm -rf tmp.csv")
