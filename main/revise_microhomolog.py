#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2021.01.14
"""define_microhomology

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
    revise_microhomology    <basename> <cutsite> <primer_strand>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script .


Author: Mengzhu LIU
Last Update:2021.01.14

"""


import os
from time import time
from docopt import docopt
import pandas as pd
import numpy as np
    
def revise_microhomology(basename=None,cutsite=None,primer_strand=None):
    
    
    print("Processing microhomology statistics...")
    os.system('mkdir results')
    
    Deletion_file = "results/" + basename+"_Deletion.tab"
    Deletion = pd.read_csv(Deletion_file, sep = '\t', index_col=False, low_memory=False)
    Deletion['cutsite_left'] = int(cutsite)
    Deletion['cutsite_right'] = int(cutsite) + 1
    if primer_strand == "+":
        Deletion.loc[:,'Off_cutsite_right'] = np.where(Deletion['Prey_start'] - Deletion['cutsite_right'] < 0, Deletion['Prey_start'] - Deletion['cutsite_right'], 0)
        Deletion.loc[:,'Off_cutsite_left'] = np.where(Deletion['cutsite_left'] - Deletion['Bait_end'] < 0, Deletion['cutsite_left'] - Deletion['Bait_end'], 0)
        Deletion['Revised_microhomology_len'] = Deletion['Microhomolog_len'] + Deletion['Off_cutsite_right'] + Deletion['Off_cutsite_left']
    #primer_strand == "-"
    else:
        Deletion.loc[:,'Off_cutsite_right'] = np.where(Deletion['Bait_start'] - Deletion['cutsite_right'] < 0, Deletion['Bait_start'] - Deletion['cutsite_right'], 0)
        Deletion.loc[:,'Off_cutsite_left'] = np.where(Deletion['cutsite_left'] - Deletion['Prey_end'] < 0, Deletion['cutsite_left'] - Deletion['Prey_end'] < 0, 0)
        Deletion['Revised_microhomology_len'] = Deletion['Microhomolog_len'] + Deletion['Off_cutsite_right'] + Deletion['Off_cutsite_left']
    
    print(Deletion['Qname'].count())
    # print(Deletion[Deletion['Revised_microhomology_len']<0])
    Deletion = Deletion[Deletion['Revised_microhomology_len']>=0]
    print(Deletion['Qname'].count())
    Deletion.to_csv("results/" + basename+"_revise_Deletion.tab", header = True, sep = '\t', index=False)
    
    #statistics
    
    DJ_count = Deletion[Deletion['Microhomolog_len'] == 0]['Qname'].count()
    M1_count = Deletion[Deletion['Microhomolog_len'] == 1]['Qname'].count()
    M2_count = Deletion[Deletion['Microhomolog_len'] == 2]['Qname'].count()
    M3_count = Deletion[Deletion['Microhomolog_len'] == 3]['Qname'].count()
    M4_count = Deletion[Deletion['Microhomolog_len'] == 4]['Qname'].count()
    M5_count = Deletion[Deletion['Microhomolog_len'] == 5]['Qname'].count()
    M6_count = Deletion[Deletion['Microhomolog_len'] == 6]['Qname'].count()
    M7_count = Deletion[Deletion['Microhomolog_len'] == 7]['Qname'].count()
    M8_count = Deletion[Deletion['Microhomolog_len'] >= 8]['Qname'].count()
    
    del_microhomology_stats_file = open("results/" + basename+"_del_microhomology_statistics.txt","w")
    del_microhomology_stats_file.write("D.J."+"\t"+str(DJ_count)+"\n")
    del_microhomology_stats_file.write("1"+"\t"+str(M1_count)+"\n")
    del_microhomology_stats_file.write("2"+"\t"+str(M2_count)+"\n")
    del_microhomology_stats_file.write("3"+"\t"+str(M3_count)+"\n")
    del_microhomology_stats_file.write("4"+"\t"+str(M4_count)+"\n")
    del_microhomology_stats_file.write("5"+"\t"+str(M5_count)+"\n")
    del_microhomology_stats_file.write("6"+"\t"+str(M6_count)+"\n")
    del_microhomology_stats_file.write("7"+"\t"+str(M7_count)+"\n")
    del_microhomology_stats_file.write(">=8"+"\t"+str(M8_count)+"\n")
    del_microhomology_stats_file.close()
    
    small_del = Deletion[Deletion['deletion_length']<=100]
    
    DJ_count = small_del[small_del['Microhomolog_len'] == 0]['Qname'].count()
    M1_count = small_del[small_del['Microhomolog_len'] == 1]['Qname'].count()
    M2_count = small_del[small_del['Microhomolog_len'] == 2]['Qname'].count()
    M3_count = small_del[small_del['Microhomolog_len'] == 3]['Qname'].count()
    M4_count = small_del[small_del['Microhomolog_len'] == 4]['Qname'].count()
    M5_count = small_del[small_del['Microhomolog_len'] == 5]['Qname'].count()
    M6_count = small_del[small_del['Microhomolog_len'] == 6]['Qname'].count()
    M7_count = small_del[small_del['Microhomolog_len'] == 7]['Qname'].count()
    M8_count = small_del[small_del['Microhomolog_len'] >= 8]['Qname'].count()
    
    small_del_microhomology_stats_file = open("results/" + basename+"_small_del_microhomology_statistics.txt","w")
    small_del_microhomology_stats_file.write("D.J."+"\t"+str(DJ_count)+"\n")
    small_del_microhomology_stats_file.write("1"+"\t"+str(M1_count)+"\n")
    small_del_microhomology_stats_file.write("2"+"\t"+str(M2_count)+"\n")
    small_del_microhomology_stats_file.write("3"+"\t"+str(M3_count)+"\n")
    small_del_microhomology_stats_file.write("4"+"\t"+str(M4_count)+"\n")
    small_del_microhomology_stats_file.write("5"+"\t"+str(M5_count)+"\n")
    small_del_microhomology_stats_file.write("6"+"\t"+str(M6_count)+"\n")
    small_del_microhomology_stats_file.write("7"+"\t"+str(M7_count)+"\n")
    small_del_microhomology_stats_file.write(">=8"+"\t"+str(M8_count)+"\n")
    small_del_microhomology_stats_file.close()
    
    large_del = Deletion[Deletion['deletion_length']>100]
    
    DJ_count = large_del[large_del['Microhomolog_len'] == 0]['Qname'].count()
    M1_count = large_del[large_del['Microhomolog_len'] == 1]['Qname'].count()
    M2_count = large_del[large_del['Microhomolog_len'] == 2]['Qname'].count()
    M3_count = large_del[large_del['Microhomolog_len'] == 3]['Qname'].count()
    M4_count = large_del[large_del['Microhomolog_len'] == 4]['Qname'].count()
    M5_count = large_del[large_del['Microhomolog_len'] == 5]['Qname'].count()
    M6_count = large_del[large_del['Microhomolog_len'] == 6]['Qname'].count()
    M7_count = large_del[large_del['Microhomolog_len'] == 7]['Qname'].count()
    M8_count = large_del[large_del['Microhomolog_len'] >= 8]['Qname'].count()
    
    large_del_microhomology_stats_file = open("results/" + basename+"_large_del_microhomology_statistics.txt","w")
    large_del_microhomology_stats_file.write("D.J."+"\t"+str(DJ_count)+"\n")
    large_del_microhomology_stats_file.write("1"+"\t"+str(M1_count)+"\n")
    large_del_microhomology_stats_file.write("2"+"\t"+str(M2_count)+"\n")
    large_del_microhomology_stats_file.write("3"+"\t"+str(M3_count)+"\n")
    large_del_microhomology_stats_file.write("4"+"\t"+str(M4_count)+"\n")
    large_del_microhomology_stats_file.write("5"+"\t"+str(M5_count)+"\n")
    large_del_microhomology_stats_file.write("6"+"\t"+str(M6_count)+"\n")
    large_del_microhomology_stats_file.write("7"+"\t"+str(M7_count)+"\n")
    large_del_microhomology_stats_file.write(">=8"+"\t"+str(M8_count)+"\n")
    large_del_microhomology_stats_file.close()
    
def main():

    start_time = time()

    args = docopt(__doc__,version='revise_microhomology 1.0')

    kwargs = {'basename':args['<basename>'],'cutsite':args['<cutsite>'],'primer_strand':args['<primer_strand>']}
    print('[PEM-Q] basename: ' + str(kwargs['basename']))
    print('[PEM-Q] cutsite: ' + str(kwargs['cutsite']))
    print('[PEM-Q] primer_strand: ' + str(kwargs['primer_strand']))
    
    
    revise_microhomology(**kwargs)

    print("\nrevise_microhomology.py Done in {}s".format(round(time()-start_time, 3)))

if __name__ == '__main__':
    main()
