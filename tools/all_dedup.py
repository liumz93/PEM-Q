#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.6.23
"""all_dedup

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
    all_dedup   <basename>

Options:
-h --help               Show this screen.
-v --version            Show version.

This script find bait break site and translocation site from bwa align file.
And also filter reads for indel define.

Input file: transloc and indel tab file / Output file: unique transloc and indel tab file

Author: Mengzhu LIU
Last Update:2019.6.24

"""
import os
from time import time
from docopt import docopt
import pandas as pd

def dedup_all_transloc(basename):
    transloc = pd.read_csv(basename + "_transloc_all.tab",sep = '\t')
    data = transloc.drop_duplicates(['Bait_start', 'Bait_end', 'Prey_rname', 'Prey_strand', 'Junction'],keep='first')
    data.to_csv(basename + "_transloc_all_unique.tab", header = True, sep = '\t', index=False)

def main():
    start_time = time()
    
    args = docopt(__doc__,version='all_dedup 1.0')
    
    kwargs = {'basename':args['<basename>']}
    print('[superCasQ] basename: ' + str(kwargs['basename']))

    dedup_all_transloc(str(kwargs['basename']))
    print("\nall_dedup.py Done in {}s".format(round(time()-start_time, 3)))
    
if __name__ == '__main__':
    main()
    