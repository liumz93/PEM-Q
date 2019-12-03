#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.9.25

import os
from time import time
from docopt import docopt
import pandas as pd

data  = pd.read_csv("MZ003_transloc_all.tab",sep = "\t")
datb = data[data.duplicated(['Sequence'],keep = False)]
datb = datb.sort_values(['Prey_rname','Prey_strand','Prey_start','Prey_end','Rname','Strand','Junction'])
datb.to_csv("dup.tab", header = True, sep = '\t', index=False)