#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Mengzhu
#Date:2019.8.16

import pandas as pd
import numpy as np
import sys, os
from time import time

data = pd.read_table("YJ052a_transloc_all.tab",sep = '\t')
datb = pd.read_table("YJ052a_vector.tab",sep = '\t')
vector_transloc = data[data['Qname'].isin(datb['Qname'])]
vector_transloc.to_csv("YJ052a_transloc_vector_all.tab",header=True,sep='\t',index=False)