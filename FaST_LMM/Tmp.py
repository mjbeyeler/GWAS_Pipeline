# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:25:33 2017

@author: micha
"""

import pandas as pd
df = pd.read_csv('../Data/freeze2.common.rel.mat', sep='\t').astype(str)
lines = pd.Series(df.columns.values)
lines = lines.str.extract('(.* )')
lines = lines.values.tolist()
import re
indices = re.compile()
