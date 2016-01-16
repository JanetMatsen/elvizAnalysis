#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="whitegrid")

def read_elviz_csv(filename):
    df = pd.read_csv(filename, sep=",", low_memory=False)
    # repalce nans with ""
    df.fillna("", inplace=True)
    return df

# Read in all the .csv files and concatenate them.
df = read_elviz_csv("Metagenome_9_HOW4.csv");
#print(df.head())

# add some useful columns
df['reads'] = df['Plus reads'] + df['Minus reads']
df['rpk'] = df['reads']/(df['Length']/1000)

# find unique values of taxonomy columns
dfgb = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])
for key in dfgb.indices.keys():
    idx = dfgb.indices[key]
    taxRows = df.iloc[idx]
    print(taxRows.head())
