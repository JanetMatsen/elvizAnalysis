#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

minRows = 10

sns.set(style="whitegrid")

def readElvizCSV(filename):
    df = pd.read_csv(filename, sep=",", low_memory=False)
    # repalce nans with ""
    df.fillna("", inplace=True)
    return df

# Read in all the .csv files and concatenate them.
df = readElvizCSV("Metagenome_9_HOW4.csv");
#print(df.head())

# add some useful columns
df['reads'] = df['Plus reads'] + df['Minus reads']
df['rpk'] = df['reads']/(df['Length']/1000)

clusterColumns = ['Reference GC', 'rpk', 'Average fold']

# find unique values of taxonomy columns
dfgb = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])
for key in dfgb.indices.keys():
    idx = dfgb.indices[key]
    taxRows = df.iloc[idx]
    if len(taxRows) < minRows:
        continue
    # normalize all dimensions to be used in clustering, e.g. GC, coverage, rpk
    # begin by zero centering the columns
    #for col in clusterColumns:
    ##    taxRows[col] -= np.mean(taxRows[col])
    ##   # normalize by the range or the variance?
    ##    # for now, variance
    ##    taxRows[col] /= np.var(taxRows[col])
