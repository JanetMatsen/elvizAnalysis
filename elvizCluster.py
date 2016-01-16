#!/usr/bin/env python
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

MIN_ROWS = 20
DENSITY = 0.4
MIN_SAMPLES = 10

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
clusterColumns = ['Average fold', 'Reference GC']

def plotClusters(pdf, originalX, title):
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = 'k'

        class_member_mask = (labels == k)

        xy = originalX[class_member_mask & core_samples_mask]
        x = xy[clusterColumns[0]]
        y = xy[clusterColumns[1]]
        plt.plot(x, y, 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

        xy = originalX[class_member_mask & ~core_samples_mask]
        plt.plot(x, y, 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=3)

        plt.xlabel(clusterColumns[0], fontsize=10)
        plt.ylabel(clusterColumns[1], fontsize=10)

    plt.title(title)
    plt.axes([originalX[clusterColumns[0]].min(), originalX[clusterColumns[0]].max(),
            originalX[clusterColumns[1]].min(), originalX[clusterColumns[1]].max()])
    pdf.savefig()
    plt.close()

# create a multiage PDF for storing the plots
with PdfPages('elvizCluster.pdf') as pdf:
    # find unique values of taxonomy columns
    dfgb = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])
    for key in dfgb.indices.keys():
        idx = dfgb.indices[key]
        taxRows = df.iloc[idx]
        if len(taxRows) < MIN_ROWS:
            continue
        # normalize all dimensions to be used in clustering, e.g. GC, coverage, rpk
        taxRowsClusterColumns = StandardScaler().fit_transform(taxRows[clusterColumns])

        db = DBSCAN(eps=DENSITY, min_samples=MIN_SAMPLES).fit(taxRowsClusterColumns)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_

        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        if n_clusters_ < 1:
            continue

        title = ', '.join(key)
        print(title)
        print('Estimated number of clusters: %d' % n_clusters_)
        plotClusters(pdf, taxRows[clusterColumns], title)
