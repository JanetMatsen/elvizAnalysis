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
DENSITY = 0.3
MIN_SAMPLES = 10

clusterColumns = ['Average fold', 'Reference GC']

s = StandardScaler()
sns.set(style="whitegrid")

# read in metadata about samples
meta_info = pd.read_csv("./data/sample_meta_info.tsv", sep="\t")

def readElvizCSV(filename):
    df = pd.read_csv(filename, sep=",", low_memory=False)
    # repalce nans with ""
    df.fillna("", inplace=True)
    return df

def plotClusters(pdf, X, title):
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)), alpha=0.6)
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # black used for noise.
            col = [0, 0, 0, .6]

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=3)

    plt.title(title)
    plt.xlabel(clusterColumns[0], fontsize=10)
    plt.ylabel(clusterColumns[1], fontsize=10)
    pdf.savefig()
    plt.close()


# read in all the .csv files and process them
elviz_files = [f for f in os.listdir("./data/") if ".csv" in f]
for f in elviz_files:
    fpdf = f.replace("csv", "pdf")

    # skip if the PDF already exists
    if os.path.isfile("./results/" + fpdf):
        continue

    print("processing file %s" % f)
    # read the dataframe from the csv
    df = readElvizCSV("./data/" + f)

    # create a multiage PDF for storing the plots
    with PdfPages('./results/' + fpdf) as pdf:
        # find unique values of taxonomy columns
        dfgb = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])
        for key in dfgb.indices.keys():
            idx = dfgb.indices[key]
            taxRows = df.iloc[idx]
            if len(taxRows) < MIN_ROWS:
                continue
            # normalize all dimensions to be used in clustering, e.g. GC, coverage, rpk
            taxRowsClusterColumns = s.fit_transform(taxRows[clusterColumns])

            db = DBSCAN(eps=DENSITY, min_samples=MIN_SAMPLES).fit(taxRowsClusterColumns)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_

            # number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

            if n_clusters_ < 1:
                continue

            title = ', '.join(key)
            #print(title)
            #print('Estimated number of clusters: %d' % n_clusters_)
            plotClusters(pdf, s.inverse_transform(taxRowsClusterColumns), title)

