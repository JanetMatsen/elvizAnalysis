#!/usr/bin/env python
import sys
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

from elviz_utils import IMPORT_DATA_TYPES, read_elviz_CSV, read_elviz_CSVs, read_pickle_or_CSVs

MIN_ROWS = 20
EPS = 0.15
MIN_SAMPLES = MIN_ROWS
MAX_AVG_FOLD = 500 # I've seen over 20k, 500 seems to be 2x typical
DATA_PICKLE = 'data.pkl' # filename of previously parsed data
DATA_DIR = './results'   # location of CSV files

CLUSTER_COLUMNS = ['Log10 Average fold', 'Reference GC']

sns.set(style="whitegrid")


def plot_clusters(pdf, X, title, labels, core_samples_mask, limits):
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)), alpha=0.6)
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # black used for noise.
            col = [0, 0, 0, .6]

        class_member_mask = (labels == k)

        # plot the core samples
        xy = X[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

        # plot those joined by extension
        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=3)

    # format the plot
    plt.title(title)
    plt.xlabel(CLUSTER_COLUMNS[0], fontsize=10)
    plt.ylabel(CLUSTER_COLUMNS[1], fontsize=10)
    plt.xlim(limits["x"])
    plt.ylim(limits["y"])
    pdf.savefig()
    plt.close()


def main():
    [ elviz_data, combined_df ] = read_pickle_or_CSVs(DATA_PICKLE, DATA_DIR)

    # setup plotting limits
    print("determining plotting limits")
    limits = { }
    # below changed in favor of fixed MAX
    # limits["x"] = [combined_df['Log10 Average fold'].min(), combined_df['Log10 Average fold'].max()]
    # fixed MAX below
    limits["x"] = [combined_df['Log10 Average fold'].min(), MAX_AVG_FOLD]
    limits["y"] = [combined_df['Reference GC'].min(), combined_df['Reference GC'].max()]

    print("normalizing data prior to clustering")
    # normalize the combined data to retrieve the normalization parameters
    scaler = StandardScaler().fit(combined_df[CLUSTER_COLUMNS])
    # serializing outputs

    print("serially processing files")
    for filename in elviz_data.keys():
        pdf_filename = filename.replace("csv", "pdf")
        # skip if the PDF already exists
        if os.path.isfile("./results/" + pdf_filename):
            print("skiping file %s" % filename)
            continue
        print("processing file %s" % filename)

        df = elviz_data[filename]

        # create a multipage PDF for storing the plots
        with PdfPages('./results/' + pdf_filename) as pdf:
            # find unique values of taxonomy columns
            dfgb = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])
            for key in dfgb.indices.keys():
                idx = dfgb.indices[key]
                tax_rows = df.iloc[idx]
                if len(tax_rows) < MIN_ROWS:
                    continue
                # normalize all dimensions to be used in clustering, e.g. GC, coverage, rpk
                # reuse the scaler we created from all of the data for the transform
                tax_rows_cluster_columns = scaler.transform(tax_rows[CLUSTER_COLUMNS])

                db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES).fit(tax_rows_cluster_columns)
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
                plot_clusters(pdf, scaler.inverse_transform(tax_rows_cluster_columns),
                        title, labels, core_samples_mask, limits)


if __name__ == "__main__":
    main()

