#!/usr/bin/env python3
"""
"""

import getopt
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import re
import seaborn as sns
import sys

from matplotlib.backends.backend_pdf import PdfPages
import scipy.spatial.distance
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn import mixture
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

from elviz_utils import IMPORT_DATA_TYPES, read_elviz_CSV, read_elviz_CSVs, read_pickle_or_CSVs

MIN_ROWS = 20
EPS = .5
MIN_SAMPLES = 4
MAX_AVG_FOLD = 500  # I've seen over 20k, 500 seems to be 2x typical
DATA_PICKLE = 'data/data.pkl'  # filename of previously parsed data
RAW_DATA_DIR = './raw_data/'  # location of CSV files
RESULTS_DIR = './results/'  # location of results output
HEURISTIC_SAMPlE_SIZE = 10
HEURISTIC_PDF = 'heuristic.pdf'
HEURISTIC_PICKLE = 'data/heuristic.pkl'
DPGMM_N_COMPONENTS = 3

CLUSTER_COLUMNS = ['Average fold', 'Reference GC']

sns.set()
sns.set(style="whitegrid")


def dbscan_heuristic(elviz_data, scaler):
    if os.path.isfile(HEURISTIC_PICKLE):
        print("reading %s for previously computed heuristic data" % HEURISTIC_PICKLE)
        with open(HEURISTIC_PICKLE, 'rb') as file:
            distances = pickle.load(file)
    else:
        print("processing full data set to compute heuristic for epsilon / N")
        distances = []
        for filename in elviz_data.keys():
            # progress monitor
            print(".", end="")
            sys.stdout.flush()
            df = elviz_data[filename]
            dfgb = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'])

            for key in dfgb.indices.keys():
                idx = dfgb.indices[key]
                tax_rows = df.iloc[idx]
                if len(tax_rows) < HEURISTIC_SAMPlE_SIZE:
                    continue
                # select out the columns of relevance
                #reduced_df = scaler.transform(df[CLUSTER_COLUMNS])
                reduced_df = pd.DataFrame(scaler.transform(df[CLUSTER_COLUMNS]), columns=CLUSTER_COLUMNS)
                # create a random sample of the rows
                random_sample = reduced_df.sample(HEURISTIC_SAMPlE_SIZE)
                # create array of complex #s using a as first columns and b from second (a + bi)
                p1 = (random_sample[CLUSTER_COLUMNS[0]] + 1j * random_sample[CLUSTER_COLUMNS[1]]).values
                p2 = (reduced_df[CLUSTER_COLUMNS[0]] + 1j * reduced_df[CLUSTER_COLUMNS[1]]).values

                # calculate all the distances, between each point in random sample
                # and all data (using an array-broadcasting trick)
                all_dists = abs(p1[..., np.newaxis] - p2)
                # sort along each row
                all_dists.sort(axis=1)
                # if we don't tolist this, it will save the entire matrix
                # with a view representation stored in distances, converting
                # this to a list sidesteps this and the associated memory
                # problem
                distances.append(all_dists[:, MIN_SAMPLES].tolist())
        with open(HEURISTIC_PICKLE, 'wb') as file:
            pickle.dump(distances, file, pickle.HIGHEST_PROTOCOL)

    distances = [item for sublist in distances for item in sublist]
    with PdfPages(RESULTS_DIR + HEURISTIC_PDF) as pdf:
        print("\ncomputing histogram and making figure")
        # plt.hist(distances, bins=50, range=(0, 2))
        distances.sort(reverse=True)
        plt.plot(distances[0:100])
        plt.title("THIS IS A PLOT TITLE, YOU BET")
        plt.xlabel("Sorted index")
        plt.ylabel("k-dist (k = %d)" % MIN_SAMPLES)
        pdf.savefig()
        plt.close()


def plot_clusters(pdf, df, title, labels, core_samples_mask, limits):
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)), alpha=0.6)
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # black used for noise.
            col = [0, 0, 0, .6]

        class_member_mask = (labels == k)

        # plot the core samples
        xy = df[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=6)

        # plot those joined by extension
        xy = df[class_member_mask & ~core_samples_mask]
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


def main(argv):
    dbscan_heuristic_mode = False
    dpgmm_mode = False
    try:
        opts, args = getopt.getopt(argv,"he")
    except getopt.GetoptError:
        print('elviz_cluster.py [-h] [-e] [-g]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('elviz_cluster.py [-h] [-e]')
            print('  -h = help, -e = run dbscan' +
                  ' epsilon heuristic plot generation code')
            print('  -g = use a DPGMM for clustering')
            sys.exit()
        elif opt == '-e':
            dbscan_heuristic_mode = True
        elif opt == '-g':
            dpgmm_mode = True

    [elviz_data, combined_df] = read_pickle_or_CSVs(DATA_PICKLE, RAW_DATA_DIR)

    # Setup plotting limits
    print("determining plotting limits")
    limits = {"x": [combined_df['Average fold'].min(), MAX_AVG_FOLD],
              "y": [combined_df['Reference GC'].min(), combined_df['Reference GC'].max()]}
    # Below changed in favor of fixed MAX
    # limits["x"] = [combined_df['Average fold'].min(), combined_df['Average fold'].max()]
    # fixed MAX below

    print("normalizing data prior to clustering")
    # normalize the combined data to retrieve the normalization parameters
    scaler = StandardScaler().fit(combined_df[CLUSTER_COLUMNS])
    # serializing outputs

    if dbscan_heuristic_mode:
        print("making DBSCAN heuristic plots")
        dbscan_heuristic(elviz_data, scaler)
        os.sys.exit()

    print("serially processing files")
    for filename in elviz_data.keys():
        pdf_filename = filename.replace("csv", "pdf")
        # skip if the PDF already exists
        if os.path.isfile(RESULTS_DIR + pdf_filename):
            print("skiping file %s" % filename)
            continue
        print("processing file %s" % filename)

        df = elviz_data[filename]

        # create a multipage PDF for storing the plots
        with PdfPages(RESULTS_DIR + pdf_filename) as pdf:
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

                if not dpgmm_mode:
                    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES)
                else:
                    db = mixture.GPGMM(n_components=DPGMM_N_COMPONENTS,
                                       covariance_type='full')
                db.fit(tax_rows_cluster_columns)
                core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                core_samples_mask[db.core_sample_indices_] = True
                labels = db.labels_

                # number of clusters in labels, ignoring noise if present.
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

                if n_clusters_ < 1:
                    continue

                title = ', '.join(key)
                plot_clusters(pdf, scaler.inverse_transform(tax_rows_cluster_columns),
                              title, labels, core_samples_mask, limits)


if __name__ == "__main__":
    main(sys.argv[1:])
