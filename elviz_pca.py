import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('TkAgg')
from elviz_utils import read_sample_info



def import_elviz_data():
    """
    Read pandas data from summarised table
    :return: pandas dataframe
    """
    return pd.read_csv("./results/reduced_data--genus_only.csv",
                       dtype={'abundance': 'float'})


def pivot_for_pca(dataframe):
    """
    pivot to the format required for scikitlearn pca
    :param dataframe:input dataframe reduced to genus only
    :return:dataframe with columns as samples and rows as genera
    """
    dataframe = dataframe.pivot(index='Genus', columns='ID',
                                values='abundance')
    # fill NA values with 0.
    return dataframe.fillna(0)


def most_abundant_genera_for_pca(data, top_percent):
    data['row_sum'] = data.sum(axis=1)
    data.head()
    # reduce to top __ %
    # find the row_sum corresponding to the top 10%.
    num_rows_to_keep = int(round(data.shape[0]*top_percent/100.))
    print(num_rows_to_keep)

    # sort by row_sum so I can take the firs number of rows.
    data.sort(columns='row_sum', ascending=False, inplace=True)
    # print data.head()

    # data.groupby('row_sum').head(num_rows_to_keep)
    del data['row_sum']
    return data.head(num_rows_to_keep)


def colnames_to_sample_info_array(dataframe):
    col_df = pd.DataFrame({'ID' : dataframe.columns.values.tolist()})
    sample_info = read_sample_info()
    return pd.merge(col_df, sample_info, how='left')
