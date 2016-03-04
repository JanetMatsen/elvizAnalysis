import pandas as pd
import numpy
import os
import pickle

# Pandas is happy when you tell it the data types 
# for each column in the raw Elviz .csv files
IMPORT_DATA_TYPES = {'datasetId': 'str',
                     'contigId': 'str',
                     'Average fold': 'float',
                     'Length': 'int',
                     'Reference GC': 'float',
                     'Covered percent': 'float',
                     'Covered bases': 'int',
                     'Plus reads': 'int',
                     'Minus reads': 'int',
                     'Median fold': 'float',
                     'Read GC': 'float',
                     'Complete Lineage': 'str',
                     'IMG scaffold_oid': 'str',
                     'Kingdom': 'str',
                     'Phylum': 'str',
                     'Class': 'str',
                     'Order': 'str',
                     'Family': 'str',
                     'Genus': 'str',
                     'Species': 'str'
                     }
IMPORT_METAINFO_TYPES = {'ID': 'str',
                         'oxy': 'str',
                         'rep': 'int',
                         'week': 'int',
                         'project': 'int'}


def read_sample_info():
    '''
    Read in sample_meta_info.tsv using particular dtypes
    '''
    return pd.read_csv('./raw_data/sample_meta_info.tsv',
                       dtype=IMPORT_METAINFO_TYPES,
                       sep='\t')


def read_elviz_CSV(filename):
    df = pd.read_csv(filename, sep=",", dtype=IMPORT_DATA_TYPES)
    # repalce nans with ""
    df.fillna("", inplace=True)
    return df


def make_directory(dirpath):
    """
    Make the directory if it doesn't exist.
    """
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)


def prepare_plot_dir(dirpath):
    """
    Add a trailing / to the plot dir if missing and creat dir if needed.

    :param dirpath: path of existing or planned directory
    :return: dirpath with / as the final character
    """
    # add a slash to filename if it wasn't provided.
    if dirpath[-1] != '/':
        dirpath += '/'
    make_directory(dirpath)
    return dirpath


def read_elviz_CSVs(directory):
    elviz_data = {}
    elviz_files = [filename for filename in os.listdir(directory) if
                   ".csv" in filename]
    for filename in elviz_files:
        print(filename)
        # read the dataframe from the csv
        df = read_elviz_CSV("./data/" + filename)
        df['Log10 Average fold'] = numpy.log10(df['Average fold'])
        elviz_data[filename] = df
    return elviz_data


def read_pickle_or_CSVs(pickle_filename, CSV_directory):
    # if the pickle data file exists containing the individual data frames
    # in a list and the combined dataframe then skip loading the CSVs 
    # individually and load the pickle
    if os.path.isfile(pickle_filename):
        print("reading %s for previously parsed data" % pickle_filename)
        with open(pickle_filename, 'rb') as file:
            elviz_data = pickle.load(file)
            combined_df = pickle.load(file)
    else:
        # OK, no pickle found, do it the hard way
        print("reading in all Elviz CSV files")
        elviz_data = read_elviz_CSVs(CSV_directory)
        # assemble the uber frame
        print("concatenating data frames prior to normalization")
        # create a combined dataframe from all the CSV files
        combined_df = pd.concat(elviz_data.values())
        # save the two new objects to a pickle for future use
        with open(pickle_filename, 'wb') as file:
            pickle.dump(elviz_data, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(combined_df, file, pickle.HIGHEST_PROTOCOL)

    return [elviz_data, combined_df]
