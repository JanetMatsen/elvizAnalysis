import pandas as pd
import math

# Pandas is happy when you tell it the data types 
# for each column in the raw Elviz .csv files
IMPORT_DATA_TYPES = {'datasetId':'int', 
                    'contigId':'str', 
                    'Average fold':'float',
                    'Length':'int', 
                    'Reference GC':'float',
                    'Covered percent':'float', 
                    'Covered bases':'int', 
                    'Plus reads':'int', 
                    'Minus reads':'int',
                    'Median fold':'float', 
                    'Read GC':'float', 
                    'Complete Lineage':'str',
                    'IMG scaffold_oid':'str',
                    'Kingdom':'str', 
                    'Phylum':'str',
                    'Class':'str', 
                    'Order':'str', 
                    'Family':'str', 
                    'Genus':'str',
                    'Species':'str'
                    }
IMPORT_METAINFO_TYPES = {'ID':'str',
                         'oxy':'str',
                         'rep':'int',
                         'week':'int',
                         'project':'int'}


def read_sample_info():
    ''' 
    Read in sample_meta_info.tsv using particular dtypes
    '''
    print 'testing'
    return pd.read_csv('./data/sample_meta_info.tsv', 
            dtype=IMPORT_METAINFO_TYPES, 
            sep='\t')


def read_elviz_CSV(filename):
    df = pd.read_csv(filename, sep=",", dtype=IMPORT_DATA_TYPES)
    # repalce nans with ""
    df.fillna("", inplace=True)
    return df


def read_elviz_CSVs(directory, logfold=True):
    elviz_data = { }
    elviz_files = [filename for filename in os.listdir(directory) if ".csv" in filename]
    for filename in elviz_files:
        print(filename)
        # read the dataframe from the csv
        df = read_elviz_CSV("./data/" + filename)
        if logfold:
            df['Log10 Average fold'] = math.log(df['Average fold'], 10)
        elviz_data[filename] = df
    return elviz_data


def read_pickle_or_CSVs(pickle_filename, CSV_directory):
    # if the pickle data file exists containing the individual data frames
    # in a list and the combined dataframe then skip loading the CSVs 
    # individually and load the pickle
    if os.path.isfile(DATA_PICKLE):
        print("reading %s for previously parsed data" % DATA_PICKLE)
        with open(DATA_PICKLE, 'rb') as file:
            elviz_data = pickle.load(file)
            combined_df = pickle.load(file)
    else:
        # OK, no pickle found, do it the hard way
        print("reading in all Elviz CSV files")
        elviz_data = read_elviz_CSVs("./data/")
        # assemble the uber frame
        print("concatenating data frames prior to normalization")
        # create a combined dataframe from all the CSV files
        combined_df = pd.concat(elviz_data.values())
        # save the two new objects to a pickle for future use
        with open(DATA_PICKLE, 'wb') as file:
            pickle.dump(elviz_data, file, pickle.HIGHEST_PROTOCOL)
            pickle.dump(combined_df, file, pickle.HIGHEST_PROTOCOL)

    return [elviz_data, combined_df]



