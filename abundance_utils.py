import os
import pandas as pd
import re
from elviz_utils import IMPORT_DATA_TYPES

def get_elviz_filenames(main_dir):
    """
    Return a list of the Elviz filenames
    :return:
    """
    search_dir = os.path.join(main_dir,  'raw_data')
    elviz_files = [f for f in os.listdir(search_dir) if '.csv' in f]
    return elviz_files


def read_elviz_csv(filename, filepath):
    """
    Read in a csv from elviz, extract the ID and project number,
    then reduce to one row per taxonomy.
    Not using Dave's elviz_utils one b/c I get "usecols" arg in
    Python 2.7

    :param filename: name of file to read
    :param filepath: path to find file at
    :return:
    """

    # The first column of some files has funny characters that can
    # be avoided by skipping import of the 0th column: range(1,20)
    # see elviz_1056076_32_HOW6.csv":
    # 'usecols' argument doesn't work in Dave's Python 3!
    df = pd.read_csv(filepath + '/' + filename, sep=",",
                     dtype=IMPORT_DATA_TYPES, usecols=range(1, 20))
    df.fillna("", inplace=True)

    return df


def reduce_elviz_to_genus_rpk(df):
    """
    Take in a raw elviz DataFrame and return a Pandas array with a sum
    of reads per kilobase at the genus level.

    :param df:
    :return:
    """

    # get rid of undesired columns because the data is big.
    # Note: this is not necessary for reduction by groupby below.
    columns_to_drop = ['datasetId', 'contigId',
                       'Average fold',
                       'Reference GC', 'Covered percent', 'Covered bases',
                       # 'Plus reads', 'Minus reads', # keep for read counting
                       'Median fold',
                       'Read GC', 'Complete Lineage', 'IMG scaffold_oid']
    for c in columns_to_drop:
        if c in df.columns:
            del df[c]

    # rename NaN at Genus to "other"  
    df.Genus.replace('', "other", inplace=True)

    # sum the reads that mapped to the plus strand and minus strand
    df['reads'] = df['Plus reads'] + df['Minus reads']

    # Sum lengths and reads.
    df = df.groupby(['Kingdom', 'Phylum', 'Class', 'Order',
                     'Family', 'Genus'])['Length', 'reads'].sum()

    # sort to put the highest ones at the top (useful for debugging)
    df.sort_values(by='reads', axis=0,
                   inplace=True, ascending=False)

    # rename our new measure of abundance: 
    df.rename(columns={'reads': 'sum of reads'},
              inplace=True)
    return df.reset_index()


def normalize_groupby(group, column):
    """
    Normalize all the sum of column to 1 for a given groupy object.

    :param group: group to normalize by.  E.g. #_HOW#
    :param column: coulumn to sum within.
    """
    fold = group[column]
    group[column] = fold / sum(fold)
    return group


def read_and_reduce_elviz_csv(filename, filepath, sample_info):
    """
    Read in a csv from elviz, extract the ID and project number, 
    then reduce to one row per taxonomy.

    :param filename: filename to raw Elviz .csv
    :param filepath: filepath where file lives
    :param sample_info: sample meta_info to use.
    :return:
    """

    # The first column of some files has funny characters that can
    # be avoided by skipping import of the 0th column: range(1,20) 
    # see elviz_1056076_32_HOW6.csv":
    df = read_elviz_csv(filename=filename, filepath=filepath)

    # reduce to one row per genus with reduce_elviz_to_genus_rpk()
    # Need to use DataFrame to convert the seris back to a 
    # dataframe if you want to add an ID column.
    df = pd.DataFrame(reduce_elviz_to_genus_rpk(df))

    # Add a column for the JGI project number.
    df['project'] = project_number_from_filename(filename)

    sample_info.set_index(['project'])
    # merge to get sample_info on
    df = df.groupby('project').apply(normalize_groupby,
                                     'sum of reads')

    df = pd.merge(df, sample_info, how='left')

    # rename column to abundance since we normalized it.
    df.rename(columns={'sum of reads': 'fraction of reads'}, inplace=True)

    # sort so most abundant is on top. 
    df.sort_values(by='fraction of reads', axis=0, ascending=False, inplace=True)

    return df


def read_and_reduce_all(filename_list, filepath, sample_info):
    """

    :param filename_list: list of raw Elviz .csv filenames
    :param filepath: location of Elviz .csv files
    :param sample_info: sample metainfo Pandas DataFrame
    :return:
    """
    # Load in a first dataframe.  Will append the rest to it. 
    dataframe = read_and_reduce_elviz_csv(filename=filename_list[0],
                                          filepath=filepath,
                                          sample_info=sample_info)
    # Loop over the rest of the files. 
    for f in filename_list[1:]:
        # read and write a csv for that sample
        df_to_add = read_and_reduce_elviz_csv(filename=f,
                                              filepath=filepath,
                                              sample_info=sample_info)
        dataframe = dataframe.append(df_to_add)
        print(f)  # prints filename

    # make sure all the samples are there
    if len(dataframe.project.unique()) != 88:
        print("Warning!  only {} samples loaded.".format(
            len(dataframe.ID.unique())))

    # make sure all the fraction of reads add up to 1 for each sample ID
    for t, d in dataframe.groupby('ID'):
        if abs(1 - d['fraction of reads'].sum()) > 0.01:
            print('warning: fraction of reads(s) may not sum to 1')

    return dataframe


def project_number_from_filename(s):
    """
    Get the project numer as an integer from JGI-provided data.

    E.g. 'elviz-contigs-1056013.csv' --> 1056013
    :param s: filename (string) to input
    :return: string of numbers for Elviz project number.
    """
    # print(s)
    return int(re.search('elviz-contigs-([0-9]+).csv', s).group(1))


def prepare_excel_writer_dict(dataframe, filepath, by_genus=False):
    """
    Prepare a dictionary of Pandas ExcelWriters for each replicate at each
    oxygen level.

    :param dataframe: aggregated fraction of reads dataframe
    :param filepath: filepath to store .xlsx files
    :param by_genus: whether data has been reduced to genus level before
    aggregating.
    :return:
    """
    writer_dict = {}
    for oxy in dataframe['oxy'].unique():
        for rep in dataframe['rep'].unique():
            if by_genus:
                filename = 'elviz--Genus_only--{}O2_rep_{}.xlsx'.format(oxy,
                                                                        rep)
            else:
                filename = 'elviz--{}O2_rep_{}.xlsx'.format(oxy, rep)
            filename = filepath + "/" + filename
            writer_dict[(oxy, rep)] = \
                pd.ExcelWriter(filename, engine='xlsxwriter')
    return writer_dict


def write_excel_files(dataframe, filepath, by_genus=False):
    """

    :param dataframe:
    :param filepath:
    :param by_genus:
    :return:
    """

    # prepare excel_writers for each condition; store in a dictionary. 
    # e.g. {('High', 4): <pandas.io.excel._XlsxWriter object at 0x11291e490>,
    #  ('Low', 1): <pandas.io.excel._XlsxWriter object at 0x114156110>, ...}
    writer_dict = prepare_excel_writer_dict(dataframe=dataframe,
                                            filepath=filepath,
                                            by_genus=by_genus)

    # The groupby returns tuples of the conditions.
    by_repl_and_week = dataframe.groupby(['rep', 'week', 'oxy'])

    # loop over these tuples. 
    for (rep, week, oxy), d in by_repl_and_week:
        # use the writer that matches the replicate:
        # print(rep, week, oxy)
        writer = writer_dict[(oxy, rep)]
        sheet_name = oxy + '_O2' + "_rep_" + str(rep) + "_week_" + str(week)
        d.reset_index()
        d.to_excel(writer, sheet_name=sheet_name, index=False)

    # close each writer.  This saves them. 
    for w_dict in writer_dict.values():
        w_dict.close()


def reduce_to_genus_only(dataframe):
    dataframe_genus = dataframe.copy()
    dataframe_genus = dataframe_genus.groupby(['ID', 'rep',
                                               'week', 'oxy',
                                               'Genus']).sum().reset_index()
    dataframe_genus.sort_values(by=['rep', 'fraction of reads'], inplace=True,
                                ascending=False)
    return dataframe_genus


def filter_by_abundance(data, abundance_column, high, low,
                        taxonomy_column='Genus'):
    """
    Return only rows where the specified taxonomy_column's set of rows have at
    least one value of abundance_column in range(low, high)

    :param data: dataframe to filter
    :param abundance_column: column to filter by.  Genus?
    :param taxonomy_column: column to grab unique values from.  Defaults to
    'Genus' for historical reasons.
    :param high: highest abundance to look for
    :param low: lowest abundance to look for
    :return: dataframe
    """
    # get a list of the taxonomy_column names that meet our criteria.
    tax_column_values_to_keep = \
        data[(data[abundance_column] <= high) &
             (data[abundance_column] >= low)][taxonomy_column].unique()
    print(tax_column_values_to_keep[0:5])
    # Return ALL rows for a taxonomy_column label if any of the rows had an
    # abundance value in the desired range.
    return data[data[taxonomy_column].isin(tax_column_values_to_keep)]
