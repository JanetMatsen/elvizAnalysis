import numpy as np
import numpy as np
import re
import pandas as pd

from elviz_utils import IMPORT_DATA_TYPES
from elviz_utils import IMPORT_METAINFO_TYPES


IMPORT_METAINFO_TYPES = {'ID':'str',
                         'oxy':'str',
                         'rep':'int',
                         'week':'int',
                         'project':'int'}


def reduce_elviz_to_genus_rpk(df):
    '''
    Take in a raw elviz DataFrame and return a Pandas array with a sum at the genus level
    '''
    # get rid of undesired columns because the data is big.
    columns_to_drop = ['datasetId', 'contigId', 
                       'Average fold', 
                       'Reference GC','Covered percent','Covered bases',
                       #'Plus reads', 'Minus reads', # keep for read counting 
                       'Median fold',
                       'Read GC','Complete Lineage', 'MG scaffold_oid']
    for c in columns_to_drop:
        if c in df.columns:
            del df[c]

    # rename NaN at Genus to "other"  
    df.Genus.replace(np.nan, "other", inplace=True)
    
    # calculate RPK: sum the reads that mapped to the plus strand and minus strand
    df['reads per kilobase'] = (df['Plus reads'] + \
                                df['Minus reads'])/(df['Length']/1000.)

    df = df.groupby(['Kingdom','Phylum','Class','Order','Family','Genus']) \
                        ['reads per kilobase'].sum().reset_index()

    # rename our new measure of abundance: 
    df.rename(columns={'reads per kilobase': 'sum of reads per kilobase'}, inplace=True)
    return df 


def norm_by_ID(group):
    ''' 
    Normalize all the sum of reads per kilobase to 1 for 
    '''
    fold = group['sum of reads per kilobase']
    group['sum of reads per kilobase'] = fold/sum(fold)
    return group


def read_and_reduce_elviz_csv(filename, filepath, sample_info):
    '''
    Read in a csv from elviz, extract the ID and project number, 
    then reduce to one row per phylogeny.
    '''

    # The first column of some files has funny characters that can
    # be avoided by skipping import of the 0th column: range(1,20) 
    # see elviz_1056076_32_HOW6.csv":
    df = pd.read_csv(filepath + '/' + filename, sep=",", 
                    dtype=IMPORT_DATA_TYPES, usecols=range(1,20))

    # reduce to one row per genus with reduce_elviz_to_genus_rpk()
    # Need to use DataFrame to convert the seris back to a 
    # dataframe if you want to add an ID column.
    df = pd.DataFrame(reduce_elviz_to_genus_rpk(df))

    # add an ID label.  (don't use 'id'; that is a Python built-in.)
    df['ID'] = id_from_filename(filename)
    #df['project'] = project_number_from_filename(filename)
    
    sample_info.set_index(['ID'])
    # merge to get sample_info on
    df = pd.merge(df, sample_info, how='left')
    df = df.groupby('ID').apply(norm_by_ID)
    
    # after norm_by_ID is applied, 'Average fold' is now a pooled number.
    # rename column to abundance since we normalized it. 
    df.rename(columns={'sum of reads per kilobase':'abundance'}, inplace=True)
    
    # sort so most abundant is on top. 
    #print df.head()
    df.sort_values(by='abundance', axis=0, ascending=False, inplace=True)
    
    return df


def read_and_reduce_all(filename_list, filepath, sample_info):
    # Load in a first dataframe.  Will append the rest to it. 
    dataframe = read_and_reduce_elviz_csv(filename = filename_list[0], 
                                             filepath=filepath,
                                             sample_info=sample_info)
    # Loop over the rest of the files. 
    for f in filename_list[1:]:
        # read and write a csv for that sample
        df_to_add = read_and_reduce_elviz_csv(filename = f, 
                                              filepath=filepath,
                                              sample_info=sample_info)
        dataframe = dataframe.append(df_to_add)
        print f  # prints filename

    # define a sample name
    dataframe['sample_name'] = 'replicate '+ \
        dataframe['rep'].astype(str) +": " + \
        dataframe['oxy'] + ' O2'

    # make sure all the samples are there
    if len(dataframe.ID.unique()) != 88:
        print "Warnng!  only {} samples loaded.".format( \
                                    len(dataframe.ID.unique()))
    
    # make sure all the abundances add up to 1 for each sample ID
    for t, d in dataframe.groupby('ID'):
            if abs(1 - d['abundance'].sum()) >0.01:
                print 'warning: abundance(s) may not sum to 1'
    
    return dataframe


def id_from_filename(s):
    return re.search('[\w]+_[0-9]+_([0-9]+_[HL]OW[0-9]+).csv', s).group(1)


def project_number_from_filename(s):
    return re.search('[\w]+_([0-9]+)_[0-9]+_[HL]OW[0-9]+.csv', s).group(1)


def prepare_excel_dictionary(dataframe):
    """
    Make a dictionary like
    {('High', 1): 'elviz_binned--HighO2_rep1.xlsx',
     ('High', 2): 'elviz_binned--HighO2_rep2.xlsx', ...}
    """
    by_repl_and_week = dataframe.groupby(['rep','week','oxy'])
    # write same dictionary in a loop
    excel_files = {}
    for ox in dataframe['oxy'].unique():
        for re in dataframe['rep'].unique():
            excel_files[(ox, re)] =  'elviz_binned--{}O2_rep{}.xlsx'.format(ox, re)
    return excel_files

def prepare_excel_writer_dict(dataframe, filepath, by_genus=False):
    writer_dict = {}
    for oxy in dataframe['oxy'].unique():
        for rep in dataframe['rep'].unique():
            #print ox, re
            if by_genus:
                filename = 'elviz--Genus_only--{}O2_rep_{}.xlsx'.format(oxy, rep)
            else:
                filename = 'elviz--{}O2_rep_{}.xlsx'.format(oxy, rep)
            filename = filepath + "/" + filename
            #print filename
            writer_dict[(oxy, rep)] = pd.ExcelWriter(filename, engine='xlsxwriter')
            #print writer_dict
    return writer_dict


def write_excel_files(dataframe, filepath, by_genus=False):

    # Initialize a dictionary of excel filenames using the oxygen and replicate #s
    #excel_filenames = prepare_excel_dictionary(dataframe = dataframe)

    # prepare excel_writers for each condition; store in a dictionary. 
    # e.g. {('High', 4): <pandas.io.excel._XlsxWriter object at 0x11291e490>,
    #  ('Low', 1): <pandas.io.excel._XlsxWriter object at 0x114156110>, ...}
    writer_dict = prepare_excel_writer_dict(dataframe = dataframe,
                                            filepath=filepath,
                                            by_genus=by_genus)
    #print 'writer_dict:', 
    #print writer_dict

    # The groupby returns tuples of the conditions. 
    by_repl_and_week = dataframe.groupby(['rep','week','oxy'])

    # loop over these tuples. 
    for (rep, week, oxy), d in by_repl_and_week:
        # use the writer that matches the replicate:
        #print rep, week, oxy
        writer = writer_dict[(oxy, rep)]
        sheet_name = oxy + '_O2' +"_rep_" + str(rep)+"_week_" + str(week)
        #print sheet_name
        d.reset_index()
        #print d.columns
        #del d['index']
        #print d.columns
        d.to_excel(writer, sheet_name=sheet_name, index=False)
        #writer.save()
        #print ""
    
    # close each writer.  This saves them. 
    for w_dict in writer_dict.values():
        w_dict.close()

def reduce_to_genus_only(dataframe):
    dataframe_genus = dataframe.copy()
    #del dataframe_genus['reads per kilobase']
    print dataframe_genus.head()
    print dataframe.columns
    dataframe_genus = dataframe_genus.groupby(['ID','rep',
                'week','oxy','Genus']).sum().reset_index()
    dataframe_genus.sort_values(by=\
        ['rep', 'abundance'], inplace=True, ascending=False)
    return dataframe_genus
