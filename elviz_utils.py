import pandas as pd

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
    return pd.read_csv('./data/sample_meta_info.tsv', sep='\t', dtype=IMPORT_METAINFO_TYPES)
