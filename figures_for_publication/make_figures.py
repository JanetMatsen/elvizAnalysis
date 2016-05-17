%run GLOBALS.pyi

import matplotlib

#import numpy as np
import os
import pandas as pd
import re
import seaborn as sns

# Component 1:
# Methylococcales, Methylophilales and Burkholderiales at the order level
# Bacteroidetes at the phylum level.
# -------------------------------------------------------
# Similar to /Users/jmatsen/programming/elvizAnalysis/ipython_notebooks
# /plots/mixed_taxonomy
# /Order-Burkholderiales_Methylophilales_Methylococcales
# --Phylum-Bacteroidetes--rep.pdf
MAJOR_PLAYERS={'Phylum':['Bacteroidetes'],
               'Order':['Burkholderiales','Methylophilales',
                        'Methylococcales']}
# todo: get an "other" bar/row


# Component 2:  Methylococcaceae
# Split the family Methylococcaceae into the following genera:
# Methylobacter, Methylosarcina, Methylovulum and Methylomonas. 
# All others should be lumped into 'other'.
# -------------------------------------------------------
# Similar to /Users/jmatsen/programming/elvizAnalysis/ipython_notebooks
# /plots/mixed_taxonomy/Genus-Methylobacter_Methylovulum_Methylomonas
# _Methylomicrobium_Methyloglobulus_Methylococcus_Methylocaldum
# _Methylosarcina--rep.pdf
# But the ones missing in this new lis had low abundance:
# Gammaproteobacteria, Methylococcales, Methylococcaceae, Methylobacter
# Gammaproteobacteria, Methylococcales, Methylococcaceae, Methylosarcina
# Gammaproteobacteria, Methylococcales, Methylococcaceae, Methylovulum
# Gammaproteobacteria, Methylococcales, Methylococcaceae, Methylomonas

METHYLOCOCCACEAE = {'Genus': ['Methylobacter', 'Methylovulum',
                              'Methylomonas', 'Methylosarcina']}



# Component 3: 
# Split the family Methylophilaceae into the genera Methytlotenera
# and Methylophilus, the rest should be 'other'.
# -------------------------------------------------------
# Note: similar figure in ipython_notebooks/plots/mixed_taxonomy
# /Genus-Methylotenera_Methylovorus_Methylophilus_Methylobacillus--rep.pdf
# But Methylobacillus has about 0 in all samples.
# Betaproteobacteria; Methylophilales; Methylophilaceae; Methylotenera
# Betaproteobacteria, Methylophilales, Methylophilaceae, Methylophilus

METHYLOPHILACEAE={'Genus':['Methylotenera', 'Methylovorus',
                           'Methylophilus', 'Methylobacillus']}
# todo: get an "other" bar/row

# Component 4:
# Split the order Burkholderiales at the 0.1% cutoff and not go to 100%.
# -------------------------------------------------------
# Note: this is about the same as /Users/jmatsen/programming/elvizAnalysis
# /ipython_notebooks/plots/mixed_taxonomy/Order-Burkholderiales
# _Methylophilales_Methylococcales--Phylum-Bacteroidetes--rep.pdf
# todo: make sure she doesn't want an "other"
BURKOLDERIALES={'Order':['Burkholderiales']}


# Component 5: 
# A two-bar graph with the Y axis set at 10% of total, representing the
# orders Bdellovibrionales and Myxococcales. Thanks!
# -------------------------------------------------------
# The Bdellovibrionaceae are a family of Proteobacteria. They include genera,
# such as Bdellovibrio and Vampirovibrio, which are unusual parasites that
# enter other bacteria.
# The myxobacteria are a group of bacteria that predominantly live in the
# soil and feed on insoluble organic substances. The myxobacteria have very
# large genomes, relative to other bacteria e.g. 9-10 million nucleotides.
PREDATORS={'Order':['Bdellovibrionales', 'Myxococcales']}


def make_figures():
    data_reduced = \
        pd.read_csv(MAIN_DIR +
                    "/results/reduced_data--all_taxonomy_remains.csv")
    taxa_dicts = [MAJOR_PLAYERS, METHYLOCOCCACEAE, METHYLOPHILACEAE,
                  BURKOLDERIALES, PREDATORS]
    for taxa_dict in taxa_dicts:

    pass

if __name__ == "__main__":
    make_figures()
