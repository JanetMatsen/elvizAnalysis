

# %run GLOBALS.py
MAIN_DIR = "../"

import sys
sys.path.append('../')

#import numpy as np
import os
import pandas as pd
import re
import seaborn as sns

import abundance_plot_utils
import elviz_utils

# Aesthetics:
# no titles (J might strip off after it is a PDF)
# italic for genera and families
# Capitalize "Other", "Fraction", "Week"
# Lets say "Replicate", not "rep"
# "O2", not "oxygen"
# On some maps, the numbers look fused, can we do something about this?

# Ideally, I would like to put them all on top of each other as panels in
# one figure, so the size of the font should be the same in all panels.

# Component 1:
# Methylococcales, Methylophilales and Burkholderiales at the order level
# Bacteroidetes at the phylum level.
# Note: "other" not mentioned.
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
# Clarifying 5/25/2016 whether this is all other Methylococcaceae,
# or other taxa in general
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
# ?? Other Methylophilaceae as Other, or all other taxa into "other"?
# !! Not Methylovorous, Methylobacillus!!
# For Methylophilaceae, I only want Methylophilus, Methylotenera and other
# -------------------------------------------------------
# Note: similar figure in ipython_notebooks/plots/mixed_taxonomy
# /Genus-Methylotenera_Methylovorus_Methylophilus_Methylobacillus--rep.pdf
# But Methylobacillus has about 0 in all samples.
# Betaproteobacteria; Methylophilales; Methylophilaceae; Methylotenera
# Betaproteobacteria, Methylophilales, Methylophilaceae, Methylophilus

# METHYLOPHILACEAE={'Genus':['Methylotenera', 'Methylovorus',
#                            'Methylophilus', 'Methylobacillus']}
METHYLOPHILACEAE={'Genus':['Methylotenera', 'Methylophilus']}

# Component 4:
# Split the order Burkholderiales at the 0.1% cutoff and not go to 100%.
# 5/25/2016 update:
# On Burkholderiales, the font is way too small.
# Maybe we should reduce the number of entries to
# Acidovorax, Comamonadaceae, other Burkholderiales.
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


def make_heatmap_for_major_players(taxa_dict):
    # load the whole data set.
    data_reduced = \
        pd.read_csv(MAIN_DIR +
                    "results/reduced_data--all_taxonomy_remains.csv")

    # Make plot w/ default settings
    abundance_plot_utils.heatmap_from_taxa_dict(
        dataframe = data_reduced,
        taxa_dict = taxa_dict,
        facet = 'rep',
        annotate = False,
        main_dir=MAIN_DIR,
        plot_dir='./plots/',
        size_spec=False, aspect_spec=False)


def make_heatmap_for_particular_family_with_other(family, taxa_dict):
    # for Methylococcaceae, which we want an "other" bar
    # Load all the data
    data_reduced = \
        pd.read_csv(MAIN_DIR +
                    "results/reduced_data--all_taxonomy_remains.csv")
    # Trim the dataframe to only that family:
    data_for_family = data_reduced[data_reduced['Family']==family]

    # Generate this plot for the subset of data pertaining to that family:
    abundance_plot_utils.heatmap_from_taxa_dict(
        dataframe = data_for_family,
        taxa_dict = taxa_dict,
        facet = 'rep',
        annotate = False,
        main_dir=MAIN_DIR,
        plot_dir='./plots/',
        size_spec=False, aspect_spec=False,
        summarise_other=True,  # <-- want an "other" bar for the family.
        check_totals_sum_to_1=False)  # <-- don't expect col totals to be 1


def heatmap_burkolderiales(taxa_dict):
    # 5/25/2016 update:
    # On Burkholderiales, the font is way too small.
    # Maybe we should reduce the number of entries to
    # Acidovorax, Comamonadaceae, other Burkholderiales.
    data_reduced = \
        pd.read_csv(MAIN_DIR +
                    "results/reduced_data--all_taxonomy_remains.csv")
    abundance_plot_utils.heatmap_all_below(
        dataframe = data_reduced,
        taxa_dict = taxa_dict,
        plot_dir = './plots/',
        low_cutoff = 0.04)


def make_heatmap_for_predators(taxa_dict):
    # same as for major players, but exclude "other" from heat map.
    # load the whole data set.
    data_reduced = \
        pd.read_csv(MAIN_DIR +
                    "results/reduced_data--all_taxonomy_remains.csv")

    # Make plot w/ default settings
    abundance_plot_utils.heatmap_from_taxa_dict(
        dataframe = data_reduced,
        taxa_dict = taxa_dict,
        facet = 'rep',
        annotate = False,
        main_dir=MAIN_DIR,
        plot_dir='./plots/',
        size_spec=False, aspect_spec=False,
        summarise_other=False)  #<--- how to keep the "other" bar off.


def make_figures():
    # Make the ./plots/ dir if needed
    elviz_utils.prepare_plot_dir(MAIN_DIR)

    # Make figure 1:
    make_heatmap_for_major_players(MAJOR_PLAYERS)
    # Make figure 2:
    make_heatmap_for_particular_family_with_other(family='Methylophilaceae',
                                                  taxa_dict=METHYLOPHILACEAE)
    # Make figure 3:
    make_heatmap_for_particular_family_with_other(family='Methylococcaceae',
                                                  taxa_dict=METHYLOCOCCACEAE)
    # Make figure 4:
    # want a different kind of plot for Burkolderiales:
    heatmap_burkolderiales(BURKOLDERIALES)

    # Make figure 5:
    make_heatmap_for_predators(PREDATORS)


if __name__ == "__main__":
    make_figures()
