import itertools
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import abundance_utils
import elviz_utils

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# Need to use LaTeX to get italic fonts.
rc('text', usetex=True)
rc('font', family='sans-serif')
# http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]


def axd_portrait(axs):
    """
    axis dictionary for portrait pages.
    """
    return {('Low', 1): axs[0, 0], ('High', 1): axs[0, 1],
           ('Low', 2): axs[1, 0], ('High', 2): axs[1, 1],
           ('Low', 3): axs[2, 0], ('High', 3): axs[2, 1],
           ('Low', 4): axs[3, 0], ('High', 4): axs[3, 1]}


def axd_landscape(axs):
    """
    axis dictionary for landscape pages.
    """
    return {('Low', 1): axs[0, 0], ('Low', 2): axs[0, 1], ('Low', 3): axs[0, 2], ('Low', 4): axs[0, 3],
            ('High', 1): axs[1, 0], ('High', 2): axs[1, 1], ('High', 3): axs[1, 2], ('High', 4): axs[1, 3]}


def add_vline_to_all_subplots(fig, x, ymin, ymax, color='#636363'):
    for ax in fig.axes:
        # It counts the x positions from 1, not from 4.
        ax.axvline(x=x, ymin=ymin, ymax=ymax, color=color)
    return fig


def plot_heatmap_genus(dataframe, high, low, oxy, rep, plot_dir):
    """
    Make a heatmap at Genus, using oganisms withing the specified abundance
    cutoffs.

    :param dataframe: dataframe to pass
    :param high: highest abundance to include genera for
    :param low: lowes abundance to include genera for
    :param oxy: oxygen tension, "Low" or "High"
    :param rep: replicate (1-4)
    :param plot_dir: directory to save plots in.
    :return:
    """
    # Correction: 'other' should be 'other and unknown'
    dataframe.ix[dataframe['Genus'] == 'other', 'Genus'] = r'other \& unknown'

    # get rid of oxygen levels and replicates if specified.
    if oxy is not 'all':
        print("keep only {} oxygen samples".format(oxy))
        dataframe = dataframe[dataframe['oxy'] == oxy]
    if rep is not 'all':
        print("keep only replicate levels:", rep)
        dataframe = dataframe[dataframe['rep'].isin(rep)]
    dataframe = abundance_utils.filter_by_abundance(
        data=dataframe,
        abundance_column='fraction of reads',
        high=high, low=low)
    dataframe['facet_replicate'] = 'replicate ' + dataframe['rep'].astype(str)

    # make height of the plot a function of the number of rows (Genera):
    num_data_rows = len(dataframe['Genus'].unique())
    plot_size = 2 + num_data_rows / 7
    plot_aspect = 2
    if num_data_rows > 6:
        plot_aspect = .85
    if num_data_rows > 9:
        plot_aspect = .65
    if num_data_rows > 9:
        plot_aspect = .6

    def facet_heatmap(data, **kws):
        """
        Used to fill the subplots with data.

        :param data:
        :param kws:
        :return:
        """

        facet_data = data.pivot(index='Genus', columns='week',
                                values='fraction of reads')
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)

    with sns.plotting_context(font_scale=7):
        g = sns.FacetGrid(dataframe, col='facet_replicate',
                          margin_titles=True,
                          size=plot_size, aspect=plot_aspect)
        g.set_xticklabels(rotation=90)

    # Create a colorbar axes
    cbar_ax = g.fig.add_axes([.94, .3, .02, .4], title='fraction \n of reads')

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0,
                        # specify vmax = max abundance seen or each will
                        # have its own scale (and you might not know it!)
                        vmax=dataframe['fraction of reads'].max(),
                        )

    g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=18)
    g.set_axis_labels('week')

    # Add space so the colorbar doesn't overlap the plot
    g.fig.subplots_adjust(right=.9)

    # add a supertitle, you bet.
    plt.subplots_adjust(top=0.80)
    supertitle = str(low) + r' $<$ fraction of reads $<$ ' + str(
        high) + ', {} oxygen'.format(oxy)
    g.fig.suptitle(supertitle, size=18)

    # write a filename and save.
    filename = oxy + "_oxygen--{0}_to_{1}_abundance".format(low, high)
    print('filename:', filename)

    plot_dir = elviz_utils.prepare_plot_dir(plot_dir)

    # save figure
    g.savefig(plot_dir + filename + '.pdf')


def subset_on_taxonomy(dataframe, taxa_level, name):
    """
    Return only rows of the datframe where the value in column taxa_level
    matches the specified name.

    :param dataframe: Pandas DataFrame with columns like 'Kingdom',
    'Phylum', 'Class', ...
    :param taxa_level: a taxagenetic label such as "Genus" or "Order"
    :param name: taxa_level name to match
    :return: subset of Pandas DataFrame matching the selection
    """
    #print(dataframe.columns)
    return dataframe[dataframe[taxa_level] == name]


def other_taxonomy_levels(level):
    """
    return the name of all the taxonomic levels *not* specified by level.

    Handy if you want to drop all the other columns.

    :param level: string string string string
    :return:
    """
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    levels.remove(level)
    return levels


def taxonomy_levels_above(taxa_level):
    """
    E.g. 'Order' --> ['Kingdom', 'Phylum', 'Class']
    """
    p_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    position_of_taxa_level = p_levels.index(taxa_level)
    return p_levels[0:position_of_taxa_level]


def taxonomy_levels_below(taxa_level):
    """
    E.g. 'Order' --> ['Family', 'Genus']
    """
    p_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    position_of_taxa_level = p_levels.index(taxa_level)
    return p_levels[position_of_taxa_level + 1:]


def sum_on_taxonomy(dataframe, taxa_level, name):
    """
    Sum children rows for taxa_level == name in dataframe; returns one row
    per sample.

    E.g. if you pass taxa_level = 'Phylum' and name = 'Bacteroidetes',
    you will get one row with columns ['taxagenetic label', 'name',
    'sum of abundances'] *per* sample.
    The sum of abundances results from groupby aggregations within the
    sample group.

    :param dataframe: dataframe with taxagenetic levels as columns
    :param taxa_level: A taxagenetic level in the set
    ['Kingdom', 'Phylum', 'Class', 'Order', 'Family']
    :param name: name for the desired value of taxa_level
    :return: Pandas DataFrame with aggregated rows based on taxa_level == name
    """

    relevant_rows = subset_on_taxonomy(dataframe=dataframe,
                                       taxa_level=taxa_level,
                                       name=name)
    # Collaps all taxonomy below:
    for col_name in other_taxonomy_levels(taxa_level):
        del relevant_rows[col_name]

    # sum all
    aggregated_rows = relevant_rows.groupby(
        [taxa_level, 'ID'])['fraction of reads'].sum()

    return aggregated_rows


def delete_rows_for_taxa(dataframe, taxa_level, taxa_name):
    """
    return a copy of the dataframe with the taxa name at the specified taxa
    level removed.

    For use in aggregating_mixed_taxonomy: allows for keeping track of the
    leftovers.

    :param dataframe: dataframe to copy and delete rows from
    :param taxa_level: taxa level to look for the name in
    :param taxa_name: rows with this name at the specified taxa_leven are
    deleted
    :return: a dataframe with less rows
    """
    df = dataframe.copy()
    return df[df[taxa_level] != taxa_name]


def collapse_unused_taxa_into_other(dataframe):
    """
    After taxa have been cherry picked out of a dataframe, we want to collapse
    the rest of the taxa into "other".

    :param dataframe: dataframe of leftovers to consider "other".
    :return: dataframe with columns: ID, abundance sum, taxonomic level,
    taxonomic name
    """
    # We want to append this result onto something like:
    #        ID  abundance sum taxonomic level   taxonomic name
    # 100_LOW12       0.084171           Order  Burkholderiales
    #print(dataframe.head())

    # All we care about is the ID and the abundance sum left.
    df = dataframe[['ID', 'fraction of reads']]
    # the groupby sum returns a Pandas series, so we ask for a dataframe:
    df = pd.DataFrame(df.groupby(['ID'])['fraction of reads'].sum())
    # Now we have something like
    #         ID  fraction of reads
    #  100_LOW12           0.207207
    #  103_HOW12           0.207381

    # Rename "fraction of reads" to "abundacne sum" to match the results this
    # dataframe will be appended to.
    df.rename(columns={'fraction of reads': 'abundance sum'},
              inplace=True)

    # Attach columns to clarify these are "other" counts.
    df['taxonomic level'] = "aggregate"
    df['taxonomic name'] = r"unknown \& other"

    df.reset_index(inplace=True)
    #print(df.head())
    return df


def aggregate_mixed_taxonomy(dataframe, taxa_dict, main_dir='./',
                             summarize_other=True, check_totals_sum_to_1=True):
    """
    Summarise abundances based on cherry-picked taxonomic abundances,
    perhaps mixed at different levels.

    Loop over the different taxonomic levels specified in a dictionary of
    taxonomic level keys and name pairs.  Reduce using sum_on_taxonomy()
    and store that result in a list.  Concatenate the lists into one DataFrame
    for return.

    In order to keep track of what has not been cherry picked, we cherry pick
    out of a copy of the dataframe, deleting those rows after they are used.
    Then what is left can be lumped into "other".  This also helps (though
    does not completely solve) the issue of an invalad taxonomy dict being
    passed as an argument.

    To get an "other" for a given taxonomic level, pass in a DataFrame that
    is already restricted to the taxonomic level you are looking for
    (e.g. Methylococcaceae).  Then the taxa you *aren't* picking out will be
    represented by "other", instead of all other taxa at all other taxonomic
    levels.

    :param dataframe: dataframe containing all the data to pick through
    :param taxa_dict: a dictionary with taxonomic levels as keys and
    names as values.  E.g. {'Phylum':['Bacteroidetes'],
    'Order':['Burkholderiales','Methylophilales', 'Methylococcales']}
    :param main_dir: directory where the data is stored.  This argument was
    added so jupyter notebooks could be run in a sub-directory.
    :param summarize_other: include an "other" row per sample?  (or omit)
    :return:
    """

    # Make a copy of the dataframe so we can sum the leftovers into an
    # "other" category.
    df = dataframe.copy()

    # First make sure all of the taxa name, value pairs are valid:
    # There is a unit test for this in:
    # unit_tests.test_abundance_plot_utils.
    # testDeleteRowsForTaxa#test_invalid_taxa_dict
    for key in taxa_dict.keys():
        for name in taxa_dict[key]:
            assert (key in df.columns), "column {} doesn't exist!".format(key)
            # Note: to check for a value in a series, use set(Series)
            assert (name in set(df[key])), \
                'Value "{}" in column "{}" does not exist. \n' \
                'Check spelling and conflict with other taxa in the ' \
                'taxa dict.'.format(name, key)

    # TODO: need to start at most general taxonomic level if you want to
    # check for errors in the taxa dict.
    # Currently no error will be thrown if you pick a Genera out, then
    # subsequently select a Family that would have included that Genera.
    # Also note that if you pass the same key twice in a dict, Python keeps
    # only the second key:value pair.  E.g.:
    # print({'Genus': ['Orcinus'], 'Genus': ['ABCD']})

    reduced_data = []
    for key in taxa_dict.keys():
        for name in taxa_dict[key]:
            # Get one row per week/oxygen condition:
            reduced_rows = sum_on_taxonomy(dataframe=df,
                                           taxa_level=key,
                                           name=name)
            # check that you got some rows.  Might be a typo if not!
            assert(reduced_rows.shape[0] > 0), \
                'found no rows for {} = "{}"'.format(key, name)

            df = delete_rows_for_taxa(dataframe=df,
                                      taxa_level=key, taxa_name=name)

            # the index needs to be dropped but it is stored below as
            # 'taxonomic level' and 'taxonomic name'
            # I haven't been able to reset_index on this series to drop the
            # index so I'm doing it this way:
            reduced_rows = reduced_rows.reset_index()
            del reduced_rows[key]
            # make a new dataframe out of it.
            reduced_data.append(
                pd.DataFrame(
                    {'taxonomic level': key,
                     'taxonomic name': name,
                     'abundance sum': reduced_rows['fraction of reads'],
                     'ID': reduced_rows['ID']}))
    # Concatenate data included in the taxa_dict
    # Has form like:
    #        ID  abundance sum taxonomic level   taxonomic name
    # 100_LOW12       0.084171           Order  Burkholderiales
    dataframe_of_keepers = pd.concat(reduced_data)
    #print(dataframe_of_keepers.head())

    # Aggregate the leftovers into an "other" column, with headers to match
    # dataframe_of_keepers
    # TODO: aggregate into "other"
    dataframe_of_leftovers = collapse_unused_taxa_into_other(df)

    if summarize_other:
        # Merge the keepers and the leftovers.
        result_df = pd.concat([dataframe_of_keepers, dataframe_of_leftovers])
        #print("merged result_df:")
        #print(result_df.head())
        # merge on the sample info.
        result_df = pd.merge(left=result_df,
                             right=elviz_utils.read_sample_info(main_dir))
        # Check that the sum of abundances for each sample is really close
        # to 1:
        if check_totals_sum_to_1:
            sample_sums = result_df.groupby('ID')['abundance sum'].sum()
            assert (sample_sums > 0.999).all()
            assert (sample_sums < 1.001).all()
            #print(sample_sums.head())
    else:
        result_df = pd.merge(left=dataframe_of_keepers,
                             right=elviz_utils.read_sample_info(main_dir))

    return result_df


def taxa_dict_to_descriptive_string(taxa_dict):
    """
    Turn a taxa_dict into a string for plot names and titles.

    :param taxa_dict: a dictionary with taxonomic levels as keys and
    names as values.  E.g. {'Phylum':['Bacteroidetes'],
    'Order':['Burkholderiales','Methylophilales', 'Methylococcales']}
    :return: a string without spaces representing concatenation of the dict.
    """
    # todo: go through highest orders of taxa first.
    # e.g. phylum looped over before order.
    desc_string = ""
    for key in taxa_dict:
        print(key)
        desc_string += key
        desc_string += '-'
        for value in taxa_dict[key]:
            desc_string += value + '_'
        desc_string = desc_string[:-1]
        desc_string += '--'
    # remove last two '--' characters
    desc_string = desc_string[:-2]
    return desc_string


def heatmap_from_taxa_dict(dataframe, taxa_dict,
                           title=False,
                           facet='rep', annotate=False,
                           summarize_other=True,
                           main_dir='./',
                           cap_facet_labels=True,
                           plot_dir='./plots/mixed_taxonomy/',
                           size_spec=False,
                           aspect_spec=False,
                           check_totals_sum_to_1=True,
                           svg=False):
    """
    Make a plot using a taxa_dict.

    The taxa_dict is used to make a summary dataframe using
    aggregate_mixed_taxonomy(), and the result is plotted.

    :param dataframe: dataframe to source all data from
    :param taxa_dict: a dictionary with taxonomic levels as keys and
    names as values.  E.g. {'Phylum':['Bacteroidetes'],
    'Order':['Burkholderiales','Methylophilales', 'Methylococcales']}
    :param facet: The rows to facet the subplots by.  Defaults to replicates,
    so weeks will be the columns.
    :param annotate: print numerical values inside each square?  (Makes big
    plots *really* big; not recommended for default use.
    :param main_dir: main dir to consider "home", so notebooks can be run in
    remote directories.
    :param summarize_other: include a bar for "other"?  (Or just don't show)
    :param plot_dir: path to save plots at, relative to main_dir
    :param size_spec: manually specify the figure size (useful when default
    is ugly)
    :param aspect_spec: manually specify the figure asepct ratio (useful when
    default is ugly
    :return: saves and returns a seaborn heat map
    """

    # Cherry pick out the rows for the specified taxa.
    # If you give conflicting taxa as input, aggregate_mixed_taxonomy() will
    # throw an error.
    plot_data = aggregate_mixed_taxonomy(
        dataframe=dataframe,
        taxa_dict=taxa_dict,
        main_dir=main_dir,
        summarize_other=summarize_other,
        check_totals_sum_to_1=check_totals_sum_to_1)

    # store the maximum abundance level.  We will need to tell all the
    # sub-heat maps to use this same maximum so they aren't each on their
    # own scale.
    max_abundance = plot_data['abundance sum'].max()

    # The data is seperated by these two variables.
    # The one not used as the facet will be used as the columns in the
    # subplot.
    if facet == 'week':

        cols_in_facet = 'rep'
    else:
        cols_in_facet = 'week'

    print('plot_data.head()')
    print(plot_data.head())

    def pivot_so_columns_are_plotting_variable(dataframe, groupby):
        return dataframe.pivot(index='taxonomic name',
                               columns=groupby,
                               values='abundance sum')

    def facet_heatmap(data, groupby, xrotation, **kws):
        """
        Used to fill the subplots with data.

        :param data: dataframe to plot
        :param groupby: column to group on
        :param xrotation: degrees to rotate x labels by
        :param kws: kewyord arguments for plotting
        :return:
        """
        # pivot only supports one column for now.
        # http://stackoverflow.com/questions/32805267/pandas-pivot-on-multiple-columns-gives-the-truth-value-of-a-dataframe-is-ambigu
        facet_data = pivot_so_columns_are_plotting_variable(
            dataframe=data, groupby=groupby)
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)
        g.set_xticklabels(rotation=xrotation)

    # todo: add a label at the bottom like "replicate" or "week"
    # currently replicate is turned into facet_replicate but should just
    # make a label that says replicate.  Week

    # Control plot aesthetics depending on facet option.
    if facet == 'week':
        xrotation = 0
        num_rows = len(plot_data['taxonomic name'].unique())
        size = 2 * 0.2*num_rows
        aspect = 1
        space_for_cbar = 0.85
        x_axis_label = 'replicate'

    else: # (facet = "rep")
        xrotation = 90
        # Calculate the size, aspect depending on the number of
        #  rows per subplot
        num_rows = len(plot_data['taxonomic name'].unique())
        size = 1 + 0.22*num_rows
        aspect = 1.5  # aspect for each sub-plot, not a single tile
        space_for_cbar = 0.85

    if size_spec:
        size = size_spec
    if aspect_spec:
        aspect = aspect_spec

    print(plot_data.head())

    if cap_facet_labels:
        if facet == "rep":
            row_var='$O_2$'
            col_var = 'Week'
            facet_var = "Replicate"
        else:
            print("not set up for facet != rep")
        plot_data = capitalize_some_column_names(plot_data)
        col_var
    else:
        facet_var = 'rep'
        row_var = 'oxy'
        col_var = 'week'

    with sns.plotting_context(font_scale=8):
        g = sns.FacetGrid(plot_data,
                          col=facet_var,
                          row=row_var,
                          size=size,
                          aspect=aspect,
                          margin_titles=True)

    # Add axes for the colorbar.  [left, bottom, width, height]
    cbar_ax = g.fig.add_axes([.92, .3, .02, .4], title='fraction \n of reads')

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax,
                        # NEED vmax = MAX ABUNDANCE or each plot will have
                        # its own color scale!
                        vmin=0, vmax=max_abundance,
                        annot=annotate,
                        groupby=col_var,
                        xrotation=xrotation)

    g.set_axis_labels(col_var)

    # add space for x label
    g.fig.subplots_adjust(bottom=0.2)

    # todo: add an x-label for each facet (I want only 1)
    # g.set_axis_labels(['x label', 'ylabel'])
    # g.fig.subplots_adjust(top=0.2)
    # g.fig.text(0.5, 0.1, s='armadillo') #, *args, **kwargs)
    # g.fig.xlabel('ardvark')

    # Add space so the colorbar doesn't overlap th plot.
    g.fig.subplots_adjust(right=space_for_cbar)
    # todo: still not enough room for
    # Order-Burkholderiales_Methylophilales_Methylococcales--
    # Phylum-Bacteroidetes--rep.pdf

    # Format the y strings in each subplot of the Seaborn grid.
    # Don't put () on the function you are c
    # Todo: make the 2nd argument a function
    y_label_formatter(g, italics_unless_other)

    supertitle = taxa_dict_to_descriptive_string(taxa_dict)
    if title:
        # TODO: they are currently being converted to LaTeX
        # add a supertitle, you bet.
        plt.subplots_adjust(top=0.80)
        g.fig.suptitle(supertitle, size=16)

    # Tight layout --> title and cbar overlap heat maps.  Boo.
    # NO: plt.tight_layout()
    g.fig.subplots_adjust(wspace=.05, hspace=.05)

    # prepare filename and save.
    plot_dir = elviz_utils.prepare_plot_dir(plot_dir)
    print("plot directory: {}".format(plot_dir))
    filepath = plot_dir + supertitle
    filepath += "--{}".format(facet)
    if annotate:
        filepath += "--annotated"
    filepath += ".pdf"
    print(filepath)
    g.fig.savefig(filepath)
    if svg:
        g.fig.savefig(filepath.rstrip("pdf") + "svg")

    return g


def label_from_taxa_colnames(name_list, taxa_name):
    """
    Return a string that compresses a list of strings into a string
    separated by ", ".  It skips over names "other" and "unknown".

    Leaves out detail once "other" and "unknown" are hit.
    So ["Oxalobacteraceae", "other"] --> "Oxalobacteraceae

    e.g. ['Burkholderiales', 'Comamonadaceae, 'other'] -->
        'Burkholderiales_Comamonadaceae_other', or
    ['Burkholderiales', NaN, 'other']

    :param args: a list
    :taxa name: level to use for "other" so you get a heatmap label like
     "other Burkolderiales" instead of "other", which could look like
     "all other taxa summed"

    :return: a string like or "Burkholderiales, Comamonadaceae",
     or "other Burkolderiales"
    """
    # TODO: doesn't make sense to use *args for this function, since we only
    #  ever pass one list.  (One list, Right?)
    # 160524 update: does seem to need to be a *args thing.
    name_string = ""
    for name in name_list:
        if name != 'other':
            if name != 'unknown':
                # if not math.isnan(name):   #np.isnan(name):
                # print('adding name {}'.format(name))
                name_string += name
                name_string += ", "
    # remove the last ", "
    if len(name_string) > 0:
        # print('length of {} is > 0'.format(name_string))
        return name_string[:-2]
    else:
        # If the string has no length, return 'other'
        # print('all fields empty.  returning "?"')
        return 'other {}'.format(taxa_name)


def italics_unless_other(name_list):
    # Examples:  ["Abcd", "Other Efgh "] -->
    # [r"\textit{Abcd}", r"Other \textit{Efgh}"]
    print(name_list)
    # lists are ordered  :)
    new_name_list = []
    for original_name in name_list:
        formatted_name = ""
        strings = original_name.split(" ")
        print("strings split on the space: {}".format(strings))
        # Loop over the words in the string.
        for i, s in enumerate(strings):
            print("single word: {}".format(s))
            # check for Other and other, unknown, and don't italicize those.
            if "ther" in s or 'nknown' in s:
                # don't italacize "other" or "Other"
                formatted_name += s
            elif "iales" in s:
                # check for "iales", which signifies Genus.
                # Don't italicize Genus name.
                formatted_name += s
            else:
                # "my label" --> "\textit{my label}"
                formatted_name += r"\textit{" + s + r"}"

            # manually add back spaces between words
            formatted_name += " "  # will need to strip off the trailing one.

        # Take off the trailing " " character from the label.
        formatted_name.rstrip(" ")
        print("{} ---> {}".format(original_name, formatted_name))
        new_name_list.append(formatted_name)
    print("new name list: {}".format(new_name_list))
    return new_name_list


def y_label_formatter(seaborn_facetgrid_plot, name_format_fun):
    #def y_label_formatter(seaborn_facetgrid_plot, string_function):
    # E.g. change all text that doesn't have "other" in it to italics.
    # Try to control y labels:
    # http://stackoverflow.com/questions/34422944/seaborn-heat-map-customize-a-label
    g = seaborn_facetgrid_plot
    print("play with the y labels now!")
    for ax in g.fig.get_axes():
        for label in ax.get_yticklabels():
            # print(label)
            if "other" in label.get_text(): #== "other Burkholderiales":
                # label.label = "abcd"
                label.set_weight("bold")
                # label.set_color("red")
                # label.set_style('italic')

        # Change the "other" labels to italics for proof of principle.
        # Uses LaTeX, which requires one of a few back-ends.
        # http://stackoverflow.com/questions/11244514/modify-tick-label-text
        original_labels = [item.get_text() for item in ax.get_yticklabels()]
        labels = name_format_fun(original_labels)
        # use these new labels in the figure.
        ax.set_yticklabels(labels)


def capitalize_some_column_names(dataframe):
    """
    Rename columns as per Mila's preferences.

    :param dataframe: dataframe to modify
    :return: dataframe with modified column names.
    """
    if "week" in dataframe.columns:
        dataframe.rename(columns={'week':'Week'}, inplace=True)
    if "oxy" in dataframe.columns:
        dataframe.rename(columns={'oxy': '$O_2$'}, inplace=True)
    if "rep" in dataframe.columns:
        dataframe.rename(columns={'rep': 'Replicate'}, inplace=True)
    return dataframe


def heatmap_all_below(dataframe, taxa_dict, plot_dir, low_cutoff=0.001,
                      cap_facet_labels=True,
                      title=False, svg=False):
    """
    Make a heatmap of all the taxa below the taxa specified in taxa_dict.

    :param dataframe: dataframe of data to harvest excerpts from
    :param taxa_dict: a dictionary with taxonomic levels as keys and
    names as values.  E.g. {'Order':['Burkholderiales']}
    :param plot_dir: path to save plots to, relative to main_dir
    :param main_dir: path to data source, etc.
    :param low_cutoff: lowest abundance to include.  A taxa must be above
    this threshold in at least one sample to be included.
    :return:
    """
    # TODO: this function has a lot of commonality with heatmap_from_taxa_dict
    # and could/should be factored.

    # grab the data for that taxa:
    # for now assume just 1 key and 1 value.
    taxa_level = list(taxa_dict.keys())[0]
    taxa_name = list(taxa_dict.values())[0][0]
    dataframe = dataframe[dataframe[taxa_level] == taxa_name]
    print(dataframe.head())

    # Columns to form a concatenated label from:
    label_cols = taxonomy_levels_below(taxa_level=taxa_level)
    print('label_cols: {}'.format(label_cols))

    # change nan cells to 'unknown'
    dataframe.fillna('unknown', inplace=True)

    # make a summary string representing the taxonomy for everything below

    def label_building_lambda(f, column_value_list, taxa_name):
        """
        Returns a lambda function to make row labels from.
        :param f: function to make a lambda out of.
        :param columns: column names to pass to function f in the lambda
        :return: function
        """
        # * means unpack the list you get from the list comprehension
        print("columns passed: {}".format(column_value_list))
        print("Use those in {}".format(f))
        # Passing a list into label_from_taxa_colnames().
        # Doing a list comprehension on columns.
        # Note that (row[col] for col in columns)) is a generator .
        # building something like label_from_taxa_colnames()
        return lambda row: f([row[col] for col in column_value_list],
                             taxa_name)
        # e.g. makes:
        # my_function([Comamonadaceae, Curvibacter]) from a row of a dataframe
        # and the specification that columns = ['Family', 'Genus']

    # TODO: use the taxa_dict to get the columns to use!
    # make a name_string per row.  It's something like
    # "Comamonadaceae, Curvibacter" or "other"
    dataframe['name_string'] = dataframe.apply(
        label_building_lambda(f=label_from_taxa_colnames,
                              column_value_list=label_cols,
                              taxa_name=taxa_name),
        axis=1)
    print("dataframe.head() for name_string:")
    print(dataframe.head())

    # reduce to only name_string rows with at least one abundance > the
    # threshold set by low_cutoff to we don't have a zillion rows:
    # todo: allow high to change?
    dataframe = \
        abundance_utils.filter_by_abundance(data=dataframe,
                                            abundance_column='fraction of '
                                                             'reads',
                                            high=1,
                                            low=low_cutoff,
                                            taxonomy_column='name_string')

    # Plot as usual, using the stuff developed above.
    # todo: factor some of this??
    def pivot_so_columns_are_plotting_variable(dataframe, groupby):
        return dataframe.pivot(index='name_string',
                               columns=groupby,
                               values='fraction of reads')

    def facet_heatmap(data, groupby, xrotation, **kws):
        """
        Used to fill the subplots with data.

        :param data: dataframe to plot
        :param groupby: column to group on
        :param xrotation:
        :param kws:
        :return:
        """
        # pivot only supports one column for now.
        # http://stackoverflow.com/questions/32805267/pandas-pivot-on-multiple-columns-gives-the-truth-value-of-a-dataframe-is-ambigu
        facet_data = pivot_so_columns_are_plotting_variable(
            dataframe=data, groupby=groupby)
        # Pass kwargs to heatmap cmap.
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)
        g.set_xticklabels(rotation=xrotation)

    # set some plotting parameters
    xrotation = 90
    # Calculate the size, aspect depending on the number of
    #  rows per subplot
    num_rows = len(dataframe['name_string'].unique())
    size = 1 + 0.22*num_rows
    aspect = 1.5  # a

    if cap_facet_labels:
        dataframe = capitalize_some_column_names(dataframe)
        facet_var = "Replicate"
        row_var='$O_2$'
        col_var = "Week"
    else:
        facet_var = 'rep'
        row_var = 'oxy'
        col_var = 'week'

    # todo: this doesn't seem to be changing the font size.  Probably isn't
    # for other plotting calls either!
    with sns.plotting_context(font_scale=40):
        g = sns.FacetGrid(dataframe,
                          col=facet_var,
                          row=row_var,
                          size=size,
                          aspect=aspect,
                          margin_titles=True)

    g.set_axis_labels(col_var)

    # Add axes for the colorbar.  [left, bottom, width, height]
    cbar_ax = g.fig.add_axes([.94, .3, .02, .4], title='fraction \n of reads')

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0,
                        # MUST SET VMAX or all of the subplots will be on
                        # their own color scale and you might not know it.
                        vmax=dataframe['fraction of reads'].max(),
                        annot=False,
                        groupby=col_var,
                        xrotation=90)


    # modify labels
    # Todo: make the 2nd argument a function
    y_label_formatter(g, italics_unless_other)

    # add space for x label
    g.fig.subplots_adjust(bottom=0.2)

    # room for colorbar (cbar)
    g.fig.subplots_adjust(right=0.85)

    # add a supertitle, you bet.
    supertitle_base = taxa_dict_to_descriptive_string(taxa_dict)
    if title:
        plt.subplots_adjust(top=0.80)
        supertitle = \
            supertitle_base + '.  Min fraction of reads cutoff = {}'.format(
                low_cutoff)
        g.fig.suptitle(supertitle, size=15)

    # Also summarize # of taxa rows being grouped together.

    # prepare filename and save.
    plot_dir = elviz_utils.prepare_plot_dir(plot_dir)
    filepath = plot_dir + supertitle_base
    filepath += "--min_{}".format(low_cutoff)
    filepath += "--{}".format('x-week')
    filepath += ".pdf"
    print(filepath)
    g.fig.savefig(filepath)

    if svg:
        g.fig.savefig(filepath.rstrip("pdf") + "svg")

    return g

def italic_if_not_other(s):
    if s == 'other':
        return s
    else:
        return '\\textit{' + s + '}'


def bar_facets_from_pivoted_df(not_pivoted_df, plot_x, plot_y, order_list,
                               color_list, filename=None, portrait=True):
    """
    General function to produce faceted plots (O2 vs rep or vice versa)
    given already-pivotd data.

    :param not_pivoted_df: a dataframe that will be pivoted before plotting
    :param plot_x: value to use as x in pivot (and plot)
    :param plot_y: value to use as y (bar rectangle) in pivot and plot
    :param order_list: order to plot columns by.  Should correspond to color list.
    :param color_list: list of hex (or RGB?) colors to use
    :param filename: filename to save the file with
    :param portrait: True if tall is desired, or False if short
    :return:
    """
    if portrait:
        fig, axs = plt.subplots(4, 2, figsize=(10,10))
        axd = axd_portrait(axs)
    else:
        fig, axs = plt.subplots(2, 4, figsize=(14,8))
        axd = axd_landscape(axs)

    for (o2, rep), df in not_pivoted_df.groupby(['oxy', 'rep']):
        ax = axd[(o2, rep)]
        ax.set_title(o2.lower() + ' $\mathrm{O_2}$' + ' replicate {}'.format(rep))

        plot_df = df.pivot(index='week', columns=plot_x, values=plot_y)
        plot_df = plot_df[order_list]
        plot_df.plot.bar(stacked=True, ax=ax, legend=False, color=color_list)

    if portrait:
        # prevent subplot overlaps: set width, height to leave between subplots.
        plt.subplots_adjust(wspace = 0.3, hspace = 0.6)
        # add legend to the upper right
        axd[('High', 1)].legend(loc=(1.05, 0))
    else:
        # prevent subplot overlaps: set width, height to leave between subplots.
        plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
        # add legend to the upper right
        axd[('Low', 4)].legend(loc=(1.05, 0.2))

    # put labels up the left side.
    for ax in axs[:, 0]:
        ax.set_ylabel('fractional abundance')

    return fig

def plot_dominant_methylotrophs(genera_df, filename=None, portrait=True):
    """
    Plot the key methanotroph and methylotroph abundances by genera.
    :param genera_df: Dataframe summarized to the genera level
    :param filename: filename to save picture with
    :param portrait: True if tall, False if wide.
    :return: matplotlib figure
    """
    # make a copy, for safety.
    dataframe = genera_df.copy()

    order_list = [# methanotrophs:
                  'Methylobacter', 'Methylovulum',
                  'Methylomonas', 'Methylosarcina',
                  # methylotrophs:
                  'Methylophilus','Methylotenera']
    order_list_italic = ['\\textit{' + G + '}' for G in order_list]

    # Makes all the legend labels italic.
    dataframe['Genus italics'] = \
        genera_df['Genus'].apply(italic_if_not_other)

    # specify colors, dark to light.
    colors = ['#810f7c', '#8856a7', '#8c96c6', '#9ebcda',    # purples
              '#006d2c', '#74c476'] # greens

    fig = bar_facets_from_pivoted_df(
        not_pivoted_df=dataframe,
        plot_x='Genus italics', plot_y='fraction of reads',
        order_list=order_list_italic, color_list=colors,
        portrait=portrait, filename=filename)
    fig = add_vline_to_all_subplots(fig, x=10.5-4, ymin=0, ymax=1, color='#636363')

    if filename is not None:
        fig.savefig(filename, bbox_inches='tight')

    return fig


def plot_bars_for_taxa_dict(dataframe, taxa_dict, order_list,
                            colors, filename, main_dir, portrait=True,
                            summarize_other=True, check_totals_sum_to_1=True):

    plot_data = aggregate_mixed_taxonomy(
        dataframe=dataframe,
        taxa_dict=taxa_dict,
        main_dir=main_dir,
        summarize_other=summarize_other,
        check_totals_sum_to_1=check_totals_sum_to_1)
    # todo: move this into abundance utils and out of jupyter/.py:
    x = 'taxonomic name'
    y = 'abundance sum'
    plot_data.ix[plot_data[x] == 'other', x] = r'other \& unknown'
    #print(plot_data['taxonomic name'].drop_duplicates())

    fig = bar_facets_from_pivoted_df(
        not_pivoted_df=plot_data,
        plot_x=x, plot_y=y,
        order_list=order_list, color_list=colors,
        portrait=portrait, filename=filename)

    fig = add_vline_to_all_subplots(fig, x=10.5-4, ymin=0, ymax=1, color='#636363')

    if filename is not None:
        fig.savefig(filename, bbox_inches='tight')

    return fig

