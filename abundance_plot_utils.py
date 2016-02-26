import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import abundance_utils
import elviz_utils


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
    # get rid of oxygen levels and replicates if specified.
    if oxy is not 'all':
        print("keep only {} oxygen samples".format(oxy))
        dataframe = dataframe[dataframe['oxy'] == oxy]
    if rep is not 'all':
        print("keep only replicate levels:", rep)
        dataframe = dataframe[dataframe['rep'].isin(rep)]
    dataframe = abundance_utils.filter_by_abundance(data=dataframe,
                                                    column='abundance',
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
                                values='abundance')
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)

    with sns.plotting_context(font_scale=7):
        g = sns.FacetGrid(dataframe, col='facet_replicate',
                          margin_titles=True,
                          size=plot_size, aspect=plot_aspect)
        g.set_xticklabels(rotation=90)

    # Create a colorbar axes
    cbar_ax = g.fig.add_axes([.94, .3, .02, .4], title='abundance')

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0,
                        # Specify the colorbar axes and limits
                        vmax=max(dataframe.abundance))

    g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=18)
    g.set_axis_labels('week')

    # Add space so the colorbar doesn't overlap the plot
    g.fig.subplots_adjust(right=.9)

    # add a supertitle, you bet.
    plt.subplots_adjust(top=0.80)
    supertitle = str(low) + ' < abundance < ' + str(
        high) + ', {} oxygen'.format(oxy)
    g.fig.suptitle(supertitle, size=18)

    # write a filename and save.
    filename = oxy + "_oxygen--{0}_to_{1}_abundance".format(low, high)
    print('filename:', filename)

    plot_dir = elviz_utils.prepare_plot_dir(plot_dir)

    # save figure
    g.savefig(plot_dir + filename + '.pdf')


def subset_on_phylogeny(dataframe, phylo_level, name):
    """
    Return only rows of the datframe where the value in column phylo_level
    matches the specified name.

    :param dataframe:
    :param phylo_level:
    :param name:
    :return:
    """
    print(dataframe.columns)
    return dataframe[dataframe[phylo_level] == name]


def other_phylogeny_levels(level):
    """
    return the name of all the phylogenetic levels *not* specified by level.

    Handy if you want to drop all the other columns.

    :param level: string string string string
    :return:
    """
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    levels.remove(level)
    return levels


def phylo_levels_above(phylo_level):
    """
    E.g. 'Order' --> ['Kingdom', 'Phylum', 'Class']
    """
    p_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    position_of_phylo_level = p_levels.index(phylo_level)
    return p_levels[0:position_of_phylo_level]


def phylo_levels_below(phylo_level):
    """
    E.g. 'Order' --> ['Family', 'Genus']
    """
    p_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    position_of_phylo_level = p_levels.index(phylo_level)
    return p_levels[position_of_phylo_level + 1:]


def phyla_below_level(dataframe, phylo_dict):
    pass


def sum_on_phylogeny(dataframe, phylo_level, name):
    # E.g. if you pass phylo_level = 'Phylum' and name = 'Bacteroidetes'
    # You will get one row with columns ['phylogenetic label', 'name',
    # 'sum of abundances'] *per* sample.
    # The sum of abundances results from groupby aggregations within the
    # sample group.
    relevant_rows = subset_on_phylogeny(dataframe=dataframe,
                                        phylo_level=phylo_level,
                                        name=name)
    # Collaps all phylogeny below:
    for col_name in other_phylogeny_levels(phylo_level):
        del relevant_rows[col_name]

    # sum all
    aggregated_rows = relevant_rows.groupby(
        [phylo_level, 'ID'])['abundance'].sum()

    return aggregated_rows


def aggregate_mixed_phylogeny(dataframe, phylo_dict):
    # Loop over the different phylogenetic levels specified.
    # Make a list of each dataframe that will be concatenated.
    reduced_data = []
    for key in phylo_dict.keys():
        for name in phylo_dict[key]:
            reduced_rows = sum_on_phylogeny(dataframe=dataframe,
                                            phylo_level=key,
                                            name=name)

            # check that you got some rows.  Might be a typo if not!
            assert(reduced_rows.shape[0] > 0), \
                'found no rows for {} = "{}"'.format(key, name)

            # the index needs to be dropped but it is stored below as
            # 'phylogenetic level' and 'phylogenetic name'
            # I haven't been able to reset_index on this series to drop the
            # index so I'm doing it this way:
            reduced_rows = reduced_rows.reset_index()
            del reduced_rows[key]
            # make a new dataframe out of it.
            reduced_data.append(
                pd.DataFrame({'phylogenetic level': key,
                              'phylogenetic name': name,
                              'abundance sum': reduced_rows['abundance'],
                              'ID': reduced_rows['ID']}))
    # Concatenate data
    dataframe = pd.concat(reduced_data)
    # merge on the sample info.
    dataframe = pd.merge(left=dataframe, right=elviz_utils.read_sample_info())

    return dataframe


def phylo_dict_to_descriptive_string(phylo_dict):
    # todo: go through highest orders of phylo first.
    # e.g. phylum looped over before order.
    desc_string = ""
    for key in phylo_dict:
        print(key)
        desc_string += key
        desc_string += '-'
        for value in phylo_dict[key]:
            desc_string += value + '_'
        desc_string = desc_string[:-1]
        desc_string += '--'
    # remove last two '--' characters
    desc_string = desc_string[:-2]
    return desc_string


def plot_across_phylogeny(dataframe, phylo_dict,
                          facet='week', annotate=True,
                          plot_dir='./plots'):

    # todo: What happens if you submit a Genus for something you also
    # submitted an order for???   For now assume the user is smarter than that.
    plot_data = aggregate_mixed_phylogeny(dataframe=dataframe,
                                          phylo_dict=phylo_dict)

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
        return dataframe.pivot(index='phylogenetic name',
                               columns=groupby,
                               values='abundance sum')

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
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)
        g.set_xticklabels(rotation=xrotation)

    # todo: add a label at the bottom like "replicate" or "week"
    # currently replicate is turned into facet_replicate but should just
    # make a label that says replicate.  Week

    # Control plot aesthetics depending on facet option.
    if facet == 'week':
        xrotation = 0
        num_rows = len(plot_data['phylogenetic name'].unique())
        size = 2 * 0.2*num_rows
        aspect = 1
        space_for_cbar = 0.85
        x_axis_label = 'replicate'

    else:
        xrotation = 90
        # Calculate the size, aspect depending on the number of
        #  rows per subplot
        num_rows = len(plot_data['phylogenetic name'].unique())
        size = 0.9 + 0.2*num_rows
        aspect = 1.2
        space_for_cbar = 0.85
        x_axis_label = 'week'
    # todo: make wider if annotate = True.

    with sns.plotting_context(font_scale=7):
        g = sns.FacetGrid(plot_data,
                          col=facet,
                          row='oxy',
                          size=size,
                          aspect=aspect,
                          margin_titles=True)

    # Add axes for the colorbar.  [left, bottom, width, height]
    cbar_ax = g.fig.add_axes([.92, .3, .02, .4], title='abundance')

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0, annot=annotate,
                        groupby=cols_in_facet,
                        xrotation=xrotation)

    g.set_axis_labels(x_axis_label)

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

    # add a supertitle, you bet.
    plt.subplots_adjust(top=0.80)
    supertitle = phylo_dict_to_descriptive_string(phylo_dict)
    g.fig.suptitle(supertitle, size=15)

    # Also summarise # of taxa rows being grouped together.

    # prepare filename and save.
    plot_dir = elviz_utils.prepare_plot_dir(plot_dir)
    filepath = plot_dir + supertitle
    filepath += "--{}".format(facet)
    filepath += ".pdf"
    print(filepath)
    g.fig.savefig(filepath)

    return g


def label_from_phylo_colnames(*args):
    """
    e.g. ['Burkholderiales', 'Comamonadaceae, 'other'] -->
        'Burkholderiales_Comamonadaceae_other', or
    ['Burkholderiales', NaN, 'other']
    """
    #print('args: {}'.format(args))
    name_string = ""
    for name in args:
        if name != 'other':
            if name != 'unknown':
                #if not math.isnan(name):   #np.isnan(name):
                #print('adding name {}'.format(name))
                name_string += name
                name_string += ", "
    # remove the last ", "
    if len(name_string) > 0:
        #print('length of {} is > 0'.format(name_string))
        return name_string[:-2]
    else:
        #print('all fields empty.  returning "?"')
        return '?'


def heatmap_all_below(dataframe, phylo_dict, plot_dir):
    # grab the data for that phylo:
    # for now assume jusst 1 key and 1 value.
    phylo_level = list(phylo_dict.keys())[0]
    phylo_name = list(phylo_dict.values())[0][0]
    dataframe = dataframe[dataframe[phylo_level] == phylo_name]
    print(dataframe.head())

    # Columns to form a concatenated label from:
    label_cols = phylo_levels_below(phylo_level=phylo_level)
    print('label_cols: {}'.format(label_cols))

    # change nan cells to 'unknown'
    dataframe.fillna('unknown', inplace=True)


    # make a summary string representing the phylogeny for everything below

    def label_building_lambda(f, columns):
        """
        Returns a lambda function to make row labels from.
        :param f: function to make a lambda out of.
        :param columns: column names to pass to function f in the lambda
        :return: function
        """
        return lambda row: f(*(row[col] for col in columns))

    # TODO: use the phylo_dict to get the columns to use!
    dataframe['name_string'] = dataframe.apply(
        label_building_lambda(f=label_from_phylo_colnames,
                              columns=label_cols), axis=1)

    # Plot as usual, using the stuff developed above.
    # todo: factor some of this??
    def pivot_so_columns_are_plotting_variable(dataframe, groupby):
        return dataframe.pivot(index='name_string',
                               columns=groupby,
                               values='abundance')

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
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)
        g.set_xticklabels(rotation=xrotation)

    with sns.plotting_context(font_scale=10):
        g = sns.FacetGrid(dataframe,
                          col='rep',
                          row='oxy',
                          size=10,
                          aspect=.5,
                          margin_titles=True)

    # Add axes for the colorbar.  [left, bottom, width, height]
    cbar_ax = g.fig.add_axes([.94, .3, .02, .4], title='abundance')

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0, annot=False,
                        groupby='week',
                        xrotation=0)

    g.set_axis_labels('week')

    # add space for x label
    g.fig.subplots_adjust(bottom=0.2)

    # room for colorbar (cbar)
    g.fig.subplots_adjust(right=0.85)

    # add a supertitle, you bet.
    plt.subplots_adjust(top=0.93)
    supertitle = phylo_dict_to_descriptive_string(phylo_dict)
    g.fig.suptitle(supertitle, size=15)

    # Also summarise # of taxa rows being grouped together.

    # prepare filename and save.
    plot_dir = elviz_utils.prepare_plot_dir(plot_dir)
    filepath = plot_dir + supertitle
    filepath += "--{}".format('x-week')
    filepath += ".pdf"
    print(filepath)
    g.fig.savefig(filepath)

    return g

