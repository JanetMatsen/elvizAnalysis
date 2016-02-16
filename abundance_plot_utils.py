import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from abundance_utils import filter_by_abundance
from abundance_utils import normalize_groupby
from elviz_utils import read_sample_info
from abundance_utils import make_directory


def plot_heatmap(dataframe, high, low, oxy, rep, plot_dir):
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
    dataframe = filter_by_abundance(data=dataframe,
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
        g.set_xticklabels(rotation=30)

    # Create a colorbar axes
    # TODO: add label
    cbar_ax = g.fig.add_axes([.92, .3, .02, .4])

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
            print(name)
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
            print('reduced_rows.head(2)')
            print(reduced_rows.head(2))
            del reduced_rows[key]
            print('reduced_rows.head(2)')
            print(reduced_rows.head(2))
            # make a new dataframe out of it.
            reduced_data.append(
                pd.DataFrame({'phylogenetic level': key,
                              'phylogenetic name': name,
                              'abundance sum': reduced_rows['abundance'],
                              'ID': reduced_rows['ID']}))
            print(reduced_data[-1].head(2))
    # Concatenate data
    dataframe = pd.concat(reduced_data)
    # merge on the sample info.
    print('dataframe.head()')
    print(dataframe.head())
    dataframe = pd.merge(left=dataframe, right=read_sample_info())

    return dataframe


def phylo_dict_to_filename(phylo_dict):
    filename = ""
    for key in phylo_dict:
        print(key)
        filename += key
        filename += '-'
        for value in phylo_dict[key]:
            filename += value + '_'
        filename += '--'
    # remove last two '--' characters
    filename = filename[:-3]
    filename += ".pdf"
    return filename


def plot_across_phylogeny(dataframe, phylo_dict, facet='week', annotate=True):

    # What happens if you submit a Genus for something you also submitted an
    # order for ???   For now assume the user is smarter than that.
    plot_data = aggregate_mixed_phylogeny(dataframe=dataframe,
                                          phylo_dict=phylo_dict)
    plot_data['facet_replicate'] = 'replicate ' + plot_data['rep'].astype(str)

    # The data is seperated by these two variables.
    # The one not used as the facet will be used as the columns in the
    # subplot.
    if facet == 'week':
        cols_in_facet = 'facet_replicate'
    else:
        cols_in_facet = 'week'

    print('plot_data.head()')
    print(plot_data.head())

    def pivot_so_columns_are_plotting_variable(dataframe, groupby):
        return dataframe.pivot(index='phylogenetic name',
                               columns=groupby,
                               values='abundance sum')

    def facet_heatmap(data, groupby, **kws):
        """
        Used to fill the subplots with data.

        :param data: dataframe to plot
        :param kws:
        :return:
        """
        # pivot only supports one column for now.
        # http://stackoverflow.com/questions/32805267/pandas-pivot-on-multiple-columns-gives-the-truth-value-of-a-dataframe-is-ambigu
        facet_data = pivot_so_columns_are_plotting_variable(
            dataframe=data, groupby=groupby)
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)
        g.set_xticklabels(rotation=30)

    # Calculate the size, aspect depending on the number of rows per subplot
    num_rows = len(plot_data['phylogenetic name'].unique())
    size = 1.5 + 0.2*num_rows
    aspect = 1

    with sns.plotting_context(font_scale=7):
        g = sns.FacetGrid(plot_data,
                          col=facet,
                          row='oxy',
                          size=size,
                          aspect=aspect,
                          margin_titles=True)

    # TODO: add label for color bar.
    cbar_ax = g.fig.add_axes([.92, .3, .02, .4])

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0, annot=annotate,
                        groupby=cols_in_facet)

    # Add space so the colorbar doesn't overlap th plot.
    g.fig.subplots_adjust(right=0.9)

    # add a supertitle, you bet.
    plt.subplots_adjust(top=0.85)
    supertitle = phylo_dict_to_filename(phylo_dict)
    g.fig.suptitle(supertitle, size=18)

    # Also summarise # of taxa rows being grouped together.

    # prepare filename and save.
    plotdir = './plots/mixed_phylogeny/'
    make_directory(plotdir)
    filepath = plotdir + supertitle
    print(filepath)
    g.fig.savefig(filepath)

    return g
