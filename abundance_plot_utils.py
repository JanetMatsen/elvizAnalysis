import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
import pandas as pd
import seaborn as sns

from abundance_utils import filter_by_abundance
from abundance_utils import normalize_groupby
from elviz_utils import read_sample_info

def plot_heatmap(data, high, low, oxy, rep, plot_dir):
    """
    Make a heatmap at Genus, using oganisms withing the specified abundance
    cutoffs.

    :param data: dataframe to pass
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
        data = data[data['oxy'] == oxy]
    if rep is not 'all':
        print("keep only replicate levels:", rep)
        data = data[data['rep'].isin(rep)]
    data = filter_by_abundance(data=data,
                               column='abundance',
                               high=high, low=low)
    data['facet_replicate'] = 'replicate ' + data['rep'].astype(str)

    # make height of the plot a function of the number of rows (Genera):
    num_data_rows = len(data['Genus'].unique())
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

        :param facet_data:
        :param kws:
        :return:
        """

        facet_data = data.pivot(index='Genus', columns='week',
                                values='abundance')
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)

    with sns.plotting_context(font_scale=7):
        g = sns.FacetGrid(data, col='facet_replicate',
                          margin_titles=True,
                          size=plot_size, aspect=plot_aspect)

        # Create a colorbar axes
    cbar_ax = g.fig.add_axes([.92, .3, .02, .4])

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0,
                        # Specify the colorbar axes and limits
                        vmax=max(data.abundance))

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
    aggregated_rows = relevant_rows.groupby([phylo_level, 'ID'])['abundance'].sum()

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
    return pd.concat(reduced_data)


def plot_across_phylogeny(dataframe, phylo_dict):

    # What happens if you submit a Genus for something you also submitted an
    # order for ???   For now assume the user is smarter than that.
    abundances = aggregate_mixed_phylogeny(dataframe=dataframe,
                                           phylo_dict=phylo_dict)

    # merge on sample_info using ID column.
    sample_info = read_sample_info()
    plot_data = pd.merge(abundances.reset_index(), sample_info)

    print('plot_data.head()')
    print(plot_data.head())

    def facet_heatmap(data, **kws):
        """
        Used to fill the subplots with data.

        :param facet_data:
        :param kws:
        :return:
        """

        facet_data = data.pivot(index='Genus', columns='week',
                                values='abundance')
        # Pass kwargs to heatmap  cmap used to be 'Blue'
        sns.heatmap(facet_data, cmap="YlGnBu", **kws)

    with sns.plotting_context(font_scale=7):
        g = sns.FacetGrid(plot_data,
                          col='rep',
                          #col='facet_replicate',
                          margin_titles=True,
                          #size=plot_size,
                          #aspect=plot_aspect
                        )

    cbar_ax = g.fig.add_axes([.92, .3, .02, .4])
    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax, vmin=0,
                        # Specify the colorbar axes and limits
                        #vmax=max(data.abundance)
                        )

    # Also summarise # of taxa rows being grouped together.

    return g

