import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from abundance_utils import filter_by_abundance


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
