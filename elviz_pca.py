import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib

matplotlib.use('TkAgg')
from elviz_utils import read_sample_info


def import_elviz_data():
    """
    Read pandas data from summarised table
    :return: pandas dataframe
    """
    return pd.read_csv("./results/reduced_data--genus_only.csv",
                       dtype={'abundance': 'float'})


def pivot_for_pca(dataframe):
    """
    pivot to the format required for scikitlearn pca
    :param dataframe:input dataframe reduced to genus only
    :return:dataframe with columns as samples and rows as genera
    """
    dataframe = dataframe.pivot(index='Genus', columns='ID',
                                values='abundance')
    # fill NA values with 0.
    return dataframe.fillna(0)


def most_abundant_genera_for_pca(data, top_percent):
    """
    Return genera rows with highest sum of abundances.

    Depreciated after variance was used to identify important genera
    :param data:
    :param top_percent:
    :return:
    """
    data['row_sum'] = data.sum(axis=1)
    # reduce to top __ %
    # find the row_sum corresponding to the top 10%.
    num_rows_to_keep = int(round(data.shape[0] * top_percent / 100.))
    print('number of rows to keep:', num_rows_to_keep)
    # print(data.head())
    # sort by row_sum so I can take the firs number of rows.
    data.sort(columns='row_sum', ascending=False, inplace=True)
    # print data.head()

    # data.groupby('row_sum').head(num_rows_to_keep)
    del data['row_sum']
    return data.head(num_rows_to_keep)


def most_variant_genera_for_pca(data, top_percent):
    # reduce to top __ %
    # find the row_sum corresponding to the top 10%.
    num_rows_to_keep = int(round(data.shape[0] * top_percent / 100.))
    print('number of rows to keep:', num_rows_to_keep)
    # print(data.head())
    # sort by row_sum so I can take the firs number of rows.
    data.sort(columns='variance', ascending=False, inplace=True)
    return data.head(num_rows_to_keep)


def colnames_to_sample_info_array(dataframe):
    col_df = pd.DataFrame({'ID':dataframe.reset_index().ID})
    sample_info = read_sample_info()
    return pd.merge(col_df, sample_info, how='left')


def sort_by_variance():
    df = import_elviz_data()
    df = pivot_for_pca(df)
    df['variance'] = df.var(axis=1)
    df.sort_values(by='variance',ascending=False)
    return df


def plot_variance(df=sort_by_variance()):
    fig, ax = plt.subplots()
    df.reset_index().variance.plot(ax=ax, kind='hist',
                                   bins=100, range=(0,0.05))
    ax.set_yscale('log')
    plt.title('distribution of variances (log scale)')
    plt.xlabel("variance for genera's abundance across all samples")
    plt.ylabel("frequency (log scale)")
    plt.savefig('distribution_of_sample-wise_variances.pdf')


def run_PCA(top_percent = 20):
    # get data
    df = sort_by_variance()

    # drop rows that don't have sum of abundances in the top %.
    pca_input = most_variant_genera_for_pca(data=df, top_percent=top_percent)
    print(pca_input.shape)

    # remove the variance column if it is there:
    if 'variance' in df.columns:
        del df['variance']

    # transpose so the genera are columns and samples are rows.
    pca_input = pca_input.transpose()

    # run scikitlearn PCS using default settings.
    pca = PCA()
    pca.fit(pca_input)

    # show how much each component contributed to variance
    print("principle components' contribution to variance:")
    variances =pca.explained_variance_ratio_
    print(variances)

    return pca_input, pca.fit(pca_input).transform(pca_input), variances

def plot_PCA_results(top_percent=20):
    # get data that's ready for PCA:
    pca_input, X_r, variances = run_PCA(top_percent=top_percent)

    # import the sample info needed to map features to sample IDs.
    plot_info = colnames_to_sample_info_array(pca_input)

    # prepare axis labels, which also serve as dataframe column names.
    x_axis_label = 'principle component 1 ({0:.0%})'.format(variances[0])
    y_axis_label = 'principle component 2 ({0:.0%})'.format(variances[1])

    # put together the transformed data and sample descriptions
    plot_data = pd.concat([pd.DataFrame({x_axis_label:X_r[:,0],
                                     y_axis_label:X_r[:,1]}),
                       plot_info], axis=1)

    # define a custom color palette using:
    # Conditions were seized at week ten, so seven early samples in
    # the original condition and four latest samples in an alternative
    # condition.
    color_palette = build_color_palette(num_items=14-4+1,
                                        weeks_before_switch=7)
    #color_palette = sns.cubehelix_palette(11, start=.5, rot=-.75)

    # update matplotlib params for bigger fonts, ticks:
    mpl.rcParams.update({
    'font.size': 16, 'axes.titlesize': 17, 'axes.labelsize': 15,
    'xtick.labelsize': 10, 'ytick.labelsize': 13,
    'font.weight': 600,
    'axes.labelweight': 600, 'axes.titleweight': 600})

    # Plot with Seaborn
    plt.figure(figsize=(12, 6))
    sns.set(style="ticks")
    g = sns.FacetGrid(plot_data, col="oxy", hue='week', palette=color_palette,
                     size=5, aspect=1)
    g = (g.map(plt.scatter, x_axis_label, y_axis_label,
               edgecolor="w", s=60).add_legend())
    g.fig.savefig('pca.pdf')


def build_color_palette(num_items, weeks_before_switch):
    num_pre_switch_colors = weeks_before_switch
    num_post_switch_colors = num_items - num_pre_switch_colors
    print('preparing colors for {} pre-oxygen-switch',
          'samples and {} post-switch samples' \
            .format(num_pre_switch_colors, num_post_switch_colors))

    # get the first colors from this pallete:
    pre_switch_colors = \
        sns.cubehelix_palette(11, start=.5, rot=-.75)[0:num_pre_switch_colors]
    print(pre_switch_colors)

    # get post-switch colors here:
    #post_switch_colors = sns.diverging_palette(220, 20, n=6)[::-1][0:num_post_switch_colors]
    post_switch_colors = \
        sns.color_palette("coolwarm", num_post_switch_colors)
    #sns.light_palette("navy", reverse=True)[0:num_post_switch_colors]
    rgb_colors = pre_switch_colors + post_switch_colors
    sns.palplot(rgb_colors)

    # check that we got the right amount
    print(num_items)
    assert(num_items == len(rgb_colors))
    print("")
    return rgb_colors


