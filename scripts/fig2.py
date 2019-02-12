import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

PLOTS_DIR = '../plots'

"""
This python script plotting summary statistics for
the dots-experiments on Amazon Mechanical Turk.
"""

data_file = pd.read_csv('../data/dots/all_dots_trimmed_anonymous.csv')


def plot_stats(df_all):
    dots = [55, 148, 403, 1097]
    views = [0, 1, 3, 9]
    colors = ['#fdae6b', '#f16913', '#d94801', '#7f2704']

    # global figure settings:
    fig, axarr = plt.subplots(nrows=4, sharex=True, figsize=(5.5, 4))
    plt.xlabel('log(estimate)', fontweight='bold')

    # plot
    lv = [[''] for i in range(5)]
    for true_dots_idx, true_dots in enumerate(dots):
        df_i = df_all[df_all['dots'] == true_dots]

        # get control data
        df_control = df_i[df_i['views'] == 0]
        for view_idx, view in enumerate(views):
            df = df_i[df_i['views'] == view]

            # if view == 0:
            #     control = df.guess.values

            # calculate the p-value using the two-sided Mann-Whitney U test
            _, pvalue = mannwhitneyu(df_control.guess.values, df.guess.values)
            if 0.01 < pvalue <= 0.05:
                lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(df))+')*   ')
            elif 0.001 < pvalue <= 0.01:
                lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(df))+')** ')
            elif pvalue <= 0.001:
                lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(df))+')***')
            else:
                lv[true_dots_idx].append('v='+str(view)+' (N='+str(len(df))+')    ')

            # calculate quartiles and deciles
            quan_min = df['guess'].apply(np.log).quantile(0.25)
            quan_max = df['guess'].apply(np.log).quantile(0.75)
            dec_min = df['guess'].apply(np.log).quantile(0.1)
            dec_max = df['guess'].apply(np.log).quantile(0.9)

            # plot the thread stats for the given number of views
            axarr[true_dots_idx].plot(np.log(np.median(df['guess'].values)),
                                 [view_idx], 'o', markersize=9, c=colors[view_idx])
            axarr[true_dots_idx].plot(np.log(np.mean(df['guess'].values)),
                                 [view_idx], 'v', markersize=7, c=colors[view_idx])
            axarr[true_dots_idx].hlines(y=view_idx, xmin=quan_min, xmax=quan_max,
                                   linewidth=4, color=colors[view_idx])
            axarr[true_dots_idx].hlines(y=view_idx, xmin=dec_min, xmax=dec_max,
                                   linewidth=2, color=colors[view_idx])

        axarr[true_dots_idx].axvline(x=np.log(true_dots),
                               linewidth=1, color='black')

        # plotting paraphernalia for the subfigure for the number of dots
        axarr[true_dots_idx].set_yticks([-0.5, 0, 1, 2, 3, 3.5])
        axarr[true_dots_idx].tick_params(axis="y", length=0)  # remove small ticks
        axarr[true_dots_idx].set_yticklabels(lv[true_dots_idx])
        axr = axarr[true_dots_idx].twinx()
        axr.set_ylabel('d='+str(true_dots), fontweight='bold')
        plt.yticks([])  # removing axis numbers with empty list

    for i in range(4):
        axarr[i].tick_params(axis="x", length=0)
        axarr[i].grid(True)

    axarr[0].set_xlim([3.7, 7.7])  # set x-axis dimensions
    plt.tight_layout()

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    output_path = os.path.join(PLOTS_DIR, 'fig2.png')

    plt.savefig(output_path, transparent=True, dpi=300)
    plt.savefig(output_path, transparent=True, dpi=300)
    plt.show()


# load the data
dataframe = pd.DataFrame(data_file)

# only plot the data for the history sessions
dataframe = dataframe[dataframe.method == 'history']
plot_stats(dataframe)
