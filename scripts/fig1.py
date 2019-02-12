import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

PLOTS_DIR = '../plots'

"""
This python script plots the summary data for
the ox-experiments on Amazon Mechanical Turk.
"""

data_file = pd.read_csv('../data/ox/all_ox_untrimmed_anonymous.csv')


def plot_stats(df_all):
    true_value = 1233
    views = [1, 3, 9]
    colors = ['#fdae6b', '#f16913', '#d94801', '#7f2704']

    # global figure settings:
    fig, axarr = plt.subplots(nrows=4, sharex=True,
                              gridspec_kw={'height_ratios': [3, 3, 3, 1]},
                              figsize=(5.5, 4))
    plt.xlabel('log(estimate)', fontweight='bold')

    # make control data
    df_control = df_all[df_all['views'] == 0]

    # calculate quartiles and deciles
    quan_min = df_control['guess'].apply(np.log).quantile(0.25)
    quan_max = df_control['guess'].apply(np.log).quantile(0.75)
    dec_min = df_control['guess'].apply(np.log).quantile(0.1)
    dec_max = df_control['guess'].apply(np.log).quantile(0.9)

    # plot control data
    axarr[3].plot(np.log(np.median(df_control.guess.values)),
                  [0], 'o', markersize=9, c=colors[0])
    axarr[3].plot(np.log(np.mean(df_control.guess.values)),  # the logarithm of the mean, not the mean of the logarithms!
                  [0], 'v', markersize=7, c=colors[0])
    axarr[3].hlines(y=0, xmin=quan_min, xmax=quan_max,
                    linewidth=4, color=colors[0])
    axarr[3].hlines(y=0, xmin=dec_min, xmax=dec_max,
                    linewidth=2, color=colors[0])
    axarr[3].axvline(x=np.log(true_value),
                     linewidth=1, color='black')

    # plotting paraphernalia for the control data
    axarr[3].set_yticks([-0.5, 0, 0.5])
    axarr[3].tick_params(axis="y", length=0)  # remove small ticks
    axarr[3].set_yticklabels(
        ['', 'v='+str(0)+' (N= ' + str(len(df_control)) + ')   ', ''])

    plt.yticks([])  # removing axis numbers with empty list

    # plot the rest of the data
    treatments = ['max', 'min', 'history']
    lv = [[''] for i in range(3)]
    for treat_idx, treat in enumerate(treatments):
        df_t = df_all[df_all['method'] == treat]
        for view_idx, view in enumerate(views):
            df = df_t[df_t['views'] == view]

            # calculate the p-value using the Mann-Whitney test and
            # plot number of stars for it
            _, pvalue = mannwhitneyu(df_control.guess.values, df.guess.values)
            if 0.01 < pvalue <= 0.05:
                lv[treat_idx].append('v='+str(view)+' (N=  ' + str(len(df)) + ')*   ')
            elif 0.001 < pvalue <= 0.01:
                lv[treat_idx].append('v='+str(view)+' (N=' + str(len(df)) + ')**  ')
            elif pvalue <= 0.001:
                lv[treat_idx].append('v='+str(view)+' (N=' + str(len(df)) + ')***')
            else:
                lv[treat_idx].append('v='+str(view)+' (N=' + str(len(df)) + ')    ')

            # calculate quartiles and deciles
            quan_min = df['guess'].apply(np.log).quantile(0.25)
            quan_max = df['guess'].apply(np.log).quantile(0.75)
            dec_min = df['guess'].apply(np.log).quantile(0.1)
            dec_max = df['guess'].apply(np.log).quantile(0.9)

            # plot the thread stats for the given number of views
            axarr[treat_idx].plot(np.log(np.median(df['guess'].values)),
                                 [view_idx], 'o', markersize=9, c=colors[view_idx+1])
            axarr[treat_idx].plot(np.log(np.mean(df['guess'].values)),
                                 [view_idx], 'v', markersize=7, c=colors[view_idx+1])
            axarr[treat_idx].hlines(y=view_idx, xmin=quan_min, xmax=quan_max,
                                   linewidth=4, color=colors[view_idx+1])
            axarr[treat_idx].hlines(y=view_idx, xmin=dec_min, xmax=dec_max,
                                   linewidth=2, color=colors[view_idx+1])

        # plotting paraphernalia for the subfigure for the treatment type
        axarr[treat_idx].axvline(x=np.log(true_value), linewidth=1, color='black')
        axarr[treat_idx].set_yticks([-0.5, 0, 1, 2, 2.5])
        axarr[treat_idx].tick_params(axis="y", length=0)  # remove small ticks
        axarr[treat_idx].set_yticklabels(lv[treat_idx])

        # write treatment type on right side:
        axr = axarr[treat_idx].twinx()
        axr.set_ylabel(str(treat), fontsize='large', fontweight='bold')
        plt.yticks([])  # removing axis numbers with empty list

    # plotting paraphernalia
    # remove small ticks and make grid
    for i in range(4):
        axarr[i].tick_params(axis="x", length=0)
        axarr[i].grid(True)
    plt.tight_layout()

    # Remember: save as pdf and transparent=True for Adobe Illustrator
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)

    output_path = os.path.join(PLOTS_DIR, 'fig1.png')

    plt.savefig(output_path, transparent=True, dpi=300)
    plt.savefig(output_path, transparent=True, dpi=300)
    plt.show()


# load the data
dataframe = pd.DataFrame(data_file)

plot_stats(dataframe)
