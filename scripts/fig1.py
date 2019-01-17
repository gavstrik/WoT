import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import style, rc
import matplotlib
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.family"] = "sans-serif"

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
    df_c = df_all[df_all['views'] == 0]

    # calculate quartiles and deciles
    quan_min = df_c['guess'].apply(np.log).quantile(0.25)
    quan_max = df_c['guess'].apply(np.log).quantile(0.75)
    dec_min = df_c['guess'].apply(np.log).quantile(0.1)
    dec_max = df_c['guess'].apply(np.log).quantile(0.9)

    # plot control data
    axarr[3].plot(np.log(np.median(df_c.guess.values)),
                  [0], 'o', markersize=9, c=colors[0])
    axarr[3].plot(np.log(np.mean(df_c.guess.values)),  # the logarithm of the mean, not the mean of the logarithms!
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
        ['', 'v='+str(0)+' (N= ' + str(len(df_c)) + ')   ', ''])
    axar = axarr[3].twinx()
    # axar.set_ylabel('control', fontsize='large', fontweight='bold')
    plt.yticks([])  # removing axis numbers with empty list

    # plot the rest of the data
    treatments = ['max', 'min', 'history']
    lv = [[''] for i in range(3)]
    for position, treat in enumerate(treatments):
        df_t = df_all[df_all['method'] == treat]
        for pos, view in enumerate(views):
            df = df_t[df_t['views'] == view]
            sessions = df.session.unique()

            # calculate the p-value using the Mann-Whitney test
            _, pvalue = mannwhitneyu(df_c.guess.values, df.guess.values)

            # plot stars according to p-value
            if pvalue <= 0.05 and pvalue > 0.01:
                lv[position].append('v='+str(view)+' (N=  ' + str(len(df)) + ')*   ')
            elif pvalue <= 0.01 and pvalue > 0.001:
                lv[position].append('v='+str(view)+' (N=' + str(len(df)) + ')**  ')
            elif pvalue <= 0.001:
                lv[position].append('v='+str(view)+' (N=' + str(len(df)) + ')***')
            else:
                lv[position].append('v='+str(view)+' (N=' + str(len(df)) + ')    ')

            # calculate quartiles and deciles
            quan_min = df['guess'].apply(np.log).quantile(0.25)
            quan_max = df['guess'].apply(np.log).quantile(0.75)
            dec_min = df['guess'].apply(np.log).quantile(0.1)
            dec_max = df['guess'].apply(np.log).quantile(0.9)

            # plot
            axarr[position].plot(np.log(np.median(df['guess'].values)),
                                 [pos], 'o', markersize=9, c=colors[pos+1])
            axarr[position].plot(np.log(np.mean(df['guess'].values)),
                                 [pos], 'v', markersize=7, c=colors[pos+1])
            axarr[position].hlines(y=pos, xmin=quan_min, xmax=quan_max,
                                   linewidth=4, color=colors[pos+1])
            axarr[position].hlines(y=pos, xmin=dec_min, xmax=dec_max,
                                   linewidth=2, color=colors[pos+1])
        axarr[position].axvline(x=np.log(true_value),
                                linewidth=1, color='black')
        axarr[position].set_yticks([-0.5, 0, 1, 2, 2.5])
        axarr[position].tick_params(axis="y", length=0)  # remove small ticks
        axarr[position].set_yticklabels(lv[position])

        # write treatment type on right side:
        axr = axarr[position].twinx()
        axr.set_ylabel(str(treat), fontsize='large', fontweight='bold')
        plt.yticks([])  # removing axis numbers with empty list


    # plotting paraphernalia
    # remove small ticks and make grid
    for i in range(4):
        axarr[i].tick_params(axis="x", length=0)
        axarr[i].grid(True)
    plt.tight_layout()
    # Remember: save as pdf and transparent=True for Adobe Illustrator
    plt.savefig('../plots/Fig1.png', transparent=True, dpi=300)
    plt.show()


# main code
df_all = pd.DataFrame(data_file)
plot_stats(df_all)
